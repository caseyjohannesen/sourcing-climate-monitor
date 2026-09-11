"""
Pull the latest daily SST anomaly from NOAA's OISST v2 High-Res dataset
(PSL/NOAA, 0.25 deg, daily since 1981) and write it in a form the
climate-risk-monitor globe can fetch directly.

Data source: https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.highres.html
Daily SST:          .../Datasets/noaa.oisst.v2.highres/sst.day.mean.<year>.nc
Daily climatology:  .../Datasets/noaa.oisst.v2.highres/sst.day.mean.ltm.1991-2020.nc

Both live in the *highres* directory, so the climatology sits on exactly the
same 0.25 deg grid as the daily field and the two subtract directly -- no
regridding, and no scipy dependency. (An earlier version of this script used
the 1-degree `noaa.oisst.v2.derived` climatology and interpolated it up, which
left 1-degree blocky steps in the anomaly wherever the climatological gradient
is steep: the Gulf Stream, Kuroshio and Agulhas. At 0.25 deg output that
artifact is plainly visible, so the highres baseline is not optional.)

The daily file is read via OPeNDAP (the "dodsC" THREDDS path), so xarray only
pulls the one time slice it needs over the network rather than downloading the
whole multi-gigabyte yearly file.

Output (default, --format png):
    sst_anomaly.png        1440x720 8-bit greyscale, one pixel per grid cell,
                           row 0 = north. Pixel value is a *quantization level*,
                           not a temperature; level 0 means no data (land/ice gap).
    sst_anomaly.meta.json  ~4 KB sidecar: date, source, grid geometry, and the
                           256-entry lookup table that turns a level back into
                           degrees C. The client reads the table rather than
                           hardcoding the scheme, so the quantization can change
                           here without touching any JavaScript.

Why a PNG rather than a bigger JSON: at native resolution the grid is 1,036,800
cells. As JSON that is ~6.3 MB (~0.84 MB gzipped); as a quantized PNG it is
~0.25 MB, and the browser decodes it natively into a Float32Array instead of
parsing a million boxed JS numbers. `--format json` still emits the old shape
for local inspection.

Usage:
    pip install xarray netCDF4 pandas numpy pillow requests
    python extract_sst_anomaly.py [--step 0.25] [--format png|json]

Note: this has been tested against the real NOAA grids via NOAA's CoastWatch
ERDDAP mirror, but PSL's own THREDDS host is not reachable from the environment
this was written in. It is worth confirming the first live run succeeds before
scheduling it unattended.
"""

import argparse
import calendar
import json
import os
import sys
import time
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import requests
import xarray as xr

SST_URL_TMPL = "https://psl.noaa.gov/thredds/dodsC/Datasets/noaa.oisst.v2.highres/sst.day.mean.{year}.nc"

# The daily climatology (1991-2020 baseline) never changes, so it's downloaded once as a
# plain file and cached locally rather than re-queried from NOAA's THREDDS/OPeNDAP service
# on every run -- that repeated-query pattern is what triggers the 429 rate-limit below.
CLIM_DOWNLOAD_URL = "https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2.highres/sst.day.mean.ltm.1991-2020.nc"
CLIM_CACHE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sst.day.mean.ltm.1991-2020.nc")

# ---------------------------------------------------------------------------
# Quantization
# ---------------------------------------------------------------------------
# One byte per cell, so the scheme has 255 usable levels (0 is reserved for
# "no data") to cover the whole anomaly range. A flat ramp cannot do both jobs:
# real fields run from about -7 C to +16 C (the warm extreme is a sea-ice-edge
# cell up near 68 N), so a 0.05 C ramp would clip at +/-6.35 and be wrong by
# 9 C in the Arctic, while a ramp wide enough for the tails would be too coarse
# in the tropics where the whole signal of interest lives.
#
# So the scale is piecewise: fine (0.05 C) across |anomaly| <= 5, coarse
# (0.5 C) in the tails out to -18 / +19. Worst-case error is 0.49 C and only in
# the far tails; RMS error over a real global field is 0.024 C. The display
# colour bands are 0.2-1.0 C wide, so the fine step is well below anything the
# map can show, and region averages over thousands of cells average it away.
NODATA_LEVEL = 0
CORE_MIN, CORE_MAX = -5.0, 5.0
CORE_STEP = 0.05
CORE_LEVEL_MIN = 27                                  # level for CORE_MIN
CORE_LEVEL_MAX = CORE_LEVEL_MIN + int(round((CORE_MAX - CORE_MIN) / CORE_STEP))   # 227
TAIL_STEP = 0.5
NEG_TAIL_LEVELS = CORE_LEVEL_MIN - 1                 # 26 levels: 1..26
POS_TAIL_LEVELS = 255 - CORE_LEVEL_MAX               # 28 levels: 228..255


def level_to_celsius(level):
    """Decode one quantization level back to degrees C (None for no data)."""
    if level == NODATA_LEVEL:
        return None
    if level < CORE_LEVEL_MIN:
        return round(CORE_MIN - (CORE_LEVEL_MIN - level) * TAIL_STEP, 4)
    if level > CORE_LEVEL_MAX:
        return round(CORE_MAX + (level - CORE_LEVEL_MAX) * TAIL_STEP, 4)
    return round(CORE_MIN + (level - CORE_LEVEL_MIN) * CORE_STEP, 4)


def build_lut():
    """The full 256-entry decode table, shipped in the sidecar so the client
    never has to know the scheme."""
    return [level_to_celsius(i) for i in range(256)]


def quantize(anomaly):
    """Encode a float anomaly array to uint8 levels. NaN -> NODATA_LEVEL.

    Returns (levels, clipped_low, clipped_high) so the caller can record in the
    sidecar whether any real value fell outside the representable range.
    """
    valid = np.isfinite(anomaly)
    v = np.where(valid, anomaly, 0.0)

    # core: nearest 0.05 C step
    core = np.round((v - CORE_MIN) / CORE_STEP) + CORE_LEVEL_MIN

    # tails: nearest 0.5 C step, at least one step out so a value past the core
    # edge never rounds back into the core and reads as exactly +/-5.
    neg_n = np.clip(np.round((CORE_MIN - v) / TAIL_STEP), 1, NEG_TAIL_LEVELS)
    pos_n = np.clip(np.round((v - CORE_MAX) / TAIL_STEP), 1, POS_TAIL_LEVELS)

    levels = np.where(v < CORE_MIN, CORE_LEVEL_MIN - neg_n,
             np.where(v > CORE_MAX, CORE_LEVEL_MAX + pos_n, core))
    levels = np.where(valid, np.clip(levels, 1, 255), NODATA_LEVEL)

    lo_edge = CORE_MIN - NEG_TAIL_LEVELS * TAIL_STEP
    hi_edge = CORE_MAX + POS_TAIL_LEVELS * TAIL_STEP
    clipped_low = int(np.sum(valid & (anomaly < lo_edge)))
    clipped_high = int(np.sum(valid & (anomaly > hi_edge)))
    return levels.astype(np.uint8), clipped_low, clipped_high


# ---------------------------------------------------------------------------
# Fetching
# ---------------------------------------------------------------------------
def retry_with_backoff(fn, attempts=4, base_delay=5, what=""):
    """Call fn() with exponential backoff on failure. Public government data
    servers rate-limit; this is meant to ride out a transient 429/5xx rather
    than fail the whole run."""
    last_err = None
    for attempt in range(1, attempts + 1):
        try:
            return fn()
        except Exception as e:
            last_err = e
            if attempt == attempts:
                break
            delay = base_delay * (2 ** (attempt - 1))
            print(f"  {what} attempt {attempt}/{attempts} failed ({e}); retrying in {delay}s...", file=sys.stderr)
            time.sleep(delay)
    raise RuntimeError(f"{what} failed after {attempts} attempts: {last_err}")


def download_with_retries(url, dest_path, attempts=5, base_delay=10, chunk_size=1024 * 1024):
    """Plain HTTPS download with retry/backoff, respecting Retry-After on 429s."""
    tmp_path = dest_path + ".part"
    last_err = None
    for attempt in range(1, attempts + 1):
        try:
            with requests.get(url, stream=True, timeout=60) as resp:
                if resp.status_code == 429:
                    retry_after = resp.headers.get("Retry-After")
                    delay = int(retry_after) if retry_after and retry_after.isdigit() else base_delay * (2 ** (attempt - 1))
                    print(f"  download got 429, waiting {delay}s (attempt {attempt}/{attempts})...", file=sys.stderr)
                    time.sleep(delay)
                    continue
                resp.raise_for_status()
                with open(tmp_path, "wb") as f:
                    for chunk in resp.iter_content(chunk_size=chunk_size):
                        f.write(chunk)
            os.replace(tmp_path, dest_path)
            return dest_path
        except Exception as e:
            last_err = e
            if attempt < attempts:
                delay = base_delay * (2 ** (attempt - 1))
                print(f"  download attempt {attempt}/{attempts} failed ({e}); retrying in {delay}s...", file=sys.stderr)
                time.sleep(delay)
    if os.path.exists(tmp_path):
        os.remove(tmp_path)
    raise RuntimeError(f"Failed to download {url} after {attempts} attempts: {last_err}")


def ensure_local_climatology(cache_path=CLIM_CACHE_PATH, force=False):
    """Return a local path to the climatology file, downloading it once if needed."""
    if os.path.exists(cache_path) and not force:
        print(f"Using cached climatology file: {cache_path}")
        return cache_path
    print(f"Downloading climatology file (one-time, ~350MB): {CLIM_DOWNLOAD_URL}")
    return download_with_retries(CLIM_DOWNLOAD_URL, cache_path)


def open_latest_sst():
    """Open the current year's file and fall back to last year's near Jan 1,
    when the new year's file may not have data yet."""
    now = datetime.now(timezone.utc)
    for year in (now.year, now.year - 1):
        url = SST_URL_TMPL.format(year=year)
        try:
            ds = retry_with_backoff(lambda: xr.open_dataset(url), what=f"open {year} SST file")
            if ds.sizes.get("time", 0) > 0:
                return ds
        except Exception as e:
            print(f"  could not open {url}: {e}", file=sys.stderr)
    raise RuntimeError("Could not open a current SST dataset for this year or last.")


def open_latest_sst_from_path(path):
    """Test/offline hook: same shape as open_latest_sst() but from a local file."""
    return xr.open_dataset(path)


# ---------------------------------------------------------------------------
# Anomaly
# ---------------------------------------------------------------------------
def ltm_index(target_date, n_steps):
    """Position of target_date's day in the climatology's synthetic time axis.

    The LTM is indexed by a placeholder year in calendar order, so we match by
    position rather than by date. PSL's daily LTM carries 365 steps -- it has no
    Feb 29 entry -- so from Mar 1 of a leap year onward, day-of-year runs one
    ahead of the LTM's index and has to be pulled back. (Feb 29 itself maps to
    Feb 28, the closest baseline there is.)
    """
    doy = target_date.dayofyear
    if n_steps < 366 and calendar.isleap(target_date.year) and doy >= 60:
        doy -= 1
    return min(doy - 1, n_steps - 1)


def match_climatology(clim_ds, target_date):
    """Select the climatology slice matching target_date's day-of-year."""
    return clim_ds.sst.isel(time=ltm_index(target_date, clim_ds.sizes["time"]))


def align_climatology(clim_slice, latest):
    """Put the climatology on the daily field's grid.

    With the highres baseline the two grids are already identical and this is a
    no-op, which is the whole point. The reindex is a safety net in case the
    wrong climatology file is ever cached: it keeps the run working (with a
    loud warning) and, unlike xarray's .interp(), needs no scipy.
    """
    same = (clim_slice.sizes.get("lat") == latest.sizes["lat"]
            and clim_slice.sizes.get("lon") == latest.sizes["lon"]
            and np.allclose(clim_slice.lat.values, latest.lat.values)
            and np.allclose(clim_slice.lon.values, latest.lon.values))
    if same:
        return clim_slice
    print(f"  WARNING: climatology grid {clim_slice.sizes.get('lat')}x{clim_slice.sizes.get('lon')} "
          f"does not match the daily grid {latest.sizes['lat']}x{latest.sizes['lon']}. "
          f"Falling back to nearest-neighbour reindex -- expect blocky anomalies. "
          f"Re-run with --refresh-climatology to pull the highres baseline.", file=sys.stderr)
    return clim_slice.reindex(lat=latest.lat, lon=latest.lon, method="nearest")


def compute_anomaly(sst_ds, clim_ds, step_deg):
    """Return (anomaly_2d_north_up, lats_desc, lons_asc, target_date).

    Rows come back north-first (row 0 = +89.875) to match image convention, and
    columns in ascending -180..180 longitude.
    """
    latest = sst_ds.sst.isel(time=-1)
    target_date = pd.Timestamp(sst_ds.time.isel(time=-1).values)

    clim_slice = align_climatology(match_climatology(clim_ds, target_date), latest)
    anomaly = latest - clim_slice

    lat_vals = latest.lat.values
    lon_vals = latest.lon.values
    native = abs(float(lat_vals[1] - lat_vals[0]))
    lat_step = max(1, int(round(step_deg / native)))
    lon_step = max(1, int(round(step_deg / abs(float(lon_vals[1] - lon_vals[0])))))
    if lat_step > 1 or lon_step > 1:
        anomaly = anomaly.isel(lat=slice(None, None, lat_step), lon=slice(None, None, lon_step))

    lats = np.asarray(anomaly.lat.values, dtype=float)
    lons = np.asarray(anomaly.lon.values, dtype=float)
    grid = np.asarray(anomaly.values, dtype=float)

    # OISST longitudes run 0..360; the map wants -180..180. Converting alone makes
    # the array wrap discontinuously (..., 179.875, -179.875, ...), so re-sort the
    # columns afterwards or they land in the wrong place.
    lons = np.where(lons > 180, lons - 360, lons)
    order = np.argsort(lons)
    lons = lons[order]
    grid = grid[:, order]

    # lat ascends south->north in the source; flip so row 0 is the north pole.
    if lats[0] < lats[-1]:
        lats = lats[::-1]
        grid = grid[::-1, :]

    return grid, lats, lons, target_date


# ---------------------------------------------------------------------------
# ENSO state
# ---------------------------------------------------------------------------
# The Oceanic Nino Index: a 3-month running mean of Nino 3.4 SST anomalies, and
# the thing ENSO phase is actually declared from. It cannot be derived from a
# single day's grid, and CPC serves no CORS headers, so the page cannot fetch it
# either -- hence pulling it here and carrying it in the sidecar.
#
# The strip's *value* is not taken from here: the page computes Nino 3.4 straight
# from the grid it already has, which is a day fresh rather than a month, and
# agrees with CPC's own weekly figure to about 0.02 C. This supplies the phase
# and trend, which need the running mean.
ONI_URL = "https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt"

# Standard CPC strength bands, applied to |ONI|.
ONI_BANDS = [(0.5, "Weak"), (1.0, "Moderate"), (1.5, "Strong"), (2.0, "Very strong")]


def classify_oni(oni):
    """(phase, strength) for an ONI value, using CPC's conventional thresholds."""
    if abs(oni) < 0.5:
        return "Neutral", ""
    phase = "El Niño" if oni > 0 else "La Niña"
    strength = "Weak"
    for lo, name in ONI_BANDS:
        if abs(oni) >= lo:
            strength = name
    return phase, strength


def parse_oni(text, keep=6):
    """Last `keep` seasons from CPC's oni.ascii.txt as [{season, oni}, ...]."""
    rows = []
    for line in text.splitlines():
        parts = line.split()
        if len(parts) != 4 or parts[0] == "SEAS":
            continue
        try:
            rows.append({"season": f"{parts[0]} {parts[1]}", "oni": float(parts[3])})
        except ValueError:
            continue
    return rows[-keep:]


def build_enso(fetch=None):
    """ENSO block for the sidecar, or None if CPC can't be reached.

    Deliberately non-fatal: a missing ONI should cost the strip its phase
    readout, not cost the whole run its SST grid.
    """
    try:
        if fetch is None:
            resp = requests.get(ONI_URL, timeout=30)
            resp.raise_for_status()
            text = resp.text
        else:
            text = fetch()
        recent = parse_oni(text)
        if len(recent) < 3:
            return None
        latest = recent[-1]
        phase, strength = classify_oni(latest["oni"])
        # Trend over the last three seasons, measured on |ONI| so that a
        # deepening La Nina reads as strengthening rather than falling.
        delta = abs(latest["oni"]) - abs(recent[-3]["oni"])
        if phase == "Neutral":
            trend = "Neutral"
        elif delta > 0.2:
            trend = "Strengthening"
        elif delta < -0.2:
            trend = "Weakening"
        else:
            trend = "Steady"
        return {
            "season": latest["season"],
            "oni": round(latest["oni"], 2),
            "phase": phase,
            "strength": strength,
            "trend": trend,
            "delta_3season": round(delta, 2),
            "recent": recent,
            "source": "NOAA CPC Oceanic Nino Index (ERSSTv5, 3-month running mean)",
        }
    except Exception as e:
        print(f"  WARNING: could not fetch the ONI ({e}); the strip will fall back "
              f"to the grid-derived value alone.", file=sys.stderr)
        return None


def build_meta(grid, lats, lons, target_date, step_deg, clipped_low, clipped_high, png_name):
    valid = np.isfinite(grid)
    vals = grid[valid]
    return {
        "date": target_date.strftime("%Y-%m-%d"),
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source": "NOAA OISST v2 High-Res, anomaly vs 1991-2020 daily climatology",
        "resolution_deg": step_deg,
        "grid": {
            # Regular grid, so the corner + step is all the client needs; shipping
            # the full lat/lon arrays would be another 30 KB for no information.
            "width": int(grid.shape[1]),
            "height": int(grid.shape[0]),
            "lat_max": round(float(lats[0]), 6),      # row 0
            "lat_min": round(float(lats[-1]), 6),     # last row
            "lon_min": round(float(lons[0]), 6),      # column 0
            "lon_max": round(float(lons[-1]), 6),
            "step_deg": round(float(abs(lats[0] - lats[1])), 6) if len(lats) > 1 else step_deg,
            "row_order": "north_to_south",
        },
        "encoding": {
            "file": png_name,
            "format": "png-l8-lut",
            "nodata_level": NODATA_LEVEL,
            # level -> degrees C. null marks no data. The client builds a
            # Float32Array straight off this, so the scheme above can change
            # without any JavaScript edit.
            "lut": build_lut(),
        },
        # ENSO phase/trend, or absent if CPC was unreachable. No fetch timestamp
        # here on purpose: the workflow diffs this file (minus generated_at) to
        # decide whether to commit, and a timestamp would force a commit daily.
        "enso": build_enso(),
        "stats": {
            "ocean_cells": int(valid.sum()),
            "min_c": round(float(vals.min()), 3),
            "max_c": round(float(vals.max()), 3),
            "mean_c": round(float(vals.mean()), 4),
            "clipped_low": clipped_low,
            "clipped_high": clipped_high,
        },
    }


def write_png(levels, path):
    from PIL import Image
    Image.fromarray(levels, mode="L").save(path, format="PNG", optimize=True)


def build_anomaly_json(grid, lats, lons, target_date, step_deg):
    """Legacy JSON shape (--format json), for local inspection and diffing.

    Note this flips back to south->north lat order, matching the original file.
    """
    rows = np.round(grid[::-1, :], 2)
    rows = np.where(np.isnan(rows), None, rows).tolist()
    return {
        "date": target_date.strftime("%Y-%m-%d"),
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source": "NOAA OISST v2 High-Res, anomaly vs 1991-2020 daily climatology",
        "resolution_deg": step_deg,
        "lat": [round(float(v), 3) for v in lats[::-1]],
        "lon": [round(float(v), 3) for v in lons],
        "anomaly_c": rows,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--step", type=float, default=0.25, help="output grid spacing in degrees (native OISST is 0.25)")
    parser.add_argument("--format", choices=["png", "json"], default="png",
                        help="png (default): quantized PNG + meta sidecar. json: the original single-file shape.")
    parser.add_argument("--out", help="output path (default sst_anomaly.png / sst_anomaly.json)")
    parser.add_argument("--sst-file", help="(testing) path to a local SST NetCDF instead of the live URL")
    parser.add_argument("--clim-file", help="(testing) path to a local climatology NetCDF instead of the live download")
    parser.add_argument("--refresh-climatology", action="store_true", help="force re-download of the cached climatology file")
    args = parser.parse_args()

    out = args.out or ("sst_anomaly.png" if args.format == "png" else "sst_anomaly.json")

    print("Opening SST dataset...")
    sst_ds = open_latest_sst_from_path(args.sst_file) if args.sst_file else open_latest_sst()

    print("Opening climatology dataset...")
    clim_path = args.clim_file or ensure_local_climatology(force=args.refresh_climatology)
    clim_ds = xr.open_dataset(clim_path)

    print("Computing anomaly grid...")
    grid, lats, lons, target_date = compute_anomaly(sst_ds, clim_ds, args.step)

    if args.format == "json":
        payload = build_anomaly_json(grid, lats, lons, target_date, args.step)
        with open(out, "w") as f:
            json.dump(payload, f)
        print(f"Wrote {out}: {payload['date']}, {len(payload['lat'])}x{len(payload['lon'])} grid")
        return

    levels, clipped_low, clipped_high = quantize(grid)
    write_png(levels, out)

    meta_path = os.path.splitext(out)[0] + ".meta.json"
    meta = build_meta(grid, lats, lons, target_date, args.step, clipped_low, clipped_high, os.path.basename(out))
    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=1)

    if clipped_low or clipped_high:
        print(f"  WARNING: {clipped_low} low / {clipped_high} high cells fell outside the "
              f"representable range and were clamped.", file=sys.stderr)

    size = os.path.getsize(out)
    print(f"Wrote {out}: {meta['date']}, {meta['grid']['width']}x{meta['grid']['height']} "
          f"({size/1e6:.2f} MB), {meta['stats']['ocean_cells']:,} ocean cells")
    print(f"Wrote {meta_path}: {os.path.getsize(meta_path)/1e3:.1f} KB")


if __name__ == "__main__":
    main()
