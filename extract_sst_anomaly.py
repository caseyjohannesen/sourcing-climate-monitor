"""
Pull the latest daily SST anomaly from NOAA's OISST v2 High-Res dataset
(PSL/NOAA, 0.25 deg, daily since 1981) and write a small JSON file the
climate-risk-monitor globe can fetch directly.

Data source: https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.highres.html
Daily SST:          .../Datasets/noaa.oisst.v2.highres/sst.day.mean.<year>.nc
Daily climatology:  .../Datasets/noaa.oisst.v2.derived/sst.day.ltm.1991-2020.nc

Both are read via OPeNDAP (the "dodsC" THREDDS path), so xarray only pulls the
one time slice it needs over the network rather than downloading the whole
multi-gigabyte yearly file.

Usage:
    pip install xarray netCDF4 pandas numpy
    python extract_sst_anomaly.py [--step 1.0] [--out sst_anomaly.json]

Note: this has been tested against synthetic NetCDF files that mirror the
real dataset's structure (see the accompanying test), but not against NOAA's
live server directly, since that network path isn't reachable from the
environment this was written in. The URLs below are taken from NOAA's own
published directory listing, but it's worth confirming the first live run
succeeds before scheduling it unattended.
"""

import argparse
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
CLIM_DOWNLOAD_URL = "https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2.derived/sst.day.ltm.1991-2020.nc"
CLIM_CACHE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sst.day.ltm.1991-2020.nc")


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
    print(f"Downloading climatology file (one-time, ~90MB): {CLIM_DOWNLOAD_URL}")
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


def match_climatology(clim_ds, target_date):
    """Select the climatology slice matching target_date's day-of-year.

    The daily LTM file is indexed by a synthetic 366-step 'time' axis in
    calendar order (Jan 1 .. Dec 31, including Feb 29). We match by position
    (day-of-year - 1) rather than by actual date, since the LTM's own year
    is a placeholder, not a real year.
    """
    doy = target_date.dayofyear
    n = clim_ds.sizes["time"]
    idx = min(doy - 1, n - 1)  # clamp for non-leap years / 365-day LTMs
    return clim_ds.sst.isel(time=idx)


def build_anomaly_json(sst_ds, clim_ds, step_deg):
    latest = sst_ds.sst.isel(time=-1)
    target_date = pd.Timestamp(sst_ds.time.isel(time=-1).values)

    clim_slice = match_climatology(clim_ds, target_date)

    # align climatology grid to the SST grid in case resolutions/labels differ slightly
    clim_slice = clim_slice.interp(lat=latest.lat, lon=latest.lon, method="nearest")

    anomaly = latest - clim_slice

    # downsample to a coarser grid for a lightweight payload
    lat_vals = latest.lat.values
    lon_vals = latest.lon.values
    lat_step = max(1, int(round(step_deg / abs(lat_vals[1] - lat_vals[0]))))
    lon_step = max(1, int(round(step_deg / abs(lon_vals[1] - lon_vals[0]))))

    anomaly_ds = anomaly.isel(lat=slice(None, None, lat_step), lon=slice(None, None, lon_step))

    lats = [round(float(v), 2) for v in anomaly_ds.lat.values]
    lons = np.array(anomaly_ds.lon.values)
    # convert lon from 0..360 (OISST convention) to -180..180 (matches the globe's projection),
    # then re-sort so columns stay in ascending longitude order -- converting alone makes the
    # array wrap discontinuously (e.g. ...,178.5, -179.5,...) which would misplace columns.
    lons_180 = np.where(lons > 180, lons - 360, lons)
    order = np.argsort(lons_180)
    lons_180 = [round(float(v), 2) for v in lons_180[order]]

    grid = np.round(anomaly_ds.values, 2)
    grid = grid[:, order]
    grid = np.where(np.isnan(grid), None, grid).tolist()

    return {
        "date": target_date.strftime("%Y-%m-%d"),
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source": "NOAA OISST v2 High-Res, anomaly vs 1991-2020 daily climatology",
        "resolution_deg": step_deg,
        "lat": lats,
        "lon": lons_180,
        "anomaly_c": grid,  # rows follow lat, columns follow lon; null = no data (land / gap)
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--step", type=float, default=1.0, help="output grid spacing in degrees")
    parser.add_argument("--out", default="sst_anomaly.json")
    parser.add_argument("--sst-file", help="(testing) path to a local SST NetCDF instead of the live URL")
    parser.add_argument("--clim-file", help="(testing) path to a local climatology NetCDF instead of the live download")
    parser.add_argument("--refresh-climatology", action="store_true", help="force re-download of the cached climatology file")
    args = parser.parse_args()

    print("Opening SST dataset...")
    sst_ds = open_latest_sst_from_path(args.sst_file) if args.sst_file else open_latest_sst()

    print("Opening climatology dataset...")
    clim_path = args.clim_file or ensure_local_climatology(force=args.refresh_climatology)
    clim_ds = xr.open_dataset(clim_path)

    print("Computing anomaly grid...")
    payload = build_anomaly_json(sst_ds, clim_ds, args.step)

    with open(args.out, "w") as f:
        json.dump(payload, f)

    print(f"Wrote {args.out}: {payload['date']}, {len(payload['lat'])}x{len(payload['lon'])} grid")


if __name__ == "__main__":
    main()
