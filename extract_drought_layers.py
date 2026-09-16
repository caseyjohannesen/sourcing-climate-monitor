"""Per-country drought indicators for the country drill-down panel.

Writes drought_country.json: one row per country (ADM0_A3 -> four values), a few
KB in total. Deliberately not gridded -- nothing renders these on the map, they
are only read as numbers inside the panel after a country is clicked, so there is
no reason to carry a raster into the repo.

The four layers, and why these specific products:

  spi3   3-month SPI, ERA5-Land          ce-GLOBAL-ERA5_LAND_DAILY-spi-90d
  spi9   9-month SPI, ERA5-Land          ce-GLOBAL-ERA5_LAND_DAILY-spi-270d
  spei3  3-month SPEI, ERA5-Land         ce-GLOBAL-ERA5_LAND_DAILY-speih-90d
  vhi    Vegetation Health Index         NOAA STAR Blended-VHP 4km

Three of those are not the products originally scoped, because the originals do
not exist as data. drought.gov's Google Cloud bucket carries GeoTIFFs for exactly
one family -- the `ce-` (Climate Engine) products. Everything else in it, the
GPCC 9-month index, the NOAA VHI layer and the SPoRT-LiS soil moisture layer
included, is published only as pre-rendered XYZ PNG tile pyramids: styled images,
with no recoverable values. So:

  * 9-month SPI is ERA5's 270-day rather than drought.gov's GPCC layer. The GPCC
    one has no GeoTIFF and is also dead -- its tiles stopped updating 2025-12-31.
  * Vegetation health is pulled from NOAA STAR directly, which is the upstream
    source drought.gov renders. Same data, real values.
  * Soil moisture is replaced by SPEI. SPoRT-LiS is a CONUS-only product (its
    tiles span lon -135..-45, lat 22..56), so per-country global values are not
    possible from it at all; no global soil-moisture raster exists in the bucket.
    SPEI is precipitation minus evaporative demand, which is the closest global
    stand-in for moisture stress, and it shares SPI's grid and standardized scale.

Scales differ and must not be conflated: SPI and SPEI are unitless standardized
anomalies centred on zero (roughly -3..+3, negative = drier than normal), while
VHI is a bounded 0-100 percentile-like index (below 40 stressed, above 60
favourable). The panel renders the two kinds differently for that reason.

Usage:
    python extract_drought_layers.py [--out drought_country.json]
"""

import argparse
import json
import os
import tempfile
import urllib.error
import urllib.request
from datetime import date, datetime, timezone

import numpy as np

# vsicurl otherwise issues a directory listing per open, which on a bucket this
# large costs more than the read itself.
os.environ.setdefault("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
os.environ.setdefault("GDAL_HTTP_TIMEOUT", "120")

import rasterio                               # noqa: E402  (after the GDAL env)
from rasterio.features import geometry_mask   # noqa: E402
from rasterio.windows import from_bounds      # noqa: E402

BOUNDARIES = "countries.geojson"
OUT = "drought_country.json"
ISO_KEY = "ADM0_A3"   # always populated, unlike ISO_A3 -- see prepare_boundaries.py

NIDIS = ("https://storage.googleapis.com/noaa-nidis-drought-gov-data"
         "/current-conditions/tile/v1")
STAR = ("https://www.star.nesdis.noaa.gov/pub/corp/scsb/wguo/data"
        "/Blended_VH_4km/geo_TIFF")


def nidis_cog(product):
    return f"{NIDIS}/ce-{product}/{product}.tif"


# `kind` drives rendering: "index" is diverging about zero, "pct" is a 0-100
# scale.
#
# Two ranges, deliberately not one:
#   `valid` is physical plausibility, used only to throw away each product's
#     nodata fill (none declare one in the GeoTIFF header). It is wide on
#     purpose -- NIDIS floors SPI at -4.00 but leaves the wet tail uncapped, and
#     real cells reach +8.2, so a tight window here would silently drop genuine
#     extremes out of the country means.
#   `clip` is the display range the 8-bit overlay quantizes across. Values
#     outside it are clamped, not dropped, and the count is recorded in the
#     sidecar.
LAYERS = [
    {"id": "spi3",  "label": "SPI (3-month)",     "kind": "index",
     "valid": (-10, 10), "clip": (-4, 4),
     "url": nidis_cog("GLOBAL-ERA5_LAND_DAILY-spi-90d"),   "source": "ERA5-Land via NIDIS"},
    {"id": "spi9",  "label": "SPI (9-month)",     "kind": "index",
     "valid": (-10, 10), "clip": (-4, 4),
     "url": nidis_cog("GLOBAL-ERA5_LAND_DAILY-spi-270d"),  "source": "ERA5-Land via NIDIS"},
    {"id": "spei3", "label": "SPEI (3-month)",    "kind": "index",
     "valid": (-10, 10), "clip": (-4, 4),
     "url": nidis_cog("GLOBAL-ERA5_LAND_DAILY-speih-90d"), "source": "ERA5-Land via NIDIS"},
    {"id": "vhi",   "label": "Vegetation health", "kind": "pct",
     "valid": (0, 100),  "clip": (0, 100),
     "url": None, "source": "NOAA STAR Blended-VHP 4km"},   # url resolved at runtime
]

# ---------------------------------------------------------------------------
# overlay grids
# ---------------------------------------------------------------------------
# The overlay PNGs are written onto OISST's exact grid rather than the rasters'
# native one. That is worth the resample for three reasons: the drought products
# only span 75N-75S and their cell centres sit half a cell off any pole-aligned
# grid, so landing them here makes the polar nodata padding fall out for free;
# the page's texture, region tool and sampling then need no per-layer geometry;
# and at 0.1 deg native each layer is 748 KB against 131 KB here, on a repo that
# takes a data commit every day.
GRID_DIR = "layers"
GRID_W, GRID_H, GRID_STEP = 1440, 720, 0.25
NODATA_LEVEL = 0        # matches extract_sst_anomaly.py: level 0 decodes to NaN
LEVELS = 255            # 1..255 carry data


def grid_meta_block():
    half = GRID_STEP / 2
    return {
        "width": GRID_W, "height": GRID_H,
        "lat_max": 90 - half, "lat_min": -90 + half,
        "lon_min": -180 + half, "lon_max": 180 - half,
        "step_deg": GRID_STEP,
        "row_order": "north_to_south",
    }


def build_lut(lo, hi):
    """256-entry decode table. Index 0 is nodata; 1..255 span [lo, hi] linearly.

    SPI and SPEI need nothing like the SST extractor's piecewise LUT: that one
    exists to spend levels on a long sea-ice tail, whereas these are already
    clipped by the provider and sit almost entirely within +/-3. A flat linear
    ramp over +/-4 gives 0.0315 per level, finer than SST's 0.05 core step.
    """
    step = (hi - lo) / (LEVELS - 1)
    return [None] + [round(lo + i * step, 6) for i in range(LEVELS)]


def quantize(values, lo, hi):
    """Float array -> uint8 levels, NaN -> NODATA_LEVEL.

    Returns (levels, clipped_low, clipped_high) so the sidecar can record how
    much real data fell outside the display range.
    """
    valid = np.isfinite(values)
    v = np.clip(np.where(valid, values, lo), lo, hi)
    levels = np.round((v - lo) / (hi - lo) * (LEVELS - 1)) + 1
    levels = np.where(valid, levels, NODATA_LEVEL)
    return (levels.astype(np.uint8),
            int(np.sum(valid & (values < lo))),
            int(np.sum(valid & (values > hi))))


def regrid(values, src_transform, src_crs):
    """Area-average `values` onto the canonical 0.25 deg grid."""
    from rasterio.transform import from_origin
    from rasterio.warp import Resampling, reproject

    dst = np.full((GRID_H, GRID_W), np.nan, dtype="float32")
    reproject(
        source=values.astype("float32"), destination=dst,
        src_transform=src_transform, src_crs=src_crs, src_nodata=np.nan,
        dst_transform=from_origin(-180.0, 90.0, GRID_STEP, GRID_STEP),
        dst_crs="EPSG:4326", dst_nodata=np.nan,
        resampling=Resampling.average,
    )
    return dst


def write_grid(layer, values, src_transform, src_crs, date, out_dir):
    """Write <id>.png and <id>.meta.json, the same pair extract_sst_anomaly.py
    emits, so the page's existing loader reads them with no special case."""
    from PIL import Image

    os.makedirs(out_dir, exist_ok=True)
    lo, hi = layer["clip"]
    dst = regrid(values, src_transform, src_crs)
    levels, clipped_low, clipped_high = quantize(dst, lo, hi)

    png_name = f"{layer['id']}.png"
    Image.fromarray(levels, mode="L").save(
        os.path.join(out_dir, png_name), format="PNG", optimize=True)

    finite = dst[np.isfinite(dst)]
    meta = {
        "date": date,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source": layer["source"],
        "layer": {"id": layer["id"], "label": layer["label"], "kind": layer["kind"]},
        "resolution_deg": GRID_STEP,
        "grid": grid_meta_block(),
        "encoding": {
            "file": png_name,
            "format": "png-l8-lut",
            "nodata_level": NODATA_LEVEL,
            "lut": build_lut(lo, hi),
        },
        "stats": {
            "cells": int(finite.size),
            "min": round(float(finite.min()), 4) if finite.size else None,
            "max": round(float(finite.max()), 4) if finite.size else None,
            "mean": round(float(finite.mean()), 4) if finite.size else None,
            "clipped_low": clipped_low,
            "clipped_high": clipped_high,
        },
    }
    with open(os.path.join(out_dir, f"{layer['id']}.meta.json"), "w", encoding="utf-8") as fh:
        json.dump(meta, fh, separators=(",", ":"))

    size = os.path.getsize(os.path.join(out_dir, png_name))
    return size, meta["stats"]


# ---------------------------------------------------------------------------
# VHI: find the newest published week
# ---------------------------------------------------------------------------
def head(url, timeout=20):
    """HEAD a URL, returning the response headers, or None if it is not there.

    Timeout is deliberately short: the VHI search below can issue two dozen of
    these, so a slow default would let one unreachable host stall the whole job
    for many minutes.
    """
    req = urllib.request.Request(url, method="HEAD")
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return resp.headers if resp.status == 200 else None
    except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, OSError):
        return None


def published_date(headers):
    """Last-Modified as a plain YYYY-MM-DD, or None.

    Worth carrying into the output because these products do not all refresh on
    the same cadence -- the 90-day SPI lands daily while the 270-day one can sit
    for weeks -- and the panel should not imply they are equally current.
    """
    raw = headers.get("Last-Modified") if headers else None
    if not raw:
        return None
    try:
        return datetime.strptime(raw, "%a, %d %b %Y %H:%M:%S %Z").date().isoformat()
    except ValueError:
        return None


def latest_vhi_url(max_lookback=8):
    """Newest Blended-VHP weekly tile, found by probing backwards.

    STAR's directory index is JavaScript-rendered, so it cannot be scraped for a
    file list; probing is what is left. The satellite token is part of the name
    and changes as platforms are retired (npp -> j01 -> j02), so each week is
    tried against all three rather than pinning one and breaking on the next
    handover. Weeks run 1..53 and the product lags real time by about two weeks.
    """
    year, week, _ = date.today().isocalendar()
    for _ in range(max_lookback):
        for sat in ("j02", "j01", "npp"):
            url = f"{STAR}/VHP.G04.C07.{sat}.P{year}{week:03d}.VH.VHI.tif"
            headers = head(url)
            if headers is not None:
                return url, f"{year}week{week:02d}", published_date(headers)
        week -= 1
        if week < 1:                      # step back across the year boundary
            year -= 1
            week = date(year, 12, 28).isocalendar()[1]   # 28 Dec is always in the last ISO week
    raise RuntimeError(f"no VHI tile found in the last {max_lookback} weeks")


# ---------------------------------------------------------------------------
# raster access
# ---------------------------------------------------------------------------
def read_whole(url, download):
    """Read a raster into memory, with its transform.

    Reading the entire global grid once and masking countries out of memory beats
    a windowed read per country here, which is the opposite of the OISST pattern
    and worth saying why: OPeNDAP windowing wins because it slices a multi-decade
    archive, whereas each of these is a single global snapshot small enough that
    177 windowed reads just pay 177 round-trips and re-fetch overlapping tiles.
    Measured over all 177 countries: 16.8s windowed against 2.4s for one whole
    read. The VHI file settles it anyway -- it is stripe-organised rather than
    tiled, so a bounding-box window would pull full-width scanlines regardless.

    `download` fetches to a temp file first. The NIDIS COGs stream fine over
    vsicurl; STAR's file is neither tiled nor overviewed and errored partway
    through a streamed read, so it is pulled whole and opened locally.
    """
    if not download:
        with rasterio.open(f"/vsicurl/{url}") as ds:
            return ds.read(1), ds.transform, ds.crs

    tmp = tempfile.NamedTemporaryFile(suffix=".tif", delete=False)
    try:
        tmp.close()
        with urllib.request.urlopen(url, timeout=300) as resp, open(tmp.name, "wb") as fh:
            while chunk := resp.read(1 << 20):
                fh.write(chunk)
        with rasterio.open(tmp.name) as ds:
            return ds.read(1), ds.transform, ds.crs
    finally:
        os.unlink(tmp.name)


def geom_bounds(geom):
    xs, ys = [], []

    def walk(c):
        if isinstance(c[0], (int, float)):
            xs.append(c[0]); ys.append(c[1])
        else:
            for part in c:
                walk(part)

    walk(geom["coordinates"])
    return min(xs), min(ys), max(xs), max(ys)


def zonal_mean(arr, transform, geom, lo, hi):
    """cos(lat)-weighted mean of `arr` inside `geom`.

    Area weighting matches nino34_index() in archive.py: on a lat/lon grid a cell
    at 60 degrees covers half the ground a cell at the equator does, so an
    unweighted mean would over-count high latitudes. Rasterising only the
    country's bounding-box window rather than the full grid keeps this cheap --
    otherwise every country would burn a 3601x1501 mask.

    all_touched=True so small countries still catch a cell; at 0.1 deg several
    are narrower than one cell across and would otherwise come back null.
    """
    win = from_bounds(*geom_bounds(geom), transform).round_offsets().round_lengths()
    # Clip to the raster: a country outside the layer's coverage (VHI stops at
    # 55S) yields an empty window rather than an out-of-range read.
    h, w = arr.shape
    r0, c0 = max(0, win.row_off), max(0, win.col_off)
    r1 = min(h, win.row_off + win.height)
    c1 = min(w, win.col_off + win.width)
    if r1 <= r0 or c1 <= c0:
        return None

    sub = arr[r0:r1, c0:c1]
    sub_tr = rasterio.windows.transform(
        rasterio.windows.Window(c0, r0, c1 - c0, r1 - r0), transform)

    inside = geometry_mask([geom], sub.shape, sub_tr, invert=True, all_touched=True)
    lats = sub_tr.f + (np.arange(sub.shape[0]) + 0.5) * sub_tr.e
    weights = np.cos(np.radians(np.clip(lats, -89.9, 89.9)))[:, None]

    ok = inside & np.isfinite(sub) & (sub >= lo) & (sub <= hi)
    if not ok.any():
        return None
    wt = np.broadcast_to(weights, sub.shape)[ok]
    if wt.sum() <= 0:
        return None
    return round(float((sub[ok] * wt).sum() / wt.sum()), 3)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--boundaries", default=BOUNDARIES)
    ap.add_argument("--grids", default=GRID_DIR, help="directory for the overlay PNG/meta pairs")
    args = ap.parse_args()

    with open(args.boundaries, encoding="utf-8") as fh:
        features = json.load(fh)["features"]
    print(f"{len(features)} countries from {args.boundaries}")

    countries = {}
    layer_meta = []

    for layer in LAYERS:
        url, stamp, download = layer["url"], None, False
        if layer["id"] == "vhi":
            url, stamp, updated = latest_vhi_url()
            download = True
            print(f"  VHI resolved to {stamp}")
        else:
            updated = published_date(head(url))

        print(f"[{layer['id']}] reading {url.rsplit('/', 1)[-1]} ...", flush=True)
        try:
            arr, transform, crs = read_whole(url, download)
        except Exception as exc:
            # One bad layer should not cost the other three: the panel already
            # renders a row as unavailable when its value is missing.
            print(f"  FAILED: {exc}")
            layer_meta.append({**{k: layer[k] for k in ("id", "label", "kind", "source")},
                               "status": "unavailable", "error": str(exc)})
            continue

        # Mask nodata fill once, up front, so the country means and the overlay
        # grid are built from exactly the same values.
        lo, hi = layer["valid"]
        arr = np.where(np.isfinite(arr) & (arr >= lo) & (arr <= hi), arr, np.nan)

        got = 0
        for f in features:
            iso = f["properties"][ISO_KEY]
            val = zonal_mean(arr, transform, f["geometry"], lo, hi)
            if val is not None:
                countries.setdefault(iso, {})[layer["id"]] = val
                got += 1
        print(f"  {got}/{len(features)} countries")

        size, gstats = write_grid(layer, arr, transform, crs, updated, args.grids)
        print(f"  overlay {GRID_W}x{GRID_H} -> {size/1024:.0f} KB"
              f" (clipped {gstats['clipped_low']} low / {gstats['clipped_high']} high)")

        meta = {k: layer[k] for k in ("id", "label", "kind", "source")}
        meta["status"] = "ok"
        if updated:
            meta["updated"] = updated
        if stamp:
            meta["week"] = stamp
        layer_meta.append(meta)
        print(f"  published {updated or 'unknown'}")

    out = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "layers": layer_meta,
        "countries": countries,
    }
    with open(args.out, "w", encoding="utf-8") as fh:
        json.dump(out, fh, separators=(",", ":"), sort_keys=True)

    size = os.path.getsize(args.out)
    print(f"Wrote {args.out}: {len(countries)} countries, {size/1024:.1f} KB")


if __name__ == "__main__":
    main()
