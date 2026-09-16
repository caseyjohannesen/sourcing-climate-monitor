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


# `kind` drives the panel's rendering: "index" is diverging about zero, "pct" is
# a 0-100 bar. `valid` clips the physically meaningful range, which also discards
# each product's nodata fill (none of them declare one in the GeoTIFF header).
LAYERS = [
    {"id": "spi3",  "label": "SPI (3-month)",     "kind": "index", "valid": (-5, 5),
     "url": nidis_cog("GLOBAL-ERA5_LAND_DAILY-spi-90d"),   "source": "ERA5-Land via NIDIS"},
    {"id": "spi9",  "label": "SPI (9-month)",     "kind": "index", "valid": (-5, 5),
     "url": nidis_cog("GLOBAL-ERA5_LAND_DAILY-spi-270d"),  "source": "ERA5-Land via NIDIS"},
    {"id": "spei3", "label": "SPEI (3-month)",    "kind": "index", "valid": (-5, 5),
     "url": nidis_cog("GLOBAL-ERA5_LAND_DAILY-speih-90d"), "source": "ERA5-Land via NIDIS"},
    {"id": "vhi",   "label": "Vegetation health", "kind": "pct",   "valid": (0, 100),
     "url": None, "source": "NOAA STAR Blended-VHP 4km"},   # url resolved at runtime
]


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
            return ds.read(1), ds.transform, ds.bounds

    tmp = tempfile.NamedTemporaryFile(suffix=".tif", delete=False)
    try:
        tmp.close()
        with urllib.request.urlopen(url, timeout=300) as resp, open(tmp.name, "wb") as fh:
            while chunk := resp.read(1 << 20):
                fh.write(chunk)
        with rasterio.open(tmp.name) as ds:
            return ds.read(1), ds.transform, ds.bounds
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
            arr, transform, _ = read_whole(url, download)
        except Exception as exc:
            # One bad layer should not cost the other three: the panel already
            # renders a row as unavailable when its value is missing.
            print(f"  FAILED: {exc}")
            layer_meta.append({**{k: layer[k] for k in ("id", "label", "kind", "source")},
                               "status": "unavailable", "error": str(exc)})
            continue

        lo, hi = layer["valid"]
        got = 0
        for f in features:
            iso = f["properties"][ISO_KEY]
            val = zonal_mean(arr, transform, f["geometry"], lo, hi)
            if val is not None:
                countries.setdefault(iso, {})[layer["id"]] = val
                got += 1
        print(f"  {got}/{len(features)} countries")

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
