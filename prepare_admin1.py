"""One-time prep: split Natural Earth 10m admin-1 (states / provinces / regions)
into one small GeoJSON per country, so the page can fetch a country's districts
on demand instead of carrying the whole 40 MB dataset.

The source file is ~40 MB and covers 251 countries. Committing it whole would be
absurd, and so would loading it in the browser; split and simplified it comes to
about 7 MB across ~250 files, of which a click costs one -- Ghana is 15 KB, and
the median country is under 20 KB. Those go in the repo alongside
countries.geojson, for the same reason that one is committed: the page then
never depends on a third-party service at runtime.

Simplification is Douglas-Peucker at 0.01 deg (~1.1 km). Each district is
simplified independently, which in principle lets a border shared by two
districts drift apart and draw as a doubled line. Measured, the slack is
0.006-0.064% of a country's area -- roughly 35 m along a shared border, well
under a pixel even zoomed right in -- so it isn't worth carrying a topology
library to avoid.

Note Natural Earth's admin-1 is not always current: Ghana is 10 regions here,
not the 16 it has had since 2018. geoBoundaries is fresher but ~18x larger per
country and would put two external services in the runtime path.

Source: Natural Earth (public domain) via the nvkelso/natural-earth-vector mirror.

Usage:
    pip install shapely
    python prepare_admin1.py [--tolerance 0.01] [--src <cached .geojson>]
"""

import argparse
import json
import os
import urllib.request
from collections import defaultdict

from shapely.geometry import shape, mapping

SRC = ("https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/"
       "geojson/ne_10m_admin_1_states_provinces.geojson")
CACHE = ".cache_ne_10m_admin_1.geojson"
OUT_DIR = "admin1"
DECIMALS = 3          # ~110 m, below the simplification tolerance anyway


def fetch_source(src_path):
    if src_path:
        return src_path
    if os.path.exists(CACHE):
        print(f"Using cached source: {CACHE}")
        return CACHE
    print(f"Downloading {SRC} (~40 MB, one time) ...")
    urllib.request.urlretrieve(SRC, CACHE)
    return CACHE


def round_geom(geom, n):
    t, c = geom["type"], geom["coordinates"]
    r = lambda ring: [[round(x, n), round(y, n)] for x, y in ring]
    if t == "Polygon":
        return {"type": t, "coordinates": [r(x) for x in c]}
    if t == "MultiPolygon":
        return {"type": t, "coordinates": [[r(x) for x in poly] for poly in c]}
    return geom


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tolerance", type=float, default=0.01,
                    help="Douglas-Peucker tolerance in degrees (default 0.01, ~1.1 km)")
    ap.add_argument("--src", help="path to an already-downloaded source geojson")
    ap.add_argument("--out", default=OUT_DIR)
    args = ap.parse_args()

    with open(fetch_source(args.src), encoding="utf-8") as fh:
        src = json.load(fh)

    by_country = defaultdict(list)
    for f in src["features"]:
        # adm0_a3 is always populated, unlike iso_a3 which is -99 for a handful
        # of countries (France and Norway among them).
        by_country[f["properties"].get("adm0_a3") or "UNK"].append(f)

    os.makedirs(args.out, exist_ok=True)
    index, total = {}, 0
    for iso, feats in sorted(by_country.items()):
        out_feats = []
        for f in feats:
            g = shape(f["geometry"]).buffer(0).simplify(args.tolerance, preserve_topology=True)
            if g.is_empty:
                continue
            p = f["properties"]
            out_feats.append({
                "type": "Feature",
                "properties": {"name": p.get("name"), "type": p.get("type_en")},
                "geometry": round_geom(mapping(g), DECIMALS),
            })
        if not out_feats:
            continue
        blob = json.dumps({"type": "FeatureCollection", "features": out_feats},
                          separators=(",", ":"), ensure_ascii=False)
        path = os.path.join(args.out, f"{iso}.geojson")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(blob)
        index[iso] = len(out_feats)
        total += len(blob.encode("utf-8"))

    # A manifest so the page can say "no district data for this country" without
    # first firing off a request that 404s.
    with open(os.path.join(args.out, "index.json"), "w", encoding="utf-8") as fh:
        json.dump(index, fh, separators=(",", ":"), sort_keys=True)

    print(f"Wrote {len(index)} country files to {args.out}/ "
          f"({total/1e6:.2f} MB, {sum(index.values())} districts)")
    biggest = sorted(index, key=lambda k: os.path.getsize(os.path.join(args.out, k + ".geojson")))[-5:]
    for iso in reversed(biggest):
        print(f"  largest: {iso} {os.path.getsize(os.path.join(args.out, iso + '.geojson'))/1e3:.0f} KB "
              f"({index[iso]} districts)")


if __name__ == "__main__":
    main()
