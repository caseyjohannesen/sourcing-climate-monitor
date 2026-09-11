"""One-time prep: download Natural Earth 110m admin-0 countries and slim it into
a small GeoJSON the globe fetches at load (countries.geojson, committed to the
repo so the page never depends on a third-party CDN at runtime).

The slim keeps only the handful of properties the UI needs and rounds coordinates
to 3 decimals (~100 m), taking the file from ~820 KB to ~200 KB (~67 KB gzipped).
Re-run only if you want to refresh the boundaries.

Source: Natural Earth (public domain) via the nvkelso/natural-earth-vector mirror.

Usage:
    python prepare_boundaries.py
"""

import json
import urllib.request

SRC = "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_110m_admin_0_countries.geojson"
OUT = "countries.geojson"
# ADM0_A3 is the key used to look up a country's admin-1 file (see
# prepare_admin1.py). It is kept alongside ISO_A3 rather than instead of it
# because Natural Earth leaves ISO_A3 as "-99" for a handful of countries --
# France, Norway, N. Cyprus, Somaliland and Kosovo -- while ADM0_A3 is always set.
KEEP = ["NAME", "ADMIN", "ISO_A3", "ADM0_A3", "CONTINENT", "LABEL_X", "LABEL_Y"]


def _round_ring(ring):
    return [[round(x, 3), round(y, 3)] for x, y in ring]


def round_coords(geom):
    t = geom["type"]
    c = geom["coordinates"]
    if t == "Polygon":
        return {"type": t, "coordinates": [_round_ring(r) for r in c]}
    if t == "MultiPolygon":
        return {"type": t, "coordinates": [[_round_ring(r) for r in poly] for poly in c]}
    return geom


def main():
    print(f"Downloading {SRC} ...")
    with urllib.request.urlopen(SRC) as resp:
        src = json.loads(resp.read().decode("utf-8"))

    out = {"type": "FeatureCollection", "features": []}
    for f in src["features"]:
        p = f["properties"]
        out["features"].append({
            "type": "Feature",
            "properties": {k: p.get(k) for k in KEEP},
            "geometry": round_coords(f["geometry"]),
        })

    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(out, fh, separators=(",", ":"), ensure_ascii=False)
    print(f"Wrote {OUT}: {len(out['features'])} countries")


if __name__ == "__main__":
    main()
