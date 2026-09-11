# sourcing-climate-monitor

A single-page 3D globe (Three.js) showing daily **sea surface temperature (SST)
anomaly** data — built as a cocoa-sourcing climate-risk view for a chocolate
manufacturer. The data is refreshed automatically each day by a GitHub Actions
workflow and served over GitHub Pages.

## How it works

- **`extract_sst_anomaly.py`** pulls the latest daily SST from NOAA's OISST v2
  High-Res dataset and subtracts the matching 1991-2020 daily climatology at the
  dataset's **native 0.25° resolution** (1440×720). Both grids come from the same
  `highres` directory, so they subtract directly with no regridding.
- The output is a pair: **`sst_anomaly.png`** (8-bit greyscale, one pixel per
  cell, row 0 = north, level 0 = land / no data) and **`sst_anomaly.meta.json`**
  (~3 KB: date, source, grid geometry, and a 256-entry lookup table turning pixel
  levels back into °C).
- **`index.html`** fetches the sidecar and then the PNG it names, decodes it into
  a `Float32Array`, and renders it onto the globe. If the fetch fails (opened
  locally, or before the first data commit) it falls back to a manual file picker.
- The view toggles between a **3D globe** and a flat **Equal Earth** map
  (Šavrič/Jenny/Patterson 2018), via `d3-geo`'s `geoEqualEarth`.
- **`.github/workflows/refresh.yml`** runs the script daily at 13:00 UTC and
  commits the refreshed pair back to the repo when the data actually changes.

### Why a PNG rather than JSON

At native resolution the grid is 1,036,800 cells. As JSON that is ~6.3 MB
(~840 KB gzipped); quantized into a PNG it is ~250 KB, and the browser decodes it
natively instead of parsing a million boxed JS numbers.

The quantization is piecewise — 0.05 °C steps across |anomaly| ≤ 5, 0.5 °C steps
in the tails out to −18/+19 °C. A flat 8-bit ramp can't do both jobs: real fields
reach +16 °C at the sea-ice edge, so a 0.05 °C ramp would clip and be wrong by
9 °C in the Arctic. Measured against a real global field, RMS round-trip error is
0.02 °C and the worst case is 0.49 °C, confined to those far tails. The display
colour bands are 0.2–1.0 °C wide, so the error sits well below anything the map
can show. `--format json` still emits the original single-file shape for local
inspection.

### Region selection

Right-click the map to drop pins; right-click the first pin again (or press
**Close region**) to close the ring and get the area-weighted mean SST anomaly
inside it. Backspace undoes a pin, Esc clears. It works identically on the globe
and the flat map, and a region drawn in one view is still there in the other.

Two details make this spherical rather than planar geometry:

- **Edges are great-circle arcs**, via `d3.geoInterpolate`. A great circle from
  (0°E, 60°N) to (90°E, 60°N) reaches **67.79°N** at its midpoint — a straight
  lon/lat interpolation says 60°N the whole way, an ~870 km error on one edge.
  The region's bounding box is therefore taken from the densified arc, never
  from the pins, or the scan would clip the bulge.
- **Point-in-polygon is done in a gnomonic frame** centred on the ring's own
  centroid. Gnomonic maps every great circle to a straight line *exactly*, so
  the spherical polygon becomes a planar one and the ordinary crossing test is
  correct — at ~50 ns per cell rather than the ~10 µs `d3.geoContains` costs,
  which matters when scanning tens of thousands of cells. It also settles the
  "which side is inside?" ambiguity the right way: inside is the side containing
  the centroid. (`d3.geoArea` on a clockwise ring returns the area of the entire
  rest of the planet — winding order is load-bearing, and this avoids depending
  on it.) The implementation is checked against `d3.geoContains` as an oracle.

Gnomonic diverges at 90° from its centre, so a ring spanning more than ~85°
from its own centroid is refused with a message rather than answered wrongly.

### The projection seam

Both views are the *same* scene. All geometry — country borders, the hover
highlight, the region outline — is authored in lon/lat and only becomes a
position at the last moment, through one matched pair:

```
project<Mode>(lon, lat, elev, out, o)   lon/lat -> render position
unproject(clientX, clientY)             screen  -> lon/lat (or null)
```

The surface itself is a lon/lat grid mesh carrying equirectangular UVs — which
is what `THREE.SphereGeometry` already is — so switching projection only moves
vertices. The UVs are identical either way, which means **the SST texture is
never resampled**; the GPU interpolates the same equirectangular raster across
either shape. Each object keeps both projections' vertex arrays, so the toggle
animates as a lerp rather than a reprojection.

`unproject` is analytic — ray-vs-sphere for the globe, ray-vs-plane plus
`equalEarth.invert()` for the flat map — rather than a mesh raycast, so pointer
handling costs the same no matter how finely the surface is tessellated. Hover,
click and the country panel are projection-agnostic: they were not modified when
the flat map was added.

## Running locally

```bash
pip install xarray netCDF4 pandas numpy requests pillow
python extract_sst_anomaly.py --out sst_anomaly.png
```

Then serve the folder (a plain file open won't allow the `fetch`, so use a
server) and visit it:

```bash
python -m http.server 8000
# open http://localhost:8000/
```

The first run downloads the ~350 MB 0.25° 1991-2020 climatology once and caches
it as `sst.day.mean.ltm.1991-2020.nc` (git-ignored).
