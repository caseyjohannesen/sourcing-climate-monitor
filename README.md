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

### Playback and the Niño 3.4 trend

A timeline under the map scrubs through the last 180 days; the sparkline in the
ENSO strip shows the Niño 3.4 index over the same window, with a marker that
tracks whatever day is on screen.

Frames live in `archive/frames/<date>.png` at 0.5°, sharing one grid and LUT from
`archive/index.json` rather than a sidecar each — 15.6 MB for 180 days against
45 MB at native, which is the only reason playback is practical in a repo that
also takes a data commit daily. Only the two small JSON files load up front;
frames are fetched on demand, cached to a bounded 40, and each displayed texture
disposes the one before it (180 undisposed would be ~180 MB of GPU memory).

The rightmost position is *live*: that day's native 0.25° grid is already loaded,
so it is shown rather than its downsampled copy. The bar says which you are
looking at, so the softer archive frames don't read as a bug.

**The series is computed at native resolution even though the frames are not.**
For 2026-09-09 the native index is +2.8877 — exactly what the browser's
`averageBox()` gives — where the 0.5° frame gives +2.8713. Scrubbing therefore
shows the archived native value, so the figure beside the sparkline always
matches the point the marker sits on.

Smoothing is a **7-day trailing mean**. Measured over 186 real days, raw daily
noise is 0.032 °C sd against a 2.9 °C seasonal swing — sub-pixel in a sparkline —
while a 90-day mean lags by 0.50 °C during a fast-developing event, which
understates risk. The window is counted in calendar days, not array positions, so
a gap shortens it rather than silently reaching further back than it claims to.

`backfill_archive.py` fills the window from history instead of waiting six months
for it to accumulate (`Backfill frame archive` workflow, manual). Cross-checked
against ERDDAP's independently-computed anomalies: correlation 0.9982 with a
near-constant offset, which is what two climatology baselines over the same water
should look like.

### The ENSO strip

The header strip is live, from two sources, because the three readouts need
different things.

The **anomaly** is computed in the browser from the grid already on screen — an
area-weighted mean of the Niño 3.4 box (5°S–5°N, 170–120°W). That costs nothing,
is a day fresh rather than a month, and is consistent with what the map shows.
Checked against NOAA CPC's own weekly Niño 3.4 figure for the same week, the two
agree to **0.02 °C** (ours +2.578, CPC +2.6).

**Phase** and **3-month trend** cannot come from a single day — they are declared
from the Oceanic Niño Index, a 3-month running mean — and CPC serves no CORS
headers, so the page cannot fetch it. `extract_sst_anomaly.py` pulls
`oni.ascii.txt`, classifies it with CPC's conventional thresholds, and leaves the
result in the sidecar. Trend is measured on |ONI| so a deepening La Niña reads as
strengthening rather than falling.

If CPC is unreachable the run still succeeds: the ENSO block is simply absent,
the anomaly still shows, and phase and trend read "—" with the reason given
rather than displaying something stale.

### Clicking a country

Clicking a country flies the view to it — rotating and dollying the globe, or
panning and zooming the flat map — loads that country's admin-1 district
boundaries, and shows its **area-averaged** air temperature, precipitation and
soil moisture. Clicking again inside a district drills down to that district's
own average. Esc steps back out one level at a time (district → country → world),
as does the back link in the panel; dragging or scrolling takes over at any point.

The averages come from Open-Meteo, which accepts many coordinates in one request
and exposes soil moisture as a `current` variable — so one request per selection
covers the whole region. Sample points are laid on a lat/lon grid inside the
polygon and **cos(lat)-weighted**, since a regular lon/lat grid over-samples
toward the poles. `timezone=GMT` is deliberate: with `timezone=auto` each point
resolves to its own local clock, so averaging a wide country would mix local
morning with local evening, and the mean would mean nothing.

This is a 25-point sample of a weather model, not an areal integral — dense for
Ghana, coarse for Russia — so the panel always states the point count rather than
letting "average" imply more than it is. Results are cached per country and per
district, and requests are strictly click-driven.

Note these files are fetched with `cache: 'default'`, deliberately not
`force-cache`: the latter returns a cached match *fresh or stale* and never
revalidates, so a browser that had loaded an earlier release would keep serving
it indefinitely. That is not hypothetical — it is exactly what broke districts
on first release, when `countries.geojson` gained the `ADM0_A3` property that
districts are keyed by and already-cached copies never picked it up.

Districts live in `admin1/<ADM0_A3>.geojson`, one file per country, generated by
`prepare_admin1.py` from Natural Earth 10m admin-1. The whole dataset is ~40 MB;
split and simplified to ~1.1 km it is ~7 MB total, of which a click costs one
file — Ghana is 15 KB, the median country under 20 KB. `admin1/index.json` is a
small manifest so the panel can say "no districts" without a request that 404s.

Note `prepare_boundaries.py` keeps `ADM0_A3` as the join key, not `ISO_A3`:
Natural Earth leaves `ISO_A3` as `-99` for France, Norway, N. Cyprus, Somaliland
and Kosovo.

The camera's near limit is per-projection. The globe must keep the camera
outside a sphere of radius `GLOBE_R`, so it stops at 1.6; the flat map has
nothing in the way and goes to 0.05. Using the globe's limit for both left every
country on the flat map barely zoomed at all.

Framing is fitted by bisection on the measured projected extent rather than from
a bounding radius. That matters in both views for different reasons: a sphere
curves away from the camera, so a large country projects smaller than its angular
radius implies (Russia framed at 0.42 of the viewport under the analytic fit,
0.70 under this one), and an elongated country's extent is not its radius either.

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
