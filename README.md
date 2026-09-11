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
