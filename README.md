# sourcing-climate-monitor

A single-page 3D globe (Three.js) showing daily **sea surface temperature (SST)
anomaly** data — built as a cocoa-sourcing climate-risk view for a chocolate
manufacturer. The data is refreshed automatically each day by a GitHub Actions
workflow and served over GitHub Pages.

## How it works

- **`extract_sst_anomaly.py`** pulls the latest daily SST from NOAA's OISST v2
  High-Res dataset, subtracts the 1991-2020 daily climatology, downsamples to a
  ~1° grid, and writes **`sst_anomaly.json`** (lat/lon arrays + a 2D anomaly
  grid; `null` = land / no data).
- **`index.html`** fetches `./sst_anomaly.json` on load and renders it onto the
  globe. If the fetch fails (opened locally, or before the first data commit) it
  falls back to a manual file picker.
- **`.github/workflows/refresh.yml`** runs the script daily at 13:00 UTC and
  commits the refreshed `sst_anomaly.json` back to the repo when it changes.

## Running locally

```bash
pip install xarray netCDF4 pandas numpy requests
python extract_sst_anomaly.py --out sst_anomaly.json
```

Then serve the folder (a plain file open won't allow the `fetch`, so use a
server) and visit it:

```bash
python -m http.server 8000
# open http://localhost:8000/
```

The first run downloads the ~90 MB 1991-2020 climatology once and caches it as
`sst.day.ltm.1991-2020.nc` (git-ignored).
