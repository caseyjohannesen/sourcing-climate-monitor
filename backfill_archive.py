"""One-time historical fill of the frame archive.

Run this once to populate ~180 days of archive rather than waiting six months
for the daily job to accumulate them. After that, extract_sst_anomaly.py
--archive extends the window one day at a time.

    python backfill_archive.py [--days 180] [--archive archive]

Meant to run on a GitHub Actions runner (see the `backfill-archive` workflow):
PSL's THREDDS/download hosts are reachable from there, and the 0.25 deg
climatology is already in the Actions cache. It is not expected to work from a
laptop behind a network that cannot reach PSL.

Why it downloads whole yearly files rather than reading day slices over OPeNDAP:
180 separate OPeNDAP reads of a global field is ~720 MB of round trips and leans
hard on a public service, where the yearly file is one request. It also keeps the
frames honest -- downsampling is a block mean over the native grid, exactly as
the daily path does it, rather than the server-side decimation a strided OPeNDAP
request would give.
"""

import argparse
import os
import sys
from datetime import date, timedelta

import numpy as np
import pandas as pd
import xarray as xr

import archive as arch
import extract_sst_anomaly as ex

YEAR_FILE_URL = "https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2.highres/sst.day.mean.{year}.nc"


def ensure_year_file(year, cache_dir="."):
    path = os.path.join(cache_dir, f"sst.day.mean.{year}.nc")
    if os.path.exists(path):
        print(f"Using cached {path}")
        return path
    url = YEAR_FILE_URL.format(year=year)
    print(f"Downloading {url} (one request, a few hundred MB) ...")
    return ex.download_with_retries(url, path)


def frames_for_window(sst_by_year, clim_ds, days, end_date=None):
    """Yield (date_str, native_anomaly) for each day in the window, newest last."""
    end = end_date or max(
        pd.Timestamp(ds.time.isel(time=-1).values).date() for ds in sst_by_year.values())
    start = end - timedelta(days=days - 1)

    day = start
    while day <= end:
        ds = sst_by_year.get(day.year)
        if ds is None:
            day += timedelta(days=1)
            continue
        try:
            latest = ds.sst.sel(time=str(day))
        except KeyError:
            print(f"  {day}: no data in the source file, skipping", file=sys.stderr)
            day += timedelta(days=1)
            continue
        if "time" in latest.dims:
            latest = latest.isel(time=0)

        ts = pd.Timestamp(day)
        clim = ex.align_climatology(ex.match_climatology(clim_ds, ts), latest)
        yield day.isoformat(), (latest - clim)
        day += timedelta(days=1)


def orient(anomaly):
    """Native anomaly -> (north-up rows, -180..180 columns), matching the daily path."""
    lats = np.asarray(anomaly.lat.values, dtype=float)
    lons = np.asarray(anomaly.lon.values, dtype=float)
    grid = np.asarray(anomaly.values, dtype=float)

    lons = np.where(lons > 180, lons - 360, lons)
    order = np.argsort(lons)
    lons, grid = lons[order], grid[:, order]
    if lats[0] < lats[-1]:
        lats, grid = lats[::-1], grid[::-1, :]
    return grid, lats, lons


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--days", type=int, default=arch.WINDOW_DAYS)
    ap.add_argument("--archive", default=arch.ARCHIVE_DIR)
    ap.add_argument("--clim-file", help="(testing) local climatology instead of the cached download")
    ap.add_argument("--sst-file", action="append",
                    help="(testing) local SST NetCDF; repeatable, one per year")
    args = ap.parse_args()

    if args.sst_file:
        sst_by_year = {}
        for path in args.sst_file:
            ds = xr.open_dataset(path)
            sst_by_year[pd.Timestamp(ds.time.isel(time=0).values).year] = ds
    else:
        today = date.today()
        years = sorted({today.year, (today - timedelta(days=args.days)).year})
        sst_by_year = {y: xr.open_dataset(ensure_year_file(y)) for y in years}

    clim_path = args.clim_file or ex.ensure_local_climatology()
    clim_ds = xr.open_dataset(clim_path)

    frames, series, meta = [], [], None
    for date_str, anomaly in frames_for_window(sst_by_year, clim_ds, args.days):
        grid, lats, lons = orient(anomaly)
        levels, _, _ = ex.quantize(grid)

        if meta is None:                       # geometry and LUT are the same every day
            meta = ex.build_meta(grid, lats, lons, pd.Timestamp(date_str),
                                 float(abs(lats[0] - lats[1])), 0, 0, "sst_anomaly.png")

        # Index from the decoded levels, not the float grid, so a backfilled point
        # is the same number the daily path and the browser would produce.
        lut = np.array([np.nan if v is None else v for v in meta["encoding"]["lut"]], dtype="float64")
        decoded = lut[levels]

        value, cells = arch.nino34_index(decoded, meta["grid"])
        factor = int(round(arch.ARCHIVE_STEP_DEG / meta["grid"]["step_deg"]))
        frame_levels, _, _ = ex.quantize(arch.downsample(decoded, factor))

        frames.append((date_str, frame_levels))
        series.append({"d": date_str, "v": None if value is None else round(value, 4)})
        print(f"  {date_str}  Nino 3.4 {value:+.3f} C ({cells:,} cells)")

    if not frames:
        raise SystemExit("No frames produced — check the source files and the window.")

    index, merged = arch.write_archive(args.archive, frames, series, meta, factor, args.days)
    total = sum(os.path.getsize(os.path.join(args.archive, f["file"])) for f in index["frames"])
    print(f"\nArchive: {len(index['frames'])} frames at {index['resolution_deg']} deg "
          f"({total/1e6:.1f} MB), series {len(merged)} points "
          f"({merged[0]['d']} .. {merged[-1]['d']})")


if __name__ == "__main__":
    main()
