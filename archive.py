"""Shared pieces of the frame archive: downsampling, the Nino 3.4 index, and the
on-disk layout. Imported by extract_sst_anomaly.py (daily) and
backfill_archive.py (one-time historical fill), so both write the same thing.

Layout:
    archive/index.json           manifest: grid geometry, the decode LUT, and the
                                 frame list. One shared LUT and grid rather than a
                                 sidecar per frame.
    archive/frames/<date>.png    0.5 deg quantized anomaly, same scheme as the
                                 full-resolution daily PNG.
    archive/nino34.json          the Nino 3.4 series -- deliberately a separate
                                 file, since it is a few KB against the frames'
                                 megabytes and the sparkline should not wait on them.

Frames are 0.5 deg rather than the native 0.25: 180 days costs 15.8 MB that way
against 45 MB at native, in a repo that already carries ~7 MB of district files
and takes a data commit every day. Playback does not need native resolution.
"""

import json
import os
from datetime import datetime, timezone

import numpy as np

ARCHIVE_DIR = "archive"
FRAMES_DIR = "frames"
WINDOW_DAYS = 180
ARCHIVE_STEP_DEG = 0.5
SMOOTHING_DAYS = 7        # see the README: 7-day kills the weather wiggle at ~zero lag

# The Nino 3.4 region, as NOAA defines it.
NINO34 = {"lat_min": -5.0, "lat_max": 5.0, "lon_min": -170.0, "lon_max": -120.0}


def nino34_index(values, grid_meta):
    """cos(lat)-weighted mean anomaly over the Nino 3.4 box.

    This mirrors averageBox() in index.html step for step -- same box, same
    row/column selection, same weighting -- because the archived series and the
    live readout have to be the same quantity or the sparkline's last point
    would disagree with the number printed beside it.

    `values` is north-up rows, ascending -180..180 columns: the same orientation
    the PNG is written in and the browser decodes to.
    """
    step = grid_meta["step_deg"]
    h, w = values.shape

    def row_for(lat):  return int(round((grid_meta["lat_max"] - lat) / step))
    def col_for(lon):  return int(round((lon - grid_meta["lon_min"]) / step))

    r0 = max(0, row_for(NINO34["lat_max"]));  r1 = min(h - 1, row_for(NINO34["lat_min"]))
    c0 = max(0, col_for(NINO34["lon_min"]));  c1 = min(w - 1, col_for(NINO34["lon_max"]))

    sub = values[r0:r1 + 1, c0:c1 + 1]
    lats = grid_meta["lat_max"] - np.arange(r0, r1 + 1) * step
    weights = np.repeat(np.cos(np.radians(lats))[:, None], sub.shape[1], axis=1)
    ok = np.isfinite(sub)
    if not ok.any():
        return None, 0
    return float((sub[ok] * weights[ok]).sum() / weights[ok].sum()), int(ok.sum())


def downsample(values, factor):
    """Block-mean a grid down by an integer factor, ignoring land.

    Averaging rather than decimating, so the archived frame keeps the 0.5 deg
    cell centres (89.875 and 89.625 average to 89.75) and does not alias. A block
    that is entirely land stays land.
    """
    h, w = values.shape
    h2, w2 = h // factor, w // factor
    blocks = values[:h2 * factor, :w2 * factor].reshape(h2, factor, w2, factor)
    ok = np.isfinite(blocks)
    total = np.where(ok, blocks, 0.0).sum(axis=(1, 3))
    count = ok.sum(axis=(1, 3))
    with np.errstate(invalid="ignore"):
        return np.where(count > 0, total / np.maximum(count, 1), np.nan)


def archive_grid_meta(native_meta, factor):
    """Grid geometry of the downsampled frame, derived from the native one."""
    g = native_meta["grid"]
    step = g["step_deg"] * factor
    return {
        "width": g["width"] // factor,
        "height": g["height"] // factor,
        # block-mean moves the edge cell centre inward by half the old step
        "lat_max": round(g["lat_max"] - g["step_deg"] * (factor - 1) / 2, 6),
        "lat_min": round(g["lat_min"] + g["step_deg"] * (factor - 1) / 2, 6),
        "lon_min": round(g["lon_min"] + g["step_deg"] * (factor - 1) / 2, 6),
        "lon_max": round(g["lon_max"] - g["step_deg"] * (factor - 1) / 2, 6),
        "step_deg": round(step, 6),
        "row_order": "north_to_south",
    }


# ---------------------------------------------------------------------------
# on-disk manifest
# ---------------------------------------------------------------------------
def _paths(root):
    return (os.path.join(root, "index.json"),
            os.path.join(root, "nino34.json"),
            os.path.join(root, FRAMES_DIR))


def load_index(root=ARCHIVE_DIR):
    path, _, _ = _paths(root)
    if os.path.exists(path):
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)
    return None


def load_series(root=ARCHIVE_DIR):
    _, path, _ = _paths(root)
    if os.path.exists(path):
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)
    return None


def trailing_mean(series, window):
    """Trailing mean over `window` days, by date rather than by position, so a
    gap in the archive shortens the window instead of silently reaching further
    back in time than it claims to."""
    out = []
    for i, row in enumerate(series):
        day = datetime.strptime(row["d"], "%Y-%m-%d").date()
        vals = []
        for j in range(i, -1, -1):
            other = datetime.strptime(series[j]["d"], "%Y-%m-%d").date()
            if (day - other).days >= window:
                break
            if series[j]["v"] is not None:
                vals.append(series[j]["v"])
        out.append(round(sum(vals) / len(vals), 4) if vals else None)
    return out


def write_archive(root, frames, series_rows, native_meta, factor, window_days=WINDOW_DAYS):
    """Write the manifest and the series, pruned to the window.

    `frames` is [(date_str, levels_uint8)] to add; existing frames on disk are
    kept unless they fall outside the window.
    """
    from PIL import Image
    index_path, series_path, frames_dir = _paths(root)
    os.makedirs(frames_dir, exist_ok=True)

    for date_str, levels in frames:
        Image.fromarray(levels, mode="L").save(
            os.path.join(frames_dir, f"{date_str}.png"), format="PNG", optimize=True)

    # everything actually on disk, so a half-finished run self-heals next time
    dates = sorted(f[:-4] for f in os.listdir(frames_dir) if f.endswith(".png"))
    keep = set(dates[-window_days:])
    for stale in set(dates) - keep:
        os.remove(os.path.join(frames_dir, f"{stale}.png"))
    dates = sorted(keep)

    index = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "window_days": window_days,
        "resolution_deg": round(native_meta["grid"]["step_deg"] * factor, 6),
        "source": native_meta.get("source"),
        "grid": archive_grid_meta(native_meta, factor),
        "encoding": {
            "format": "png-l8-lut",
            "nodata_level": native_meta["encoding"]["nodata_level"],
            "lut": native_meta["encoding"]["lut"],
        },
        "frames": [{"date": d, "file": f"{FRAMES_DIR}/{d}.png"} for d in dates],
    }
    with open(index_path, "w", encoding="utf-8") as fh:
        json.dump(index, fh, separators=(",", ":"))

    # Series is kept to the same window, and merged by date so a re-run or an
    # overlapping backfill updates in place rather than duplicating.
    existing = {r["d"]: r for r in (load_series(root) or {}).get("series", [])}
    for row in series_rows:
        existing[row["d"]] = row
    merged = [existing[d] for d in sorted(existing) if d in keep or not keep]
    merged = merged[-window_days:]
    smoothed = trailing_mean(merged, SMOOTHING_DAYS)
    for row, s in zip(merged, smoothed):
        row["s"] = s

    with open(series_path, "w", encoding="utf-8") as fh:
        json.dump({
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "region": NINO34,
            "smoothing_days": SMOOTHING_DAYS,
            "note": "v = daily index, s = trailing mean over smoothing_days",
            "series": merged,
        }, fh, separators=(",", ":"))

    return index, merged
