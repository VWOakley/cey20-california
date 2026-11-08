"""Recompute point_metrics.csv from the published NetCDF and the reference cache, and compare."""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np

from . import metrics as M
from . import reference as R

MONTH_HOURS = [R.reference_hours(2001, m) for m in range(1, 13)]  # 365-day reference year


def _as_float(values) -> np.ndarray:
    return np.ma.filled(np.ma.asarray(values).astype(np.float32), np.nan)


def load_product(nc_path: Path) -> dict[str, np.ndarray]:
    import netCDF4

    with netCDF4.Dataset(nc_path) as ds:
        out = {k: _as_float(ds[k][:]) for k in ("temperature_c", "dewpoint_c", "surface_pressure_pa")}
        for k in ("latitude", "longitude"):
            out[k] = np.asarray(ds[k][:], dtype=np.float64)
        out["product_point_index"] = np.asarray(ds["product_point_index"][:])
        out["center_in_california"] = np.asarray(ds["center_in_california"][:]).astype(bool)
        out["cell_fraction_in_california"] = np.asarray(ds["cell_fraction_in_california"][:], dtype=np.float64)
    return out


def point_metrics(nc_path: Path, cache: Path) -> list[dict]:
    prod = load_product(nc_path)
    product_complete = M.complete_points(prod["temperature_c"], prod["dewpoint_c"], prod["surface_pressure_pa"])
    n = prod["temperature_c"].shape[1]
    reference_complete = np.ones(n, dtype=bool)
    years_counts = {k: [] for k in M.COUNT_KEYS}
    monthly_rows, start = [], 0
    for m in range(1, 13):
        md = R.load_month(cache, m)
        reference_complete &= M.complete_points(md["temperature_c"], md["dewpoint_c"])
        c = M.counts(md["temperature_c"], md["dewpoint_c"])
        for k in M.COUNT_KEYS:
            years_counts[k].append(c[k])
        stop = start + MONTH_HOURS[m - 1]
        monthly_rows += M.monthly_metrics(m, prod["temperature_c"][start:stop], prod["dewpoint_c"][start:stop],
                                          md["temperature_c"], md["dewpoint_c"])
        start = stop
    years_counts = {k: np.stack(v) for k, v in years_counts.items()}
    annual_rows = M.annual_metrics(prod["temperature_c"], prod["dewpoint_c"], years_counts,
                                   product_complete & reference_complete)
    rows = []
    for p in range(n):
        lon = float(prod["longitude"][p])
        row = {"point_index": int(prod["product_point_index"][p]), "latitude": float(prod["latitude"][p]),
               "longitude": lon - 360.0 if lon > 180.0 else lon,
               "center_in_california": bool(prod["center_in_california"][p]),
               "cell_fraction_in_california": round(float(prod["cell_fraction_in_california"][p]), 6),
               "product_complete": bool(product_complete[p]), "reference_complete": bool(reference_complete[p])}
        for a in (r for r in annual_rows if r["point"] == p):
            for f in ("h_selected", "h_ref", "signed_h", "abs_rel_error", "within_5pct"):
                row[f"{a['metric']}_{f}"] = a[f]
        mon = [r for r in monthly_rows if r["point"] == p]
        row.update(M.point_monthly_summary(mon))
        for r in mon:
            row[f"d_temperature_c_m{r['month']:02d}"] = r["d_temperature_c"]
            row[f"d_rh_pct_m{r['month']:02d}"] = r["d_rh_pct"]
        rows.append(row)
    return rows


def write_csv(path: Path, rows: list[dict]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _parse(value: str):
    if value in ("True", "False"):
        return value == "True"
    try:
        return float(value)
    except ValueError:
        return value


def summarize(published_csv: Path, rows: list[dict]) -> list[str]:
    """Per column: number of points that differ and the largest absolute difference (flags: count only)."""
    with open(published_csv, encoding="utf-8") as handle:
        published = list(csv.DictReader(handle))
    lines = []
    for key in published[0]:
        if key == "cell_fraction_in_california":
            continue
        n, worst = 0, 0.0
        for pub, new in zip(published, rows):
            a, b = _parse(pub[key]), new.get(key)
            if isinstance(a, bool) or isinstance(b, bool) or not isinstance(a, float):
                n += a != b and str(b) != pub[key]
            elif a != b:
                n += 1
                worst = max(worst, abs(a - float(b)))
        if n:
            lines.append(f"{key}: {n} points, max |difference| {worst:.3g}" if worst else f"{key}: {n} points")
    return lines


def compare(published_csv: Path, rows: list[dict], rel_tol: float = 0.0) -> list[str]:
    """Differences between the published point_metrics.csv and recomputed rows (empty list: identical).

    The cell fraction comes from the boundary geometry and is not recomputed here: the published CSV
    rounds the float64 value to 6 decimals, the NetCDF stores it as float32, so the two may differ by one
    unit of the 6th decimal. Every metric is compared exactly unless ``rel_tol`` > 0.
    """
    with open(published_csv, encoding="utf-8") as handle:
        published = list(csv.DictReader(handle))
    problems = []
    if len(published) != len(rows):
        return [f"row count {len(published)} != {len(rows)}"]
    if list(published[0]) != list(rows[0]):
        problems.append("column order or names differ")
    for i, (pub, new) in enumerate(zip(published, rows)):
        for key, raw in pub.items():
            a, b = _parse(raw), new.get(key)
            if key == "cell_fraction_in_california":
                if not (isinstance(a, float) and abs(a - float(b)) <= 1e-6 + 1e-12):
                    problems.append(f"row {i} {key}: published {raw} recomputed {b!r}")
            elif isinstance(a, float) and isinstance(b, (int, float)) and not isinstance(b, bool):
                if not (a == b or (rel_tol and math.isclose(a, b, rel_tol=rel_tol, abs_tol=0.0))):
                    problems.append(f"row {i} {key}: published {raw} recomputed {b!r}")
            elif a != b and str(b) != raw:
                problems.append(f"row {i} {key}: published {raw} recomputed {b!r}")
    return problems
