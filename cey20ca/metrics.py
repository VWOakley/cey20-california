"""Descriptive metrics of the CEY20 v2 California reference year against the 20 actual years.

These functions reproduce ``point_metrics.csv`` of the Zenodo dataset (doi:10.5281/zenodo.24016356)
with the same definitions and the same numerical types:

* Relative humidity: water-relative Magnus form of the CEY20 builder,
  e_s(T) = 6.112 * exp(17.67 T / (T + 243.5)) hPa, RH = clip(100 e_s(Td) / e_s(T), 0, 100), in float32.
* Exposure counts per month: hours with RH >= 80 %, RH >= 90 % and T - Td <= 2 degC.
* Annual: H_sel = count in the reference year; H_ref = mean over the 20 years of the annual count;
  abs_rel_error = |H_sel - H_ref| / H_ref; within_5pct = abs_rel_error <= 0.05, or H_sel == 0 where
  H_ref == 0 (then abs_rel_error is "N/A").
* Monthly: d = monthly mean of the reference year minus the 20-year mean of the monthly means, for
  temperature (degC) and RH (percentage points); within = |d| <= 0.5 degC or <= 2 percentage points.
* A point is UNKNOWN for annual metrics if any product value (T, Td, p) or any reference value (T, Td)
  is not finite, and for a month if any value of T or Td of that month is not finite. Missing values
  are never counted as zero.

Arrays: product values (hours, points) float32; reference values (years, hours, points) float32.
"""

from __future__ import annotations

import numpy as np

ANNUAL_REL_TOL = 0.05
MONTHLY_T_TOL_C = 0.5
MONTHLY_RH_TOL_PCT = 2.0
COUNT_KEYS = ("hours_rh_ge_80", "hours_rh_ge_90", "hours_dewpoint_depression_le_2c")
UNKNOWN = "UNKNOWN"


def saturation_vapor_pressure_hpa(temperature_c: np.ndarray) -> np.ndarray:
    return 6.112 * np.exp((17.67 * temperature_c) / (temperature_c + 243.5))


def relative_humidity(temperature_c: np.ndarray, dewpoint_c: np.ndarray) -> np.ndarray:
    saturation = saturation_vapor_pressure_hpa(temperature_c)
    vapor = saturation_vapor_pressure_hpa(dewpoint_c)
    return np.clip(vapor / saturation * 100.0, 0.0, 100.0)


def counts(t: np.ndarray, td: np.ndarray) -> dict[str, np.ndarray]:
    """Exposure counts over the hour axis (axis 1): arrays (n, hours, points) -> (n, points)."""
    rh = relative_humidity(t, td).astype(np.float32)
    spread = t - td
    return {"hours_rh_ge_80": np.count_nonzero(rh >= 80.0, axis=1).astype(np.float64),
            "hours_rh_ge_90": np.count_nonzero(rh >= 90.0, axis=1).astype(np.float64),
            "hours_dewpoint_depression_le_2c": np.count_nonzero(spread <= 2.0, axis=1).astype(np.float64)}


def complete_points(*arrays: np.ndarray) -> np.ndarray:
    """(points,) bool: every value finite in every array; the point axis is the last axis."""
    ok = None
    for a in arrays:
        f = np.isfinite(a).reshape(-1, a.shape[-1]).all(axis=0)
        ok = f if ok is None else ok & f
    return ok


def annual_metrics(product_t: np.ndarray, product_td: np.ndarray, years_counts: dict,
                   complete: np.ndarray) -> list[dict]:
    """product_*: (8760, points); years_counts[k]: (months, years, points) monthly counts of the 20 years."""
    sel = counts(product_t[None, :, :], product_td[None, :, :])
    rows = []
    for k in COUNT_KEYS:
        h_ref = years_counts[k].sum(axis=0).mean(axis=0)
        for p in range(product_t.shape[1]):
            if not complete[p]:
                rows.append({"point": p, "metric": k, "h_selected": UNKNOWN, "h_ref": UNKNOWN,
                             "signed_h": UNKNOWN, "abs_rel_error": UNKNOWN, "within_5pct": UNKNOWN})
                continue
            h_sel = float(sel[k][0, p])
            signed = h_sel - float(h_ref[p])
            rel = abs(signed) / float(h_ref[p]) if h_ref[p] > 0 else None
            ok = (rel <= ANNUAL_REL_TOL) if rel is not None else (h_sel == 0)
            rows.append({"point": p, "metric": k, "h_selected": h_sel, "h_ref": float(h_ref[p]),
                         "signed_h": signed, "abs_rel_error": rel if rel is not None else "N/A",
                         "within_5pct": bool(ok)})
    return rows


def monthly_metrics(month: int, product_t: np.ndarray, product_td: np.ndarray,
                    years_t: np.ndarray, years_td: np.ndarray) -> list[dict]:
    """product_*: (hours, points) of the month; years_*: (years, hours, points)."""
    complete = complete_points(product_t, product_td, years_t, years_td)
    rh_p = relative_humidity(product_t, product_td)
    rh_y = relative_humidity(years_t, years_td)
    ref_t = np.mean(np.mean(years_t, axis=1), axis=0)
    ref_rh = np.mean(np.mean(rh_y, axis=1), axis=0)
    rows = []
    for p in range(product_t.shape[1]):
        if not complete[p]:
            rows.append({"point": p, "month": month, "d_temperature_c": UNKNOWN, "d_rh_pct": UNKNOWN,
                         "within_t": UNKNOWN, "within_rh": UNKNOWN})
            continue
        d_t = float(np.mean(product_t[:, p]) - ref_t[p])
        d_rh = float(np.mean(rh_p[:, p]) - ref_rh[p])
        rows.append({"point": p, "month": month, "d_temperature_c": d_t, "d_rh_pct": d_rh,
                     "within_t": abs(d_t) <= MONTHLY_T_TOL_C, "within_rh": abs(d_rh) <= MONTHLY_RH_TOL_PCT})
    return rows


def point_monthly_summary(mon: list[dict]) -> dict:
    if any(r["d_temperature_c"] == UNKNOWN for r in mon):
        return {"max_abs_d_temperature_c": UNKNOWN, "max_abs_d_rh_pct": UNKNOWN,
                "months_within_t": UNKNOWN, "months_within_rh": UNKNOWN}
    return {"max_abs_d_temperature_c": max(abs(r["d_temperature_c"]) for r in mon),
            "max_abs_d_rh_pct": max(abs(r["d_rh_pct"]) for r in mon),
            "months_within_t": sum(r["within_t"] for r in mon),
            "months_within_rh": sum(r["within_rh"] for r in mon)}
