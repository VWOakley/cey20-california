"""Unit tests of the metric definitions, the point mapping and the reference packing."""

from __future__ import annotations

import numpy as np
import pytest

from cey20ca import compute, metrics as M, reference as R


def test_relative_humidity_magnus_and_clip():
    assert M.relative_humidity(np.float32(20.0), np.float32(20.0)) == pytest.approx(100.0)
    assert M.relative_humidity(np.float32(10.0), np.float32(12.0)) == 100.0  # Td > T is clipped
    rh = M.relative_humidity(np.float32(25.0), np.float32(15.0))
    expected = 6.112 * np.exp(17.67 * 15 / (15 + 243.5)) / (6.112 * np.exp(17.67 * 25 / (25 + 243.5))) * 100
    assert rh == pytest.approx(expected, rel=1e-6)


def test_counts_thresholds_are_inclusive():
    t = np.array([[[10.0], [10.0], [10.0], [10.0]]], dtype=np.float32)
    td = np.array([[[10.0], [8.0], [9.5], [-10.0]]], dtype=np.float32)
    c = M.counts(t, td)
    assert c["hours_dewpoint_depression_le_2c"][0, 0] == 3  # spreads 0, 2, 0.5 are <= 2
    assert c["hours_rh_ge_90"][0, 0] == 2 and c["hours_rh_ge_80"][0, 0] == 3


def test_annual_metrics_zero_reference_rule():
    product_t = np.full((10, 2), 30.0, dtype=np.float32)
    product_td = np.full((10, 2), -20.0, dtype=np.float32)       # dry: zero hours
    zero = np.zeros((1, 20, 2))
    years_counts = {k: zero.copy() for k in M.COUNT_KEYS}
    years_counts["hours_rh_ge_80"][0, :, 1] = 10.0                  # point 1: reference 10 h, selected 0 h
    rows = M.annual_metrics(product_t, product_td, years_counts, np.array([True, True]))
    rh80 = [r for r in rows if r["metric"] == "hours_rh_ge_80"]
    assert rh80[0]["abs_rel_error"] == "N/A" and rh80[0]["within_5pct"] is True
    assert rh80[1]["abs_rel_error"] == 1.0 and rh80[1]["within_5pct"] is False and rh80[1]["signed_h"] == -10.0


def test_incomplete_points_are_unknown_not_zero():
    product_t = np.full((24, 1), 10.0, dtype=np.float32)
    product_t[3, 0] = np.nan
    product_td = np.full((24, 1), 5.0, dtype=np.float32)
    years = np.full((20, 24, 1), 10.0, dtype=np.float32)
    rows = M.monthly_metrics(1, product_t, product_td, years, years - 5)
    assert rows[0]["d_temperature_c"] == M.UNKNOWN
    assert M.point_monthly_summary(rows)["months_within_t"] == M.UNKNOWN


def test_native_mapping_rounds_half_to_even():
    rows, cols = R.native_rows_cols(np.array([37.75, 37.25]), np.array([-122.25, -122.0]))
    assert rows.tolist() == [522, 528]          # 522.5 -> 522, 527.5 -> 528
    assert cols.tolist() == [2378, 2380]        # 2377.5 -> 2378, (360 - 122.0) / 0.1 = 2380


def test_packing_uses_float32_arithmetic():
    kelvin = np.array([250.015], dtype=np.float32)   # float32 and float64 arithmetic round differently
    assert R.pack_celsius_x100(kelvin)[0] == -2313
    assert int(np.rint((kelvin.astype(np.float64) - 273.15) * 100)[0]) == -2314
    with pytest.raises(ValueError):
        R.pack_celsius_x100(np.array([np.nan], dtype=np.float32))


def test_reference_hours_drop_29_february():
    assert R.reference_hours(2008, 2) == 672 and R.reference_hours(2007, 2) == 672
    assert R.reference_hours(2010, 1) == 744 and sum(R.reference_hours(2001, m) for m in range(1, 13)) == 8760


def _cds_like_file(path, year, month, lats, lons, kelvin):
    import netCDF4

    hours = R.calendar.monthrange(year, month)[1] * 24
    with netCDF4.Dataset(path, "w") as ds:
        ds.createDimension("valid_time", hours)
        ds.createDimension("latitude", lats.size)
        ds.createDimension("longitude", lons.size)
        start = int(np.datetime64(f"{year:04d}-{month:02d}-01T00:00:00", "s").astype(np.int64))
        ds.createVariable("valid_time", "i8", ("valid_time",))[:] = start + 3600 * np.arange(hours)
        ds.createVariable("latitude", "f8", ("latitude",))[:] = lats
        ds.createVariable("longitude", "f8", ("longitude",))[:] = lons
        for name, offset in (("t2m", 0.0), ("d2m", -5.0)):
            v = ds.createVariable(name, "f4", ("valid_time", "latitude", "longitude"))
            v[:] = np.broadcast_to(kelvin + offset, (hours, lats.size, lons.size))


def test_sample_month_from_cds_layout(tmp_path):
    lats = np.round(np.arange(38.0, 37.0, -0.1), 1)
    lons = np.round(np.arange(-122.5, -121.5, 0.1), 1)
    path = tmp_path / "m.nc"
    _cds_like_file(path, 2008, 2, lats, lons, np.float32(290.0))
    rows, cols = R.native_rows_cols(np.array([37.75]), np.array([-122.25]))
    packed = R.sample_month(path, rows, cols, 2008, 2)
    assert packed["t2m"].shape == (672, 1)
    assert packed["t2m"][0, 0] == R.pack_celsius_x100(np.array([290.0], dtype=np.float32))[0]
    assert packed["d2m"][0, 0] == R.pack_celsius_x100(np.array([285.0], dtype=np.float32))[0]


def test_compare_reports_differences(tmp_path):
    rows = [{"point_index": 1, "x": 0.5, "flag": True, "e": "N/A"}]
    published = tmp_path / "p.csv"
    compute.write_csv(published, rows)
    assert compute.compare(published, rows) == []
    assert compute.compare(published, [{"point_index": 1, "x": 0.6, "flag": True, "e": "N/A"}])


def test_summarize_counts_points_and_largest_difference(tmp_path):
    rows = [{"point_index": 1, "x": 0.5, "flag": True}, {"point_index": 2, "x": 1.0, "flag": False}]
    published = tmp_path / "p.csv"
    compute.write_csv(published, rows)
    assert compute.summarize(published, rows) == []
    new = [{"point_index": 1, "x": 0.75, "flag": False}, {"point_index": 2, "x": 1.25, "flag": False}]
    assert compute.summarize(published, new) == ["x: 2 points, max |difference| 0.25", "flag: 1 points"]
