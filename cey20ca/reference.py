"""Reference data: the 20 actual ERA5-Land years 2006-2025 at the 711 dataset points.

The published metrics were computed against ERA5-Land hourly 2 m temperature and dew point stored as
int16 hundredths of a degree Celsius. This module rebuilds that input with the same rules from the
Copernicus Climate Data Store (CDS). The CDS re-encodes area requests, so a few packed values differ by
0.01 degC from those of the global files (see README, "Expected result"):

1. Each dataset point (a 0.25 degree grid point) takes the nearest native ERA5-Land 0.1 degree grid point:
   row = rint((90 - lat) / 0.1), column = rint((lon mod 360) / 0.1) (NumPy round-half-to-even).
2. Values are packed in float32 arithmetic: rint((float32(K) - float32(273.15)) * float32(100)) as int16.
   On a 48-hour global sample of January 2010 this reproduced the packed values used for the dataset at
   all 354,592 land points; the same formula in float64 differs for about 0.08 % of the values.
3. February uses its first 672 hours (1-28 February); 29 February of leap years is not used.

Downloads need a CDS account and ``~/.cdsapirc`` (see https://cds.climate.copernicus.eu/how-to-api).
One request per year and month for a small California box (two variables, about 30-60 MB each).
"""

from __future__ import annotations

import calendar
import hashlib
from pathlib import Path

import numpy as np

YEARS = tuple(range(2006, 2026))
DATASET = "reanalysis-era5-land"
VARIABLES = {"t2m": "2m_temperature", "d2m": "2m_dewpoint_temperature"}
GRID_STEP = 0.1


def native_rows_cols(lat: np.ndarray, lon: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Nearest native ERA5-Land 0.1 degree row (from 90N) and column (from 0E) of each point."""
    lat = np.asarray(lat, dtype=np.float64)
    lon = np.asarray(lon, dtype=np.float64)
    rows = np.rint((90.0 - lat) / GRID_STEP).astype(np.int64)
    cols = np.mod(np.rint(np.mod(lon, 360.0) / GRID_STEP).astype(np.int64), 3600)
    return rows, cols


def request_area(rows: np.ndarray, cols: np.ndarray) -> list[float]:
    """CDS area [north, west, south, east] (degrees, west/east in -180..180) covering the native points."""
    lon = cols * GRID_STEP
    lon = np.where(lon >= 180.0, lon - 360.0, lon)
    return [round(float(90.0 - rows.min() * GRID_STEP), 1), round(float(lon.min()), 1),
            round(float(90.0 - rows.max() * GRID_STEP), 1), round(float(lon.max()), 1)]


def pack_celsius_x100(kelvin: np.ndarray) -> np.ndarray:
    """int16 hundredths of degC, float32 arithmetic as in the dataset's source packing; NaN -> error."""
    values = np.asarray(kelvin, dtype=np.float32)
    if not np.isfinite(values).all():
        raise ValueError("non-finite reference value at a dataset point")
    scaled = (values - np.float32(273.15)) * np.float32(100.0)
    return np.rint(np.clip(scaled, -32767, 32767)).astype(np.int16)


def reference_hours(year: int, month: int) -> int:
    """Hours of the month used for the reference (February: 672, also in leap years)."""
    return (28 if month == 2 else calendar.monthrange(year, month)[1]) * 24


def cache_path(cache: Path, year: int, month: int) -> Path:
    return Path(cache) / f"era5land_t2m_d2m_{year}_{month:02d}.npz"


def download_month(year: int, month: int, area: list[float], target: Path) -> None:
    import cdsapi

    request = {"variable": list(VARIABLES.values()), "year": str(year), "month": f"{month:02d}",
               "day": [f"{d:02d}" for d in range(1, calendar.monthrange(year, month)[1] + 1)],
               "time": [f"{h:02d}:00" for h in range(24)], "area": area,
               "data_format": "netcdf", "download_format": "unarchived"}
    cdsapi.Client(quiet=True).retrieve(DATASET, request, str(target))


def _index(values: np.ndarray, wanted: np.ndarray, label: str) -> np.ndarray:
    diff = np.abs(np.asarray(values, dtype=np.float64)[None, :] - wanted[:, None])
    idx = diff.argmin(axis=1)
    if np.any(diff[np.arange(wanted.size), idx] > 1e-4):
        raise ValueError(f"CDS file does not contain all native {label} values of the dataset points")
    return idx


def sample_month(path: Path, rows: np.ndarray, cols: np.ndarray, year: int, month: int) -> dict[str, np.ndarray]:
    """Packed reference values (hours, points) from one downloaded CDS NetCDF month."""
    import netCDF4

    with netCDF4.Dataset(path) as ds:
        time_name = "valid_time" if "valid_time" in ds.variables else "time"
        times = np.asarray(ds[time_name][:], dtype=np.int64)
        start = int(np.datetime64(f"{year:04d}-{month:02d}-01T00:00:00", "s").astype(np.int64))
        hours = reference_hours(year, month)
        expected = start + 3600 * np.arange(calendar.monthrange(year, month)[1] * 24, dtype=np.int64)
        if not np.array_equal(times, expected):
            raise ValueError(f"{path.name}: time axis is not the complete hourly month")
        lat_idx = _index(ds["latitude"][:], 90.0 - rows * GRID_STEP, "latitude")
        lon_values = np.mod(np.asarray(ds["longitude"][:], dtype=np.float64), 360.0)
        lon_idx = _index(lon_values, cols * GRID_STEP, "longitude")
        out = {}
        for short in VARIABLES:
            var = ds[short]
            var.set_auto_maskandscale(True)
            data = np.ma.filled(var[:hours, :, :], np.nan)
            out[short] = pack_celsius_x100(data[:, lat_idx, lon_idx])
        return out


def build_cache(lat: np.ndarray, lon: np.ndarray, cache: Path, downloads: Path, years=YEARS,
                workers: int = 4) -> None:
    """Download (if missing) and pack every year-month; existing cache files are kept.

    ``workers`` CDS requests run at the same time (each request waits in the CDS queue for minutes).
    A download is written to a temporary name and renamed when complete, so an interrupted run resumes.
    Months that fail (for example a rejected CDS request) are listed at the end; a second run retries them.
    """
    from concurrent.futures import ThreadPoolExecutor

    rows, cols = native_rows_cols(lat, lon)
    area = request_area(rows, cols)
    Path(cache).mkdir(parents=True, exist_ok=True)
    Path(downloads).mkdir(parents=True, exist_ok=True)

    def one(year: int, month: int) -> None:
        target = cache_path(cache, year, month)
        if target.exists():
            return
        raw = Path(downloads) / f"era5land_t2m_d2m_{year}_{month:02d}.nc"
        if not raw.exists():
            print(f"CDS request {year}-{month:02d} area {area}", flush=True)
            part = raw.with_suffix(".nc.part")
            download_month(year, month, area, part)
            part.replace(raw)
        packed = sample_month(raw, rows, cols, year, month)
        digest = hashlib.sha256(raw.read_bytes()).hexdigest()
        part = target.with_suffix(".npz.part")
        with open(part, "wb") as handle:  # a file handle keeps numpy from appending ".npz"
            np.savez_compressed(handle, t2m=packed["t2m"], d2m=packed["d2m"], source_sha256=digest)
        part.replace(target)
        print(f"packed {year}-{month:02d}", flush=True)

    with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
        futures = {(y, m): pool.submit(one, y, m) for y in years for m in range(1, 13)}
    failed = {k: f.exception() for k, f in futures.items() if f.exception() is not None}
    for (year, month), error in failed.items():
        print(f"failed {year}-{month:02d}: {type(error).__name__}: {error}", flush=True)
    if failed:
        raise RuntimeError(f"{len(failed)} months failed; run the command again to retry them")


def load_month(cache: Path, month: int, years=YEARS) -> dict[str, np.ndarray]:
    """Reference (years, hours, points) in degC, float32 = int16 * 0.01 as in the original reader."""
    t, td = [], []
    for year in years:
        with np.load(cache_path(cache, year, month)) as data:
            t.append(np.asarray(data["t2m"], dtype=np.float32) * 0.01)
            td.append(np.asarray(data["d2m"], dtype=np.float32) * 0.01)
    return {"temperature_c": np.stack(t), "dewpoint_c": np.stack(td)}
