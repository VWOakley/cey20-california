# CEY20 v2 California exposure reference year: tools

Dataset: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24016356.svg)](https://doi.org/10.5281/zenodo.24016356)
Code: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24073640.svg)](https://doi.org/10.5281/zenodo.24073640)

Code to download, verify and use the dataset **CEY20 v2 exposure reference year: California
(0.25 degree), dataset v1.0.0** ([doi:10.5281/zenodo.24016356](https://doi.org/10.5281/zenodo.24016356)),
and to recompute its descriptive metrics (`point_metrics.csv`) independently from ERA5-Land.

The data stay on Zenodo; this repository contains no data files.

## What the dataset is

An hourly 365-day reference year of temperature, dew point, surface pressure and relative humidity for
the 711 ERA5-Land land points (0.25 degree) whose cells overlap California (663 with the cell centre
inside the state). Each calendar month is one actual source year (2006-2025) chosen per point by the
CEY20 selection; month interfaces of T, Td and p are smoothed over 6 h on each side.

It is an exposure reference year for moisture and thermal screening, **not an energy TMY**: it contains
no radiation and no wind. It is reanalysis data, not station data; coastal, urban and mountain
microclimates are not resolved. A county file is one grid point near the county centroid, not a county
average. See the dataset's `README.md` for the full description and limitations.

## Quick start

```bash
python -m pip install -r requirements.txt
python -m cey20ca download            # Zenodo files into ./data, MD5 and SHA-256 verified
jupyter notebook notebooks/getting-started.ipynb
```

The notebook is also on Kaggle:
[CEY20 California Exposure Year: Getting Started](https://www.kaggle.com/code/okleywadim/cey20-california-exposure-year-getting-started).

## Recompute the metrics from ERA5-Land

```bash
python -m cey20ca reference           # 240 CDS requests (2006-2025, monthly, California box)
python -m cey20ca metrics             # recompute point_metrics.csv and compare with the published file
```

`reference` needs a free account at the Copernicus Climate Data Store and an API key in `~/.cdsapirc`
([how to](https://cds.climate.copernicus.eu/how-to-api)); you must accept the ERA5-Land licence on the
CDS website. It downloads 2 m temperature and dew point (about 20 MB per month) and builds the
reference with the same rules as the dataset:

- each dataset point takes the nearest native ERA5-Land 0.1 degree grid point
  (`rint((90 - lat) / 0.1)`, `rint((lon mod 360) / 0.1)`, NumPy round-half-to-even);
- values are stored as int16 hundredths of a degree Celsius, computed in float32:
  `rint((float32(K) - float32(273.15)) * float32(100))`;
- February uses its first 672 hours; 29 February is not used.

The CDS queues each request; expect several hours for all 240 months. `--workers` (default 4) sets how
many requests run at once. If the CDS rejects requests (state `rejected`, for example while other
requests of the same account are running), use `--workers 2` or `1`. Failed months are listed at the
end, and running the command again downloads only the missing months.

`metrics` then applies the definitions in [`cey20ca/metrics.py`](cey20ca/metrics.py) (relative humidity
from the water-relative Magnus form of the CEY20 builder, exposure counts, annual and monthly
comparisons) and compares the result with the published `point_metrics.csv`, listing per column how
many points differ and by how much.

### Observed CDS result on 8 November 2026: small differences

The dataset was built from the global ERA5-Land GRIB files. The checked CDS area requests, in both
NetCDF and GRIB formats, returned re-encoded fields whose values differed from the global files by up
to 0.00023 K. After applying the dataset's packing rule, about 1 % of temperature values and 0.5 % of
dew-point values differed by 0.01 °C in one checked month.

Measured on 8 November 2026 with all 240 months from the CDS:

| Published column | Points that differ (of 711) | Largest difference |
|---|---|---|
| `*_h_ref`, `*_signed_h` (mean annual hours of the 20 years) | 394-425 per metric | 0.45 h per year |
| `*_abs_rel_error` | 394-425 per metric | 0.00144 (dimensionless; equivalent to 0.144 percentage points) |
| `d_temperature_c_mMM`, `max_abs_d_temperature_c` | up to all | 0.00016 °C |
| `d_rh_pct_mMM`, `max_abs_d_rh_pct` | up to all | 0.00044 percentage points |
| `months_within_t` | 1 (point index 141804) | 1 month (July: -0.50003 °C published, -0.49993 °C recomputed, limit 0.5 °C) |
| `*_within_5pct`, `*_h_selected`, `months_within_rh`, completeness flags | none | – |

In the full CDS run checked on 8 November 2026, `metrics` ended with exit code 1 because of these
differences. Recomputed from the packed values of the global files, the same code reproduces the
published `point_metrics.csv` for all 711 points, with `cell_fraction_in_california` compared to a
tolerance of 1e-6 because it is stored as float32 in the NetCDF; all other columns match exactly. Those
packed values are not distributed.

These metrics describe how the reference year compares with the 20 actual years of the same
reanalysis. They are not a validation against observations.

## Repository layout

| Path | Content |
|---|---|
| `cey20ca/download.py` | Zenodo download, MD5 (Zenodo) and SHA-256 (`manifest.json`) checks |
| `cey20ca/reference.py` | CDS download, point mapping and packing of the 20 reference years |
| `cey20ca/metrics.py` | Metric definitions |
| `cey20ca/compute.py` | Recomputation of `point_metrics.csv` and comparison |
| `notebooks/getting-started.ipynb` | Loading the NetCDF and a county file, two example plots |
| `tests/` | Unit tests of the definitions and the packing rule |

The CEY20 month selection itself is not part of this repository.
It is described at https://www.getclimetry.com/methodology#cey20 and in the CEY20 study package
([doi:10.5281/zenodo.24091411](https://doi.org/10.5281/zenodo.24091411)).

## Citation

Please cite the dataset:

> Okley, W. (2026): CEY20 v2 exposure reference year: California (0.25 degree), dataset v1.0.0.
> Zenodo. https://doi.org/10.5281/zenodo.24016356

The dataset is described in the data descriptor (preprint):

> Okley, W. (2026): CEY20 v2 exposure reference year for California: an hourly temperature and humidity
> dataset at 0.25° from ERA5-Land 2006–2025. Data descriptor, version 1.0. Zenodo.
> https://doi.org/10.5281/zenodo.24073751

To cite this code, use doi:10.5281/zenodo.24073640 (all versions; v1.0.0: doi:10.5281/zenodo.24073641).
`CITATION.cff` contains its citation metadata.

## Licences

- Code: MIT (see `LICENSE`).
- Dataset: CC BY 4.0 (on Zenodo). Contains modified Copernicus Climate Change Service information
  (2026). Neither the European Commission nor ECMWF is responsible for any use that may be made of the
  Copernicus information or data it contains.
- ERA5-Land downloaded with `reference` is subject to the CDS licence for that dataset
  (doi:10.24381/cds.e2161bac).
