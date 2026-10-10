# Whale Movement Predictability

Building a model / formula to quantify how predictable whale movement is, using
environmental covariates (e.g. IMOS pulls, SRS ocean colour Chl-a as a proxy for
prey / productivity).

## Repo structure

```
data-raw/     # download scripts and original inputs (never commit big rasters)
data/         # processed, gridded outputs
R/            # functions (PDE, suitability, prey proxy)
analysis/     # notebooks / scripts
```

- Raw rasters / NetCDF files are git-ignored. Anything in `data-raw/` must be
  re-creatable from a download script so the IMOS pulls are reproducible.
- Reusable code goes in `R/` as functions; `analysis/` calls those functions.
- Pipeline: plan is an R `{targets}` pipeline (or plain scripts run in order).

## Setup (R / RStudio)

NetCDF reading on Windows needs system NetCDF libraries. Setup steps TBD.

## Data: `data-raw/Humpback_East_Australia_data.xlsx`

Two sheets.

**`Humpback`** (33,406 sighting records, OBIS 28,571 + GBIF 4,835; lat -45 to -9, lon 142 to 170)

| column | notes |
|---|---|
| source_database | OBIS or GBIF |
| record_id | unique per record |
| date, year, month | dates run 1770 to 2026; ~900 records have no date, ~760 no year |
| latitude, longitude | decimal degrees |
| individual_count | blank for ~14,400 records (count not recorded); max 544 |
| depth_m | mostly blank; where present it is 0 |
| basis_of_record | 6 types (e.g. HumanObservation) |
| dataset, institution | provenance |

**`Water`** (62,190 rows, monthly, 2-degree grid cells, 2002 to 2026; lat -45 to -11, lon 143 to 169)

| column | notes |
|---|---|
| year, month | |
| cell_centre_lat, cell_centre_lon | grid cell centre |
| sea_surface_temp_C | ~3% missing |
| chlorophyll_mg_per_m3 | Chl-a, prey proxy; ~23% missing (cloud / gaps) |
| salinity_psu_SODA_to2015 | SODA reanalysis, only available to 2015; ~44% missing |

Things to handle in the model: sightings before 2002 have no matching water
data, sightings are presence-only (no absences), and effort is uneven over time.

## Extra environmental data (`data/`, built by `data-raw/0*.py`)

| file | what | coverage |
|---|---|---|
| `soi_monthly.csv` | Southern Oscillation Index (NOAA CPC, standardized) | 1951 to 2026 |
| `bathymetry_cells.csv` | depth, depth SD, shelf fraction, distance to 200 m shelf edge per 2-degree cell (ETOPO) | static |
| `currents_cells.csv` | IMOS GSLA geostrophic currents (u, v, speed, EKE) per cell-month | 2002 to 2015 |

Currents notes: sampled every 15th day (about 2 per month; 2002 uses every 5th
day), so `n_days` is small and `eke` is unreliable where `n_days` is 1. Years after
2015 were not downloaded (chlorophyll only goes to 2015); extend with
`python data-raw/03_download_currents.py 2002 2025` (slow, about 8 min per year).
Run from the repo root; needs pandas, requests, xarray, netCDF4.

## Modelling table: `data/cell_month_2002_2015.csv`

Built by `data-raw/04_build_cell_month_table.py`. One row per 2-degree cell and month
(35,420 rows, 2002 to 2015). Contains the Water columns, `sst_anom_C` (vs. the cell's
2002-2015 monthly mean), `chl_log10`, previous-month `sst_lag1_C` / `chl_log10_lag1`,
SOI, bathymetry, currents, and the whale columns `n_sightings`, `n_individuals`
(only records that gave a count; see `n_with_count`), `present` and effort proxies
`effort_year`, `effort_year_month`.

Zeros in `n_sightings` mean "no sighting recorded", not confirmed absence. Effort is
very uneven (e.g. 3,611 sightings in 2008 vs ~150 to 650 in most years), so include
effort in models.
