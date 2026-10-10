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

## Pipeline (run in order, from the repo root)

1. `analysis/03_download_southern_water.R` — merges the IMOS/AODN southern SST
   extension (downloaded by `analysis/fetch_sst_southern.py` into
   `data-raw/sst_southern_2003_2022.csv`) with the Water sheet →
   `data/water_extended.rds` (lat -67 to -11, 2003-2022 southern coverage).
2. `analysis/04_add_southern_sightings.R` — appends Southern Ocean humpback
   records (`data-raw/humpback_southern_ocean_obis.csv`, OBIS API pull,
   lat -65 to -45, lon 80-180) to the sightings → `data-raw/humpback_all.rds`.
3. `analysis/01_join_sightings_water.R` — joins sightings (2002+) to the
   nearest water cell per year-month (1.5° tolerance, recovers coastal
   records), generates 10x background points → `data/sightings_with_env.rds`,
   `data/background_with_env.rds`.
4. `analysis/02_fit_suitability.R` — seasonal presence-background GAMs
   (breeding Jun-Sep, feeding Nov-Feb; background capped to presence SST ±1 °C
   and lat ±2°) → `data/thermal_envelopes.csv`, `data/suitability_curves.csv`,
   `data/suitability_models.rds`, `analysis/thermal_envelopes.png`.

`analysis/01b_diagnostics.R` is an optional assumption-check pass (SST
contrast, Chl-a missingness, effort bias) — run any time after step 3.

## Fitted thermal envelopes (2026-10-10 run)

| season   | n presences | T_opt (°C) | sigma_T (°C) |
|----------|-------------|------------|--------------|
| breeding | 12,486      | 21.3       | 4.6          |
| feeding  | 8,566       | -0.9       | 6.1          |

Feeding T_opt sits at the cold data edge: interpret the feeding envelope as a
broad cold-water envelope (suitability ~1 below ~12 °C, 0 above ~17 °C).

## Data provenance

- `data-raw/sst_southern_2003_2022.csv`: IMOS SRS GHRSST L3S-1m/day (monthly
  AVHRR SST composites), AODN THREDDS NetcdfSubset, lat -67 to -39,
  lon 143 to 170, aggregated to the 2° grid. 240 months, 51,263 cell-months,
  no gaps. For 2023+ use IMOS/SRS/SST/ghrsst/L4/RAMSSA (daily, gap-free).
- `data-raw/humpback_southern_ocean_obis.csv`: OBIS API v3,
  *Megaptera novaeangliae*, lat -65 to -45, lon 80-180 (7,824 records;
  7,551 in Nov-Feb feeding season, 6,909 post-2003).
- Salinity (SODA, to 2015) is unused; Chl-a retained only for feeding-season
  fits (breeding coverage ~47% missing).
