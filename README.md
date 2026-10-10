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
