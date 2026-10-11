# Whale Movement Predictability

Building a model / formula to quantify the predictability of whale movement patterns on average, using
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

## Satellite tag data

`data-raw/tags/` holds the public tag dataset (Andrews-Goff et al. 2023, Movebank,
doi:10.5441/001/1.294, CC0): 50 east Australian humpbacks (breeding stock E1) tagged
2008 to 2010, 30,388 locations (13,311 smoothed state-space-model positions, 17,077 raw
Argos fixes). `data-raw/05_make_tag_excel.py` builds `data/humpback_satellite_tags.xlsx`
(sheets: Read_me, Whales, Tracks_modelled, Tracks_raw_argos). 48 whales have modelled
tracks. Please cite the dataset and paper when using it.

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

## Tag migration vs IMOS environment (`analysis/05_tag_migration_by_year.py`)

Matches every smoothed satellite-tag position to IMOS data and compares the tagging campaigns.

| step | script | output |
|---|---|---|
| IMOS monthly SST + SST front strength at each tag position | `data-raw/06_tag_imos_sst.py` | `data/tag_points_env.csv` |
| IMOS GSLA sea level / geostrophic currents (every 7th day, nearest day within 3 days, lat -60 to -10 only) | `data-raw/07_tag_imos_currents.py` | `data/tag_points_currents.csv` |
| Per-whale southbound metrics and figures | `analysis/05_tag_migration_by_year.py` | `data/tag_migration_metrics.csv`, `analysis/tag_*.png` |

The tags are four separate campaigns, not four comparable years: Oct-Nov 2008 (15 whales off Eden,
southbound), Jun 2009 (13, Queensland, northbound), Feb-Mar 2010 (14, tagged in Antarctic waters),
Oct 2010 (6, Queensland, southbound, tracked to Jul 2011). Only 2008 vs 2010 southbound can be
compared, and they started 11 degrees of latitude apart.

Findings so far (suggestive only): both groups crossed 40S in early November (5 Nov 2008, 10 Nov 2010);
the 2010 whales reached 60S about two weeks earlier (23 Nov vs 8 Dec) and moved faster between 40S and
60S (about 112 vs 92 km/day), but that rests on 3 whales in 2010 vs 8 in 2008. SST along the route was
almost identical (about 9.3 vs 9.5 C on the 40-60S leg). Front-strength and current curves are noisy
for 2010 because few whales contribute to each latitude bin.

Notes: SST below -1.9 C is treated as ice/retrieval artefact; about 12% of points have no SST (ice,
cloud). Chlorophyll is not included (IMOS only provides daily MODIS files, mostly cloud-covered at
these latitudes). Tracks cross the dateline (longitude is shown as 0-360).

## What drives tag movement, by group (`analysis/06_tag_movement_drivers.py`)

Adds chlorophyll (`data-raw/08_tag_chlorophyll.py`: ESA CCI monthly ocean colour via NOAA ERDDAP,
`data/tag_points_chl.csv`, 87% of points matched; gaps are polar night / ice in Apr-May 2010) and
asks what is associated with how fast whales travel, within groups = tagging campaign x latitude zone
(N of 45S, 45-60S, S of 60S).

Method: one position per whale per day; daily travel speed; drivers z-scored within group (effect per
+1 SD): SST anomaly and chlorophyll anomaly (both relative to the latitude/month norm, because raw SST is
99% explained by latitude and season), SST front strength, southward current (IMOS GSLA, N of 60S), plus latitude as a
control. Speed uses a linear mixed model with whale as random intercept; "tortuous path" (path
straightness < 0.5 over 4 days) uses a logistic GEE with whale clusters.
Outputs: `data/tag_driver_models.csv`, `data/tag_daily_steps.csv`, `analysis/tag_driver_heatmap.png`,
`analysis/tag_driver_binned.png`.

Main results (speed model): higher chlorophyll goes with slower travel for the 2008 whales in all three
zones (-0.24, -0.10, -0.18 log-speed per SD); a stronger southward current goes with faster travel in
three groups (+0.08 to +0.14); SST anomaly is positive for the 2008 whales N of 45S only; SST fronts show no
consistent effect; the Antarctic feeding groups show nothing significant. 11 of 46 non-latitude effects have
intervals excluding 0 (about 24%), more than chance alone but unadjusted for multiple tests.

Cautions: associations, not causes. Daily positions are autocorrelated and only 3-15 whales feed each group
(65-284 whale-days), so intervals are likely too narrow. The tortuous-path model is unstable (huge latitude
coefficients, likely near-separation); treat it as exploratory. Monthly SST/chlorophyll cannot see
day-scale prey patches, and the binned plot pools all latitudes, so it is confounded by latitude.

## References

See [docs/references.md](docs/references.md) for papers, datasets and similar studies. Add to it as you find more.
