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
