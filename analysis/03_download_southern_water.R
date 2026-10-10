# 03_download_southern_water.R
# Extend the water grid south: monthly SST for lat -60 to -38, lon 143 to 170,
# 2003-2022, from IMOS/AODN. Output format matches the Water sheet so scripts
# 01/02 run unchanged on the combined dataset.
#
# DATA PROVENANCE (verified 2026-10-10 against the live AODN THREDDS catalogue):
#   Product : IMOS/SRS/SST/ghrsst/L3S-1m/day  (monthly multi-sensor AVHRR SST,
#             skin temperature, one file per month per year, 2003-2022)
#   Server  : https://thredds.aodn.org.au/thredds/
#   Access  : NetcdfSubset (ncss/grid) server-side subset:
#             var=sea_surface_temperature, north=-38, south=-60, west=143,
#             east=170, horizStride=10 (~2 km), accept=netcdf4
#   For 2023-2026 use IMOS/SRS/SST/ghrsst/L4/RAMSSA (daily, gap-free,
#   RAMSSA_09km-AUS) -- aggregate daily files to monthly means.
#
# The actual download is done by the companion Python script
# (analysis/fetch_sst_southern.py) because Windows R/ncdf4 OPeNDAP is fragile.
# It writes data-raw/sst_southern_2003_2022.csv, which THIS script merges with
# the existing Water sheet.
#
# Run order: python analysis/fetch_sst_southern.py   (once, ~15 min)
#            then source this file in R.
#
# Output: data/water_extended.rds

library(readxl)
library(dplyr)

RAW   <- "data-raw/Humpback_East_Australia_data.xlsx"
SST_S <- "data-raw/sst_southern_2003_2022.csv"

water <- read_excel(RAW, sheet = "Water") |>
  mutate(year = as.integer(year), month = as.integer(month))

south <- read.csv(SST_S) |>
  transmute(year = as.integer(year), month = as.integer(month),
            cell_centre_lat, cell_centre_lon,
            sea_surface_temp_C = sst,
            chlorophyll_mg_per_m3 = NA_real_,
            salinity_psu_SODA_to2015 = NA_real_)

# keep original rows where they exist; add southern rows only where new
south_new <- anti_join(south, water,
                       by = c("year", "month", "cell_centre_lat", "cell_centre_lon"))

extended <- bind_rows(water, south_new) |>
  arrange(year, month, cell_centre_lat, cell_centre_lon)

cat("Original water rows:", nrow(water), "\n")
cat("Southern rows added:", nrow(south_new), "\n")
cat("Extended total:", nrow(extended), "\n")
cat("Lat coverage now:", min(extended$cell_centre_lat), "to",
    max(extended$cell_centre_lat), "\n")

saveRDS(extended, "data/water_extended.rds")
cat("Saved data/water_extended.rds\n")
