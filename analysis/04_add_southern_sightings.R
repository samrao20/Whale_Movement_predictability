# 04_add_southern_sightings.R
# Append Southern Ocean humpback records (OBIS, lat -45..-65) to the sightings
# used by the pipeline, formatted to match the Humpback sheet.
#
# Input : data-raw/humpback_southern_ocean_obis.csv  (OBIS API pull, 2026-10-10)
#         data-raw/Humpback_East_Australia_data.xlsx
# Output: data-raw/humpback_all.rds  (original + southern records combined)
#
# Run BEFORE 01_join_sightings_water.R; that script reads data-raw/humpback_all.rds
# when present (patch below adds the switch).

library(readxl)
library(dplyr)

orig <- read_excel("data-raw/Humpback_East_Australia_data.xlsx", sheet = "Humpback")

south <- read.csv("data-raw/humpback_southern_ocean_obis.csv") |>
  transmute(
    source_database  = "OBIS-SouthernOcean",
    record_id        = occurrenceID,
    date             = as.Date(substr(eventDate, 1, 10)),
    year             = suppressWarnings(as.numeric(year)),
    month            = suppressWarnings(as.numeric(month)),
    latitude         = decimalLatitude,
    longitude        = decimalLongitude,
    individual_count = suppressWarnings(as.numeric(individualCount)),
    depth_m          = NA_real_,
    basis_of_record  = basisOfRecord,
    dataset          = datasetName,
    institution      = institutionCode
  )

cat("Southern records:", nrow(south),
    "| with year+month:", sum(!is.na(south$year) & !is.na(south$month)), "\n")
cat("Feeding-season (Nov-Feb):", sum(south$month %in% c(11,12,1,2), na.rm = TRUE), "\n")

combined <- bind_rows(orig, south)

# de-duplicate: same lat/lon/date appearing in both extracts
before <- nrow(combined)
combined <- combined |>
  distinct(latitude, longitude, date, .keep_all = TRUE)
cat("De-duplicated:", before, "->", nrow(combined), "\n")

saveRDS(combined, "data-raw/humpback_all.rds")
cat("Saved data-raw/humpback_all.rds\n")
