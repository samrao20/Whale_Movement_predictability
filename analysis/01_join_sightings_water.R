# 01_join_sightings_water.R
# Join humpback sightings (2002+) to monthly water-grid covariates (SST, Chl-a)
# and generate background (pseudo-absence) points for the suitability model.
#
# Inputs : data-raw/Humpback_East_Australia_data.xlsx
# Outputs: data/sightings_with_env.rds   (presences, 1 row per sighting)
#          data/background_with_env.rds  (background points, 10x presences)

library(readxl)
library(dplyr)

RAW <- "data-raw/Humpback_East_Australia_data.xlsx"
MAX_DIST_DEG <- 1.5   # max distance to nearest water cell (~165 km)

# Use combined sightings (original + Southern Ocean) when script 04 has run.
if (file.exists("data-raw/humpback_all.rds")) {
  hump <- readRDS("data-raw/humpback_all.rds")
  cat("Using combined sightings (incl. Southern Ocean)\n")
} else {
  hump <- read_excel(RAW, sheet = "Humpback")
  cat("Using original Humpback sheet only\n")
}

# Use the extended southern grid (from script 03) if it has been built;
# otherwise fall back to the original Water sheet.
if (file.exists("data/water_extended.rds")) {
  water <- readRDS("data/water_extended.rds")
  cat("Using extended water grid (includes Southern Ocean extension)\n")
} else {
  water <- read_excel(RAW, sheet = "Water")
  cat("Using original Water sheet (no southern extension found)\n")
}

# --- Clean sightings ---------------------------------------------------------
sightings <- hump |>
  filter(year >= 2002, !is.na(latitude), !is.na(longitude),
         !is.na(month), month >= 1, month <= 12) |>
  mutate(year = as.integer(year), month = as.integer(month))

# --- Join water covariates ---------------------------------------------------
water <- water |>
  mutate(year = as.integer(year), month = as.integer(month)) |>
  select(year, month, cell_centre_lat, cell_centre_lon,
         sst = sea_surface_temp_C, chl = chlorophyll_mg_per_m3)

# Nearest water-cell join per year-month: recovers coastal sightings whose
# 2-degree snap cell is on land. Exact-cell matches get distance 0.
nearest_join <- function(sight_df, water_df) {
  out <- list()
  for (ym in unique(paste(sight_df$year, sight_df$month))) {
    s <- sight_df[paste(sight_df$year, sight_df$month) == ym, ]
    w <- water_df[paste(water_df$year, water_df$month) == ym, ]
    if (nrow(w) == 0) next
    d <- outer(s$latitude, w$cell_centre_lat, `-`)^2 +
         outer(s$longitude, w$cell_centre_lon, `-`)^2
    idx <- apply(d, 1, which.min)
    dist <- sqrt(d[cbind(seq_len(nrow(s)), idx)])
    keep <- dist <= MAX_DIST_DEG
    s2 <- s[keep, ]
    s2$cell_lat <- w$cell_centre_lat[idx[keep]]
    s2$cell_lon <- w$cell_centre_lon[idx[keep]]
    s2$sst <- w$sst[idx[keep]]
    s2$chl <- w$chl[idx[keep]]
    s2$dist_to_cell <- dist[keep]
    out[[length(out) + 1]] <- s2
  }
  bind_rows(out)
}

cat("Joining sightings to nearest water cell (tolerance", MAX_DIST_DEG, "deg)...\n")
sightings_env <- nearest_join(sightings, water)

cat("Sightings kept:", nrow(sightings_env), "of", nrow(sightings),
    sprintf("(%.0f%%)", 100 * nrow(sightings_env) / nrow(sightings)), "\n")
cat("Median dist to cell (deg):", round(median(sightings_env$dist_to_cell), 2), "\n")

# --- Background points -------------------------------------------------------
# Sample available environment: same grid cells, same years/months, weighted by
# how often each cell-month appears in the water data (effort proxy v1).
set.seed(42)
n_bg <- nrow(sightings_env) * 10

cell_months <- water |>
  filter(!is.na(sst)) |>
  distinct(year, month, cell_centre_lat, cell_centre_lon)

background <- cell_months |>
  slice_sample(n = min(n_bg, nrow(cell_months))) |>
  left_join(water, by = c("year", "month",
                          "cell_centre_lat", "cell_centre_lon")) |>
  rename(cell_lat = cell_centre_lat, cell_lon = cell_centre_lon)

# --- Save --------------------------------------------------------------------
dir.create("data", showWarnings = FALSE)
saveRDS(sightings_env, "data/sightings_with_env.rds")
saveRDS(background,    "data/background_with_env.rds")

cat("Saved data/sightings_with_env.rds and data/background_with_env.rds\n")

# Quick sanity check: seasonal signal in sighting latitude (migration pulse)
print(
  sightings_env |>
    group_by(month) |>
    summarise(n = n(), median_lat = median(latitude), .groups = "drop")
)
