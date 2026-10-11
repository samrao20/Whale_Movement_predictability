# 01b_diagnostics.R
# Check that the joined data meets the assumptions of the presence-background
# suitability model BEFORE fitting. Prints a PASS/WARN/FAIL line per check.
#
# Run after: analysis/01_join_sightings_water.R

library(dplyr)
library(ggplot2)

sightings  <- readRDS("data/sightings_with_env.rds")
background <- readRDS("data/background_with_env.rds")

verdict <- function(ok, msg) cat(if (ok) "PASS " else "WARN ", msg, "\n", sep = "")

season_of <- function(m) case_when(
  m %in% 6:9             ~ "breeding",
  m %in% c(11, 12, 1, 2) ~ "feeding",
  TRUE                   ~ "transit")

pres <- sightings |>
  filter(!is.na(sst)) |>
  mutate(season = season_of(month))

cat("================ 1. SAMPLE SIZES ================\n")
n_season <- pres |> count(season)
print(n_season)
for (s in c("breeding", "feeding")) {
  n <- n_season$n[n_season$season == s]
  verdict(n >= 100, sprintf("%s presences = %d (need >= 100 for a stable GAM)", s, n))
}

cat("\n================ 2. SST CONTRAST (presences vs background) ================\n")
# The model can only estimate a thermal envelope if presences and background
# span a similar SST range AND differ in distribution within it.
bg <- background |>
  filter(!is.na(sst)) |>
  mutate(season = season_of(month))

sst_summary <- bind_rows(
  pres |> mutate(type = "presence"),
  bg   |> mutate(type = "background", sst = sst)
) |>
  filter(season %in% c("breeding", "feeding")) |>
  group_by(season, type) |>
  summarise(min = min(sst), q25 = quantile(sst, .25), median = median(sst),
            q75 = quantile(sst, .75), max = max(sst), n = n(), .groups = "drop")
print(as.data.frame(sst_summary))

ks <- pres |> filter(season == "breeding") |> pull(sst)
verdict(length(unique(round(ks, 1))) >= 8,
        sprintf("breeding presences span %d distinct SST values (need >= 8)",
                length(unique(round(ks, 1)))))

cat("\n-- presences at range EDGE (risk of boundary extrapolation) --\n")
edge <- pres |>
  group_by(season) |>
  summarise(frac_warmest_1C = mean(sst > max(sst) - 1),
            frac_coldest_1C = mean(sst < min(sst) + 1), .groups = "drop")
print(as.data.frame(edge))

cat("\n================ 3. CHL-A MISSINGNESS BY SEASON ================\n")
chl_miss <- pres |>
  group_by(season) |>
  summarise(pct_chl_missing = round(100 * mean(is.na(chl)), 1), .groups = "drop")
print(as.data.frame(chl_miss))

cat("\n================ 4. TEMPORAL DISTRIBUTION (effort bias) ================\n")
by_decade <- pres |>
  mutate(decade = floor(year / 10) * 10) |>
  count(decade)
print(as.data.frame(by_decade))
recent_frac <- mean(pres$year >= 2010)
verdict(recent_frac > 0.5,
        sprintf("%.0f%% of presences are post-2010 (recent-dominated is fine for a climate model, but check old records aren't distorting)", 100 * recent_frac))

cat("\n================ 5. SPATIAL CLUMPING / DUPLICATES ================\n")
dups <- pres |>
  count(year, month, cell_lat = round(latitude), cell_lon = round(longitude)) |>
  filter(n > 1)
cat("grid cells with >1 presence in same month:", nrow(dups),
    "(", round(100 * nrow(dups) / nrow(pres), 1), "% of records)\n")

cat("\n================ 6. PRESENCE LOCATIONS vs WATER GRID COVERAGE ================\n")
# Presences on land or outside water cells got NA sst and were dropped above
dropped <- sightings |> filter(is.na(sst)) |> nrow()
total   <- nrow(sightings)
verdict(dropped / total < 0.1,
        sprintf("%d of %d sightings (%.1f%%) had no matching water cell (land/edge) and were dropped",
                dropped, total, 100 * dropped / total))

cat("\n================ 7. SST RANGE OF PRESENCES ================\n")
rng <- pres |> group_by(season) |>
  summarise(sst_min = round(min(sst), 1), sst_max = round(max(sst), 1),
            span = round(max(sst) - min(sst), 1), .groups = "drop")
print(as.data.frame(rng))
verdict(all(rng$span[rng$season %in% c("breeding", "feeding")] >= 5),
        "each season's presences span >= 5 C (enough contrast to fit an envelope)")

# --- Visual check ------------------------------------------------------------
p <- ggplot(bind_rows(pres |> mutate(type = "presence"),
                      bg   |> mutate(type = "background")) |>
              filter(season %in% c("breeding", "feeding")),
            aes(sst, fill = type)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~season, scales = "free_y") +
  labs(title = "SST distributions: presences vs available background",
       x = "SST (C)", y = "density") +
  theme_minimal()
ggsave("analysis/diagnostics_sst_density.png", p, width = 8, height = 4)
cat("\nSaved analysis/diagnostics_sst_density.png -- open it to inspect overlap\n")
