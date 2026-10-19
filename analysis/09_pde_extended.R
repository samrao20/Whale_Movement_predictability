# ============================================================================
# 09_pde_extended.R
# Extended-domain PDE: -29..-65 (full traverse), drift + currents + ARS.
#
# v(x,t) = tag drift(month) + IMOS current(x,month)   [currents only -29..-45,
#          zero south of -45: GSLA product doesn't reach]
# Forcing: climatological SST/Chl from data/water_extended.rds (lon 145-156).
# Chl is entirely missing south of -55 -> median fallback (ARS diffusion is
# near-D0 there anyway; whales are resident on feeding grounds).
#
# Validation in-script: model centroid by month (year 3) vs observed sighting
# centroid by month from data-raw/humpback_all.rds (all records, year >= 2002,
# corridor lon 145-156 OR southern-ocean pull lon 80-180).
#
# Run: setwd() to repo root, then source("analysis/09_pde_extended.R")
# ============================================================================

library(ggplot2)

if (!file.exists("data/water_extended.rds"))
  stop("Working directory is not the repo root. Set it via Files pane -> ",
       "More -> Set As Working Directory. Current getwd(): ", getwd())

water <- readRDS("data/water_extended.rds")
strip <- water[water$cell_centre_lon >= 145 & water$cell_centre_lon <= 156, ]

grid <- seq(-29, -65, by = -2)
N <- length(grid)
dx_km <- 2 * 111

chl_fill <- median(strip$chlorophyll_mg_per_m3, na.rm = TRUE)

sst <- chl <- matrix(NA_real_, N, 12)
a_sst <- aggregate(sea_surface_temp_C ~ cell_centre_lat + month,
                   data = strip, FUN = mean, na.rm = TRUE)
strip_chl <- strip[!is.na(strip$chlorophyll_mg_per_m3), ]
a_chl <- aggregate(chlorophyll_mg_per_m3 ~ cell_centre_lat + month,
                   data = strip_chl, FUN = mean, na.rm = TRUE)
for (i in seq_len(N)) for (m in 1:12) {
  j <- which(a_sst$cell_centre_lat == grid[i] & a_sst$month == m)
  if (length(j) == 1) sst[i, m] <- a_sst$sea_surface_temp_C[j]
  j <- which(a_chl$cell_centre_lat == grid[i] & a_chl$month == m)
  if (length(j) == 1) chl[i, m] <- a_chl$chlorophyll_mg_per_m3[j]
}
chl[is.na(chl)] <- chl_fill
for (m in 1:12) {                       # interpolate any missing SST along lat
  col <- sst[, m]
  if (any(is.na(col)))
    sst[, m] <- approx(seq_len(N), col, xout = seq_len(N), rule = 2)$y
}
if (any(is.na(sst))) stop("sst still has NA after imputation")

# ---- currents climatology (north of -45 only) --------------------------------
cur <- read.csv("data/currents_cells.csv")
cur <- cur[cur$cell_centre_lon >= 145 & cur$cell_centre_lon <= 156, ]
cur_clim <- aggregate(v_ms ~ cell_centre_lat + month, data = cur,
                      FUN = mean, na.rm = TRUE)
cur_kmh <- matrix(0, N, 12)
for (i in seq_len(N)) for (m in 1:12) {
  j <- which(cur_clim$cell_centre_lat == grid[i] & cur_clim$month == m)
  if (length(j) == 1) cur_kmh[i, m] <- -cur_clim$v_ms[j] * 3.6   # + poleward
}

# ---- parameters / step -------------------------------------------------------
track_params <- if (file.exists("data/track_params.rds")) readRDS("data/track_params.rds") else NULL
D0 <- if (is.null(track_params)) 10 else track_params$D0
h  <- 0.4

pde_step <- function(C, v, D, dt) {
  vf <- 0.5 * (v[-N] + v[-1])
  up <- ifelse(vf > 0, C[-N], C[-1])
  net <- vf * up - 0.5 * (D[-N] + D[-1]) * (C[-1] - C[-N]) / dx_km
  dC <- numeric(N)
  dC[1] <- -net[1] / dx_km
  dC[N] <-  net[N - 1] / dx_km
  dC[2:(N - 1)] <- -(net[2:(N - 1)] - net[1:(N - 2)]) / dx_km
  C + dt * dC
}

drift_obs <- c("1" = -1.39, "2" = -1.06, "3" = -0.88, "4" = -0.87,
               "5" =  1.06, "6" =  4.71, "7" =  1.90, "8" =  1.99,
               "9" =  0.00, "10" = -2.55, "11" = -2.55, "12" = -2.87)

hours_per_month <- 365 / 12 * 24
dt <- 0.5
n_steps <- round(36 * hours_per_month / dt)

C <- numeric(N); C[(N - 2):N] <- 1      # pulse on -61..-65 (March, far south)
t_months <- 3
snaps <- list()
for (s in seq_len(n_steps)) {
  m  <- t_months %% 12
  m0 <- floor(m) %% 12 + 1
  m1 <- m0 %% 12 + 1
  w  <- m - floor(m)
  drift <- (1 - w) * drift_obs[as.character(m0)] + w * drift_obs[as.character(m1)]
  v <- rep(-drift, N) + (1 - w) * cur_kmh[, m0] + w * cur_kmh[, m1]
  D <- D0 / (1 + ((1 - w) * chl[, m0] + w * chl[, m1]) / h)
  C <- pde_step(C, v, D, dt)
  t_months <- t_months + dt / hours_per_month
  if (s %% round(hours_per_month / dt) == 0)
    snaps[[length(snaps) + 1]] <- data.frame(month = t_months, lat = grid, density = C)
}
out <- do.call(rbind, snaps)
out$year <- floor(out$month / 12) + 1
# t_months accumulates floating-point drift (1.00000000001571 not 1), so the
# calendar month must be rounded, not truncated, before grouping
out$month_of_year <- round(out$month) %% 12
out$month_of_year[out$month_of_year == 0] <- 12

mass <- tapply(out$density, out$month, sum)
cat(sprintf("Mass conservation: %.3f .. %.3f (initial 3.000)\n", min(mass), max(mass)))
cat(sprintf("Min density: %.5f\n", min(out$density)))

# ---- validation vs sightings ----------------------------------------------------
all_s <- readRDS("data-raw/humpback_all.rds")
lat_col <- intersect(c("decimalLatitude", "latitude", "lat"), names(all_s))[1]
lon_col <- intersect(c("decimalLongitude", "longitude", "lon"), names(all_s))[1]
mo_col  <- intersect(c("month"), names(all_s))[1]
yr_col  <- intersect(c("year"), names(all_s))[1]
keep <- all_s[[yr_col]] >= 2002 &
  ((all_s[[lon_col]] >= 145 & all_s[[lon_col]] <= 156 & all_s[[lat_col]] >= -45) |
     (all_s[[lon_col]] >= 80  & all_s[[lon_col]] <= 180 & all_s[[lat_col]] < -45))
sight <- all_s[keep, ]
cat(sprintf("Sightings used for validation: %d records (year>=2002)\n", nrow(sight)))

obs_cent <- tapply(sight[[lat_col]], sight[[mo_col]], mean)
yr3 <- out[out$year == 3, ]
mod_cent <- tapply(seq_len(nrow(yr3)), yr3$month_of_year, function(idx) {
  d <- yr3$density[idx]
  sum(yr3$lat[idx] * d) / sum(d)
})
cmp <- data.frame(month = 1:12,
                  model = as.numeric(mod_cent[as.character(1:12)]),
                  observed = as.numeric(obs_cent[as.character(1:12)]))
cmp$err <- cmp$model - cmp$observed
cat("\nMonth  model  observed  err (deg lat)\n")
print(round(cmp, 1))
cat(sprintf("\nMean err %.1f, mean |err| %.1f, corr %.2f\n",
            mean(cmp$err), mean(abs(cmp$err)),
            cor(cmp$model, cmp$observed, use = "complete.obs")))

p1 <- ggplot(yr3, aes(lat, density, group = month_of_year, colour = factor(month_of_year))) +
  geom_line() + scale_colour_viridis_d(name = "month") +
  labs(title = "Year 3: density along extended corridor (-29..-65)",
       x = "latitude", y = "density (relative)") +
  theme_minimal()
print(p1)
ggsave("analysis/pde_extended_year3.png", p1, width = 8, height = 5)

p2 <- ggplot(out, aes(month, lat, fill = density)) +
  geom_raster() + scale_fill_viridis_c() +
  labs(title = "Hovmoller: extended domain, 3 years",
       x = "months since start", y = "latitude") +
  theme_minimal()
print(p2)
ggsave("analysis/pde_extended_hovmoller.png", p2, width = 8, height = 5)

p3 <- ggplot(cmp, aes(month)) +
  geom_line(aes(y = model, colour = "model"), linewidth = 1) +
  geom_line(aes(y = observed, colour = "observed"), linewidth = 1) +
  scale_x_continuous(breaks = 1:12) +
  labs(title = "Monthly centroid: model vs sightings", y = "latitude", colour = NULL) +
  theme_minimal()
print(p3)
ggsave("analysis/pde_extended_centroids.png", p3, width = 8, height = 4)

saveRDS(list(output = out, centroid_comparison = cmp), "data/pde_extended_output.rds")
cat("\nSaved data/pde_extended_output.rds + three PNGs in analysis/\n")
