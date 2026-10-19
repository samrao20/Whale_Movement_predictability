# ============================================================================
# 07_pde_currents.R
# Script 05 + IMOS GSLA currents in the advection term.
#
#   v(x,t) = v_drift(month) + v_current(x,month)
#
# Currents: data/currents_cells.csv (IMOS GSLA u/v per cell-month, 2002-2015,
# climatologised here by lat x month). v_ms is + northward; grid velocity is
# + poleward, so v_current = -v_ms * 3.6 km/h.
#
# Domain note: currents data only covers -45..-11, so the coastal strip here
# is -29..-45 (currents-blended); south of -45 the drift-only solution stands.
#
# Run: setwd() to repo root, then source("analysis/07_pde_currents.R")
# ============================================================================

library(ggplot2)

# Fail loudly if the working directory is wrong (the #1 cause of
# "no rows to aggregate": read.csv fails, a stale empty object gets reused)
if (!file.exists("data/water_extended.rds"))
  stop("Working directory is not the repo root.\n",
       "  In RStudio: Files pane -> navigate into Whale_Movement_predictability ->\n",
       "  More -> Set As Working Directory, then re-source this script.\n",
       "  Current getwd(): ", getwd())

# ---- 1. Forcing climatology (same as 05) ------------------------------------
water <- readRDS("data/water_extended.rds")
strip <- water[water$cell_centre_lon >= 145 & water$cell_centre_lon <= 156, ]

grid <- seq(-29, -45, by = -2)          # currents coverage cuts domain here
N <- length(grid)
dx_km <- 2 * 111

clim <- aggregate(cbind(sea_surface_temp_C, chlorophyll_mg_per_m3) ~
                    cell_centre_lat + month, data = strip, FUN = mean, na.rm = TRUE)
chl_fill <- median(clim$chlorophyll_mg_per_m3, na.rm = TRUE)

sst <- matrix(NA_real_, N, 12); chl <- matrix(NA_real_, N, 12)
for (i in seq_len(N)) for (m in 1:12) {
  row <- clim[clim$cell_centre_lat == grid[i] & clim$month == m, ]
  if (nrow(row) == 1) {
    sst[i, m] <- row$sea_surface_temp_C
    chl[i, m] <- ifelse(is.na(row$chlorophyll_mg_per_m3), chl_fill,
                        row$chlorophyll_mg_per_m3)
  } else chl[i, m] <- chl_fill
}
for (m in 1:12) {
  col <- sst[, m]
  if (any(is.na(col)))
    sst[, m] <- approx(seq_len(N), col, xout = seq_len(N), rule = 2, ties = mean)$y
}
chl[is.na(chl)] <- chl_fill
if (any(is.na(sst))) stop("sst still has NA after imputation")

# ---- 2. Current climatology ---------------------------------------------------
cur <- read.csv("data/currents_cells.csv")
cur <- cur[cur$cell_centre_lon >= 145 & cur$cell_centre_lon <= 156, ]
cur_clim <- aggregate(v_ms ~ cell_centre_lat + month, data = cur,
                      FUN = mean, na.rm = TRUE)
cur_mat <- matrix(0, N, 12)
for (i in seq_len(N)) for (m in 1:12) {
  row <- cur_clim[cur_clim$cell_centre_lat == grid[i] & cur_clim$month == m, ]
  cur_mat[i, m] <- if (nrow(row) == 1) row$v_ms else 0
}
cur_kmh <- -cur_mat * 3.6      # + poleward, km/h

# ---- 3. Parameters -------------------------------------------------------------
track_params <- if (file.exists("data/track_params.rds")) readRDS("data/track_params.rds") else NULL
v_max <- if (is.null(track_params)) 4.9 else track_params$v_max
D0    <- if (is.null(track_params)) 10 else track_params$D0
h     <- 0.4

hours_per_month <- 365 / 12 * 24
dt <- 0.25
total_months <- 36
n_steps <- round(total_months * hours_per_month / dt)

# ---- 4. PDE step (upwind advection + centred diffusion, zero-flux ends) ------
pde_step <- function(C, v, D) {
  vf <- 0.5 * (v[-N] + v[-1])
  up <- ifelse(vf > 0, C[-N], C[-1])
  adv_flux <- vf * up
  Df <- 0.5 * (D[-N] + D[-1])
  dif_flux <- -Df * (C[-1] - C[-N]) / dx_km
  net <- adv_flux + dif_flux
  # zero-flux boundaries: end faces carry no flux, so boundary cells keep the
  # one interior face they own. Freezing dC at boundary cells (as script 05
  # does) silently leaks mass through those faces (~2/3 over 3 years here).
  dC <- numeric(N)
  dC[1] <- -net[1] / dx_km
  dC[N] <-  net[N - 1] / dx_km
  dC[2:(N - 1)] <- -(net[2:(N - 1)] - net[1:(N - 2)]) / dx_km
  C + dt * dC
}

# ---- 5. Time loop --------------------------------------------------------------
drift_obs <- c("1" = -1.39, "2" = -1.06, "3" = -0.88, "4" = -0.87,
               "5" =  1.06, "6" =  4.71, "7" =  1.90, "8" =  1.99,
               "9" =  0.00, "10" = -2.55, "11" = -2.55, "12" = -2.87)

C <- numeric(N); C[(N-2):N] <- 1
t_months <- 3
snaps <- list()

for (s in seq_len(n_steps)) {
  m  <- t_months %% 12
  m0 <- floor(m) %% 12 + 1
  m1 <- m0 %% 12 + 1
  w  <- m - floor(m)

  drift <- (1 - w) * drift_obs[as.character(m0)] + w * drift_obs[as.character(m1)]
  cur_v <- (1 - w) * cur_kmh[, m0] + w * cur_kmh[, m1]
  v <- rep(-drift, N) + cur_v

  P <- (1 - w) * chl[, m0] + w * chl[, m1]
  D <- D0 / (1 + P / h)

  C <- pde_step(C, v, D)
  t_months <- t_months + dt / hours_per_month

  if (s %% round(hours_per_month / dt) == 0)
    snaps[[length(snaps) + 1]] <- data.frame(month = t_months, lat = grid, density = C)
}

out <- do.call(rbind, snaps)
out$year <- floor(out$month / 12) + 1
out$month_of_year <- round(out$month) %% 12   # fp drift in t_months
out$month_of_year[out$month_of_year == 0] <- 12

# ---- 6. Checks + plots -----------------------------------------------------------
mass <- tapply(out$density, out$month, sum)
cat(sprintf("Mass conservation: %.3f .. %.3f (initial 3.000)\n", min(mass), max(mass)))
cat(sprintf("Min density: %.5f\n", min(out$density)))

yr3 <- out[out$year == 3, ]
centroid <- tapply(seq_len(nrow(yr3)), yr3$month_of_year, function(idx) {
  d <- yr3$density[idx]; la <- yr3$lat[idx]
  if (sum(d) <= 0) return(NA_real_)
  sum(la * d) / sum(d)
})
cat("\nYear-3 monthly centroid (latitude), drift + currents:\n")
print(round(centroid, 1))

p1 <- ggplot(yr3, aes(lat, density, group = month_of_year, colour = factor(month_of_year))) +
  geom_line() + scale_colour_viridis_d(name = "month") +
  labs(title = "Year 3: density along corridor (drift + IMOS currents)",
       x = "latitude", y = "density (relative)") +
  theme_minimal()
print(p1)
ggsave("analysis/pde_1d_currents_year3.png", p1, width = 8, height = 5)

p2 <- ggplot(out, aes(month, lat, fill = density)) +
  geom_raster() + scale_fill_viridis_c() +
  labs(title = "Hovmoller: drift + currents", x = "months since start", y = "latitude") +
  theme_minimal()
print(p2)
ggsave("analysis/pde_1d_currents_hovmoller.png", p2, width = 8, height = 5)

saveRDS(out, "data/pde_1d_currents_output.rds")
cat("\nSaved data/pde_1d_currents_output.rds + two PNGs in analysis/\n")
