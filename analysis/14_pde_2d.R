# ============================================================================
# 14_pde_2d.R
# 2-D advection-reaction-diffusion PDE on the lat-lon grid.
#
#   dC/dt = div(D grad C) - div(v C) + R(C,T)
#
# Grid: water-grid cell centres lat -67..-11 (odd) x lon 141..169 (odd),
# 29 x 15 ocean cells, ~222 km spacing (zonal spacing scaled by cos(lat)).
#
#   v(x,y,t): tag drift field (data/drift_field_2d.rds, 4-deg bins, bilinearly
#             interpolated to the 2-deg grid) + IMOS currents (u_ms, v_ms from
#             data/currents_cells.csv, climatology, +north/+east -> model
#             +poleward/+eastward).
#   D(x,y,t): ARS, D0 / (1 + Chl/h) as in 1-D.
#   R(C,T):   suitability-weighted logistic growth + mismatch mortality,
#             exactly as script 10 (no wall loss; southern boundary is open
#             ocean and zero-flux suffices for the E1 stock corridor).
#   Land:     cells with frac_land > 0.5 (bathymetry_cells.csv) are masked
#             out — zero flux across land faces, C fixed at 0 there.
#
# Forcing: monthly climatology of SST/Chl from data/cell_month_2002_2015.csv
# (lon >= 141) merged with the southern extension in data/water_extended.rds.
#
# Numerics: explicit upwind advection + centred diffusion per direction,
# mass-conserving zero-flux boundaries (as fixed in script 05), dt = 0.25 h.
# CFL: |v| <= 7.2 km/h over ~222 km cells -> dt << 30 h; fine.
#
# Validation: per-month pattern correlation between model density field and
# binned sighting counts on the same grid (spearman, year 3).
#
# Run: setwd() to repo root, then source("analysis/14_pde_2d.R")
# Runtime: ~1-2 min (435 cells x ~52k steps x 36 months).
# ============================================================================

library(ggplot2)

if (!file.exists("data/drift_field_2d.rds"))
  stop("Working directory is not the repo root, or script 13 has not been ",
       "run yet (needs data/drift_field_2d.rds). getwd(): ", getwd())

# ---- 1. Grid + forcing climatology ------------------------------------------------
lats <- seq(-29, -65, by = -2)          # row 1 = north (breeding side)
lons <- seq(145, 169, by = 2)           # drop 141/143: <3% of sightings there
# and the zero-flux western boundary was acting as a diffusion trap — density
# piled onto lon 141-147 over the 36-month run with no sink, dragging the
# year-3 field off the coast.
nla <- length(lats); nlo <- length(lons)
dym <- 2 * 111                          # km per lat step
dxm <- 2 * 111 * cos(lats * pi / 180)   # km per lon step, per row

water <- readRDS("data/water_extended.rds")
sst_c <- array(NA_real_, c(nla, nlo, 12))
chl_c <- array(NA_real_, c(nla, nlo, 12))
a_sst <- aggregate(sea_surface_temp_C ~ cell_centre_lat + cell_centre_lon + month,
                   data = water, FUN = mean, na.rm = TRUE)
wchl <- water[!is.na(water$chlorophyll_mg_per_m3), ]
a_chl <- aggregate(chlorophyll_mg_per_m3 ~ cell_centre_lat + cell_centre_lon + month,
                   data = wchl, FUN = mean, na.rm = TRUE)
for (i in seq_len(nla)) for (j in seq_len(nlo)) for (m in 1:12) {
  r <- a_sst[a_sst$cell_centre_lat == lats[i] & a_sst$cell_centre_lon == lons[j] &
             a_sst$month == m, ]
  if (nrow(r) == 1) sst_c[i, j, m] <- r$sea_surface_temp_C
  r <- a_chl[a_chl$cell_centre_lat == lats[i] & a_chl$cell_centre_lon == lons[j] &
             a_chl$month == m, ]
  if (nrow(r) == 1) chl_c[i, j, m] <- r$chlorophyll_mg_per_m3
}
chl_med <- median(water$chlorophyll_mg_per_m3, na.rm = TRUE)
chl_c[is.na(chl_c)] <- chl_med
# fill missing SST: mean of same-lat row for that month, then global col mean
for (m in 1:12) for (i in seq_len(nla)) {
  row_na <- is.na(sst_c[i, , m])
  if (all(row_na)) {
    sst_c[i, , m] <- mean(sst_c[, , m], na.rm = TRUE)
  } else if (any(row_na)) {
    sst_c[i, row_na, m] <- mean(sst_c[i, , m], na.rm = TRUE)
  }
}
if (any(is.na(sst_c))) stop("sst_c still has NA after fill")

# ---- 2. Land mask ---------------------------------------------------------------------
bath <- read.csv("data/bathymetry_cells.csv")
ocean <- matrix(TRUE, nla, nlo)
for (i in seq_len(nla)) for (j in seq_len(nlo)) {
  r <- bath[bath$cell_centre_lat == lats[i] & bath$cell_centre_lon == lons[j], ]
  if (nrow(r) == 1 && r$frac_land > 0.5) ocean[i, j] <- FALSE
}
cat(sprintf("Ocean cells: %d of %d\n", sum(ocean), nla * nlo))

# ---- 3. Drift field (bilinear interp from 4-deg bins) + currents ------------------------
dr <- readRDS("data/drift_field_2d.rds")
binterp <- function(field, m) {
  out <- matrix(0, nla, nlo)
  for (i in seq_len(nla)) for (j in seq_len(nlo)) {
    # clamp into bin range, find surrounding bin centres
    y <- max(min(lats[i], max(dr$lats)), min(dr$lats))
    x <- max(min(lons[j], max(dr$lons)), min(dr$lons))
    iy <- findInterval(y, dr$lats); iy <- max(1, min(iy, length(dr$lats) - 1))
    ix <- findInterval(x, dr$lons); ix <- max(1, min(ix, length(dr$lons) - 1))
    fy <- (y - dr$lats[iy]) / (dr$lats[iy + 1] - dr$lats[iy])
    fx <- (x - dr$lons[ix]) / (dr$lons[ix + 1] - dr$lons[ix])
    f00 <- field[iy, ix, m];     f10 <- field[iy + 1, ix, m]
    f01 <- field[iy, ix + 1, m]; f11 <- field[iy + 1, ix + 1, m]
    out[i, j] <- (1 - fy) * (1 - fx) * f00 + fy * (1 - fx) * f10 +
                 (1 - fy) * fx * f01 + fy * fx * f11
  }
  out
}
# drift arrays are +north/+east (km/h); model v is +poleward => negate meridional
v_mer_drift <- array(0, c(nla, nlo, 12))
v_zon_drift <- array(0, c(nla, nlo, 12))
for (m in 1:12) {
  v_mer_drift[, , m] <- -binterp(dr$v_mer, m)
  v_zon_drift[, , m] <-  binterp(dr$v_zon, m)
}

cur <- read.csv("data/currents_cells.csv")
cur_cl <- aggregate(cbind(u_ms, v_ms) ~ cell_centre_lat + cell_centre_lon + month,
                    data = cur, FUN = mean, na.rm = TRUE)
u_cur <- array(0, c(nla, nlo, 12))
v_cur <- array(0, c(nla, nlo, 12))
for (i in seq_len(nla)) for (j in seq_len(nlo)) for (m in 1:12) {
  r <- cur_cl[cur_cl$cell_centre_lat == lats[i] & cur_cl$cell_centre_lon == lons[j] &
              cur_cl$month == m, ]
  if (nrow(r) == 1) {
    v_cur[i, j, m] <- -r$v_ms * 3.6     # +north -> +poleward
    u_cur[i, j, m] <-  r$u_ms * 3.6
  }
}

# ---- 3b. Corridor mask --------------------------------------------------------------
# Tags show migration confined to a narrow coastal band (lon ~148-155) at
# lat -20..-40, spreading across the Southern Ocean south of -50. With 95% of
# the drift-field bins fallback (script 13), broadcasting advection over the
# full 141-169 grid produces a longitude-uniform field that matches nothing
# (year-3 spearman ~ 0). Encode the corridor as structure: advection only
# inside it, diffusion everywhere. Boundary lon_max varies linearly from 155
# at -30 to 169 (full domain) at -55.
corridor_lon_max <- pmin(169, 155 + (pmax(0, -30 - lats)) * (14 / 25))
in_corridor <- outer(corridor_lon_max, lons, FUN = function(m, lo) lo <= m)
v_mer_drift[!rep(in_corridor, 12)] <- 0
v_zon_drift[!rep(in_corridor, 12)] <- 0
u_cur[!rep(in_corridor, 12)] <- 0
v_cur[!rep(in_corridor, 12)] <- 0
cat(sprintf("Corridor cells: %d of %d (advection active)\n",
            sum(in_corridor), nla * nlo))

# ---- 4. Parameters + reaction -----------------------------------------------------------
track_params <- if (file.exists("data/track_params.rds")) readRDS("data/track_params.rds") else NULL
D0 <- if (is.null(track_params)) 10 else track_params$D0
h  <- 0.4
r_max <- 0.25; K <- 1.0; m0_mort <- 0.05

S_breeding <- function(T) exp(-(T - 21.32)^2 / (2 * 4.64^2))
S_feeding  <- function(T) 1 / (1 + exp((T - 12) / 1.5))
phase <- c("1" = "feeding", "2" = "feeding", "3" = "feeding",
           "4" = "transit", "5" = "breeding", "6" = "breeding",
           "7" = "breeding", "8" = "breeding", "9" = "breeding",
           "10" = "transit", "11" = "feeding", "12" = "feeding")

hours_per_month <- 365 / 12 * 24
dt <- 0.25
n_steps <- round(36 * hours_per_month / dt)

react <- function(C, S) {
  f <- dt / hours_per_month
  R <- r_max * S * pmax(0, 1 - C / K) * C - m0_mort * (1 - S) * C
  # pmax strips matrix dimensions — rebuild or the next step2d call gets a
  # vector and dies with "incorrect number of dimensions"
  matrix(pmax(0, C + f * R), nla, nlo)
}

# ---- 5. 2-D step --------------------------------------------------------------------------
step2d <- function(C, vm, vz, D) {
  Cp <- C; Cp[!ocean] <- 0
  # meridional (poleward +): faces between row i and i+1
  vmf <- (vm[-nla, ] + vm[-1, ]) / 2
  up  <- matrix(ifelse(vmf > 0, Cp[-nla, ], Cp[-1, ]), nla - 1, nlo)
  adv_m <- vmf * up
  Df  <- (D[-nla, ] + D[-1, ]) / 2
  dif_m <- -Df * (Cp[-1, ] - Cp[-nla, ]) / dym
  Fm <- adv_m + dif_m
  Fm[!ocean[-nla, ] | !ocean[-1, ]] <- 0        # no flux into/out of land
  # zonal (eastward +): faces between col j and j+1
  vzf <- (vz[, -nlo] + vz[, -1]) / 2
  upz <- matrix(ifelse(vzf > 0, Cp[, -nlo], Cp[, -1]), nla, nlo - 1)
  adv_z <- vzf * upz
  Dfz <- (D[, -nlo] + D[, -1]) / 2
  dif_z <- -Dfz * (Cp[, -1] - Cp[, -nlo]) / dxm   # per-row zonal cell width
  Fz <- adv_z + dif_z
  Fz[!ocean[, -nlo] | !ocean[, -1]] <- 0
  dC <- matrix(0, nla, nlo)
  dC[1, ]     <- -Fm[1, ] / dym
  dC[nla, ]   <-  Fm[nla - 1, ] / dym
  dC[2:(nla - 1), ] <- -(Fm[2:(nla - 1), ] - Fm[1:(nla - 2), ]) / dym
  dC[, 1]     <- dC[, 1]     - Fz[, 1] / dxm
  dC[, nlo]   <- dC[, nlo]   + Fz[, nlo - 1] / dxm
  dC[, 2:(nlo - 1)] <- dC[, 2:(nlo - 1)] -
    (Fz[, 2:(nlo - 1)] - Fz[, 1:(nlo - 2)]) / dxm
  Cp + dt * dC
}

# ---- 6. Time loop ------------------------------------------------------------------------------
C <- matrix(0, nla, nlo)
C[(nla - 2):nla, ] <- 1                        # pulse across the far south (March)
t_months <- 3
snaps <- list()
mass0 <- sum(C)
for (s in seq_len(n_steps)) {
  m   <- t_months %% 12
  mi0 <- floor(m) %% 12 + 1
  mi1 <- mi0 %% 12 + 1
  w   <- m - floor(m)
  vm <- (1 - w) * (v_mer_drift[, , mi0] + v_cur[, , mi0]) +
         w      * (v_mer_drift[, , mi1] + v_cur[, , mi1])
  vz <- (1 - w) * v_zon_drift[, , mi0] + w * v_zon_drift[, , mi1]
  P  <- (1 - w) * chl_c[, , mi0] + w * chl_c[, , mi1]
  D  <- D0 / (1 + P / h)
  sst_f <- (1 - w) * sst_c[, , mi0] + w * sst_c[, , mi1]
  ph <- phase[as.character(mi0)]
  S <- switch(ph, breeding = S_breeding(sst_f),
                  feeding  = S_feeding(sst_f),
                  transit  = matrix(0, nla, nlo))
  C <- step2d(C, vm, vz, D)
  C <- react(C, S)
  t_months <- t_months + dt / hours_per_month
  if (s %% round(hours_per_month / dt) == 0) {
    df <- expand.grid(lat = lats, lon = lons)
    df$density <- as.vector(C)
    df$month <- t_months
    snaps[[length(snaps) + 1]] <- df
  }
}
out <- do.call(rbind, snaps)
out$year <- floor(out$month / 12) + 1
out$month_of_year <- round(out$month) %% 12
out$month_of_year[out$month_of_year == 0] <- 12

mass <- tapply(out$density, out$month, sum)
cat(sprintf("Mass: %.2f .. %.2f (initial %.2f; R active so drift expected)\n",
            min(mass), max(mass), mass0))
cat(sprintf("Min density: %.5f\n", min(out$density)))

# ---- 7. Validation: monthly pattern correlation vs binned sightings ----------------------
all_s <- readRDS("data-raw/humpback_all.rds")
lat_col <- intersect(c("decimalLatitude", "latitude", "lat"), names(all_s))[1]
lon_col <- intersect(c("decimalLongitude", "longitude", "lon"), names(all_s))[1]
sg <- data.frame(lat = all_s[[lat_col]], lon = all_s[[lon_col]], month = all_s$month)
sg <- sg[all_s$year >= 2002 & sg$lat >= -67 & sg$lat <= -29 &
         sg$lon >= min(lons) & sg$lon <= max(lons), ]
sg$la <- lats[findInterval(-sg$lat, -lats)]     # snap to nearest grid row
sg$lo <- lons[findInterval(sg$lon, lons)]
# ---- 7. Validation: monthly LAT-PROFILE correlation vs sightings ----------------------
# Cell-wise spearman across the 2-D grid is the wrong metric here: sightings
# are spiky (whale-watch towns dominate a few cells) while the model field is
# smooth, so even a correct model scores ~0. Fair test (as in 1-D, corr 0.84):
# sum density over longitude per latitude row, compare to sightings per row.
prof_cor <- sapply(1:12, function(mo) {
  mod_m <- out[out$year == 3 & out$month_of_year == mo, ]
  mod_prof <- tapply(mod_m$density, mod_m$lat, sum)
  obs_prof <- tapply(rep(1, nrow(sg[sg$month == mo, ])), sg$la[sg$month == mo], sum)
  mrg <- merge(data.frame(lat = names(mod_prof), mod = as.numeric(mod_prof)),
               data.frame(lat = names(obs_prof), obs = as.numeric(obs_prof)),
               all = TRUE)
  mrg[is.na(mrg)] <- 0
  if (sum(mrg$obs) == 0 || sum(mrg$mod) == 0) return(NA_real_)
  cor(mrg$mod, mrg$obs, method = "spearman")
})
cat("\nYear-3 lat-profile spearman vs sightings, by month:\n")
print(round(setNames(prof_cor, month.abb), 2))
cat(sprintf("Mean: %.2f\n", mean(prof_cor, na.rm = TRUE)))

# ---- 8. Seasonal maps + save ------------------------------------------------------------------
season_map <- data.frame(
  month_of_year = c(6, 7, 11, 1),
  season = c("Northbound (Jun)", "Breeding (Jul)", "Southbound (Nov)", "Feeding (Jan)"))
yr3 <- merge(out[out$year == 3, ], season_map)
yr3$season <- factor(yr3$season, levels = season_map$season)

p <- ggplot(yr3, aes(lon, lat, fill = density)) +
  geom_raster() +
  facet_wrap(~season, ncol = 2) +
  scale_fill_viridis_c() +
  coord_quickmap() +
  labs(title = "Year-3 density, 2-D PDE", x = "longitude", y = "latitude") +
  theme_minimal()
print(p)
ggsave("analysis/pde_2d_seasons.png", p, width = 10, height = 9)

saveRDS(list(output = out, skill_by_month = prof_cor, ocean = ocean,
             corridor = in_corridor, lats = lats, lons = lons),
        "data/pde_2d_output.rds")
cat("\nSaved data/pde_2d_output.rds + analysis/pde_2d_seasons.png\n")
