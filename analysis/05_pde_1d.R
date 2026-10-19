# ============================================================================
# 05_pde_1d.R
# 1-D along-coast advection-diffusion PDE for east-coast humpbacks.
#
#   dC/dt = d/dx( D(x,t) dC/dx ) - d/dx( v(x,t) C )
#
# x = latitude along a coastal strip (2-degree water-grid cell centres,
# odd integers, QLD at index 1 down to the Southern Ocean).
#
# Terms:
#   v(x,t)  seasonal migration drift: direction set by the migration calendar
#           (northbound May-Aug toward breeding grounds, southbound Oct-Feb
#           toward feeding grounds); magnitude = v_max * tanh(alpha * |dS/dx|)
#           with S the fitted seasonal thermal suitability envelope.
#   D(x,t)  Area-Restricted-Search diffusion: D0 / (1 + P/h), P = Chl-a.
#
# Prototyped in Python (analysis/pde_proto.py) before this R port; key lessons
# baked in:
#   * velocity must be sign-consistent along the corridor per season,
#     otherwise the wave tears apart where dS/dx changes sign;
#   * alpha must be O(1e3) because dS/dx ~ 1e-3 per km on a 2-degree grid;
#   * boundaries are zero-flux (advective flux included) to conserve mass.
#
# Run: setwd() to the repo root, then source("analysis/05_pde_1d.R")
# Packages: none beyond base R + ggplot2 (install.packages("ggplot2") once).
# ============================================================================

library(ggplot2)

# ---- 1. Coastal-strip forcing climatology ---------------------------------
water <- readRDS("data/water_extended.rds")
strip <- water[water$cell_centre_lon >= 145 & water$cell_centre_lon <= 156, ]

grid <- seq(-29, -55, by = -2)          # decreasing latitude = poleward index
N <- length(grid)
dx_km <- 2 * 111                        # ~222 km between cell centres

clim <- aggregate(cbind(sea_surface_temp_C, chlorophyll_mg_per_m3) ~
                    cell_centre_lat + month, data = strip, FUN = mean, na.rm = TRUE)

chl_fill <- median(clim$chlorophyll_mg_per_m3, na.rm = TRUE)  # background for
# southern cells where the merged product has no Chl-a

sst <- matrix(NA_real_, N, 12)
chl <- matrix(NA_real_, N, 12)
for (i in seq_len(N)) {
  for (m in 1:12) {
    row <- clim[clim$cell_centre_lat == grid[i] & clim$month == m, ]
    if (nrow(row) == 1) {
      sst[i, m] <- ifelse(is.na(row$sea_surface_temp_C),
                          NA_real_, row$sea_surface_temp_C)
      chl[i, m] <- ifelse(is.na(row$chlorophyll_mg_per_m3),
                          chl_fill, row$chlorophyll_mg_per_m3)
    } else {
      chl[i, m] <- chl_fill
    }
  }
}

# impute any missing SST cells: interpolate along latitude within each month,
# then along month within each latitude as a fallback — NA must not enter the
# PDE loop or it poisons S -> v -> C.
for (m in 1:12) {
  col <- sst[, m]
  if (any(is.na(col)))
    sst[, m] <- approx(seq_len(N), col, xout = seq_len(N),
                       rule = 2, ties = mean)$y
}
for (i in seq_len(N)) {
  rowv <- sst[i, ]
  if (any(is.na(rowv)))
    sst[i, ] <- approx(1:12, rowv, xout = 1:12, rule = 2, ties = mean)$y
}
if (any(is.na(sst))) stop("sst still has NA after imputation — check the climatology")
chl[is.na(chl)] <- chl_fill

# ---- 2. Suitability envelopes (from script 02 GAM fits) -------------------
# breeding: Gaussian peak at 21.3 C, sigma 4.6
S_breeding <- function(T) exp(-(T - 21.32)^2 / (2 * 4.64^2))
# feeding: saturating cold-plateau (~1 below ~12 C, 0 above ~17 C)
S_feeding  <- function(T) 1 / (1 + exp((T - 12) / 1.5))

S_b <- S_breeding(sst)
S_f <- S_feeding(sst)

dSdx <- function(S) {
  g <- numeric(length(S))
  n <- length(S)
  g[2:(n - 1)] <- (S[3:n] - S[1:(n - 2)]) / (2 * dx_km)
  g[1] <- (S[2] - S[1]) / dx_km
  g[n] <- (S[n] - S[n - 1]) / dx_km
  g
}

# ---- 3. Parameters ---------------------------------------------------------
# Movement parameters from the E1 tag tracks (analysis/06_track_parameters.R,
# Andrews-Goff et al. 2023, Movebank doi:10.5441/001/1.294).
track_params <- if (file.exists("data/track_params.rds")) readRDS("data/track_params.rds") else NULL
v_max <- if (is.null(track_params)) 4.9 else track_params$v_max   # km/h directed drift
alpha <- 3000       # tanh gain on |dS/dx| (dS/dx ~ 1e-3 per km on this grid)
D0    <- if (is.null(track_params)) 10 else track_params$D0       # km^2/h; resident-phase diffusivity
h     <- 0.4        # Holling half-saturation (mg/m3 Chl-a)

hours_per_month <- 365 / 12 * 24
dt <- 0.25          # hours (CFL: v*dt/dx << 1)
total_months <- 36
n_steps <- round(total_months * hours_per_month / dt)

# ---- 4. One explicit step: upwind advection + centred diffusion ------------
pde_step <- function(C, v, D) {
  vf <- 0.5 * (v[-N] + v[-1])                    # face velocities
  up <- ifelse(vf > 0, C[-N], C[-1])             # upwind density at each face
  adv_flux <- vf * up
  Df <- 0.5 * (D[-N] + D[-1])
  dif_flux <- -Df * (C[-1] - C[-N]) / dx_km
  net <- adv_flux + dif_flux                     # poleward flux across faces
  # zero-flux boundaries: end faces carry no flux, so boundary cells keep the
  # one interior face they own. Freezing dC at the boundary cells (previous
  # version) leaked ~2/3 of total mass through those faces over 3 years.
  dC <- numeric(N)
  dC[1] <- -net[1] / dx_km
  dC[N] <-  net[N - 1] / dx_km
  dC[2:(N - 1)] <- -(net[2:(N - 1)] - net[1:(N - 2)]) / dx_km
  C + dt * dC
}

# ---- 5. Time loop (monthly forcing, linear interpolation) ------------------
# Advection: uniform seasonal drift from the E1 tag tracks (median meridional
# speed by month). Prototyping (analysis/pde_proto.py) showed a tanh(alpha*|dS/dx|)
# spatial modulation is counterproductive on this coarse grid — it pins the
# migration front against the southern boundary each June because dS/dx fades
# there. The tags themselves show whales crossing thermally-flat water at
# 3-7 km/h during migration, i.e. the drive is seasonal/hormonal, not a local
# gradient response. Uniform drift + ARS diffusion reproduces the observed
# traverse (-43 Jun-Aug / -55 Dec-Mar) and is the version to build on.
# The suitability envelopes still enter via D (ARS) and, later, R(C,T).
drift_obs <- c("1" = -1.39, "2" = -1.06, "3" = -0.88, "4" = -0.87,
               "5" =  1.06, "6" =  4.71, "7" =  1.90, "8" =  1.99,
               "9" =  0.00, "10" = -2.55, "11" = -2.55, "12" = -2.87)
# drift_obs: km/h, + = northward. Grid velocity + = poleward, so v = -drift.

C <- numeric(N); C[(N-2):N] <- 1      # initial pulse on the southern feeding
                                      # grounds (March, ~-53 to -55)
t_months <- 3                         # start in March (southern limit, per tags)
snaps <- list()

for (s in seq_len(n_steps)) {
  m  <- t_months %% 12
  m0 <- floor(m) %% 12 + 1            # R 1-indexed month
  m1 <- m0 %% 12 + 1
  w  <- m - floor(m)

  drift <- (1 - w) * drift_obs[as.character(m0)] +
           w       * drift_obs[as.character(m1)]   # km/h, + = northward
  v <- rep(-drift, N)                              # uniform along corridor

  P <- (1 - w) * chl[, m0] + w * chl[, m1]
  D <- D0 / (1 + P / h)

  C <- pde_step(C, v, D)
  t_months <- t_months + dt / hours_per_month

  if (s %% round(hours_per_month / dt) == 0)          # monthly snapshot
    snaps[[length(snaps) + 1]] <- data.frame(month = t_months, lat = grid,
                                             density = C)
}

out <- do.call(rbind, snaps)
out$year  <- floor(out$month / 12) + 1
out$month_of_year <- out$month %% 12

# ---- 6. Checks + plots ------------------------------------------------------
mass <- tapply(out$density, out$month, sum)
cat(sprintf("Mass conservation: sum(C) ranges %.3f .. %.3f (initial %.3f)\n",
            min(mass), max(mass), 3.0))
cat(sprintf("Min density: %.5f (should be >= 0)\n", min(out$density)))

yr3 <- out[out$year == 3, ]
centroid <- tapply(seq_len(nrow(yr3)), yr3$month_of_year, function(idx) {
  d <- yr3$density[idx]; la <- yr3$lat[idx]
  if (sum(d) <= 0) return(NA_real_)
  sum(la * d) / sum(d)
})
cat("\nYear-3 monthly density centroid (latitude):\n")
print(round(centroid, 1))

p1 <- ggplot(out[out$year == 3, ],
             aes(x = lat, y = density, group = month_of_year,
                 colour = factor(month_of_year))) +
  geom_line() +
  scale_colour_viridis_d(name = "month") +
  labs(title = "Year 3: whale density along the east-coast corridor",
       x = "latitude", y = "density (relative)") +
  theme_minimal()
print(p1)
ggsave("analysis/pde_1d_year3_monthly.png", p1, width = 8, height = 5)

p2 <- ggplot(out, aes(x = month, y = lat, fill = density)) +
  geom_raster() +
  scale_fill_viridis_c() +
  labs(title = "Hovmoller: 3-year along-coast density",
       x = "months since start", y = "latitude") +
  theme_minimal()
print(p2)
ggsave("analysis/pde_1d_hovmoller.png", p2, width = 8, height = 5)

saveRDS(out, "data/pde_1d_output.rds")
cat("\nSaved data/pde_1d_output.rds and two PNGs in analysis/\n")
