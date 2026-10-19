# ============================================================================
# 10_pde_reaction.R
# Script 09 + the demographic term R(C,T). Full equation, per framework PDF:
#
#   dC/dt = d/dx( D dC/dx ) - d/dx( v C ) + R(C, T)
#
# R is drafted here as three independently switchable mechanisms (flags in
# section 3). C is RELATIVE density (not absolute abundance), so R is written
# as fractional rates (per month) applied to C, not absolute births:
#
#   (1) growth on the grounds:  R+ = r_max * S(x, season) * (1 - C/K)
#       Suitability-weighted logistic gain. S = seasonal thermal envelope
#       (breeding peak 21.3 C / feeding cold-plateau), K = local carrying
#       density (in units of C, default 1 = local peak density). Gains accrue
#       where conditions suit the current life-history phase.
#
#   (2) mismatch mortality:     R- = -m0 * (1 - S) * C
#       Loss weighted by unsuitability: whales die off the grounds, not on
#       them. (A uniform -m0*C was tested first: it kills the population
#       within 2 years because the wave spends months in transit where S ~ 0,
#       so losses compound with no compensating gain.)
#
#   (3) boundary loss at the southern wall: extra mortality m_wall applied in
#       the southernmost WALL_CELLS cells. Represents whales dispersing
#       west/east out of the E. Australia stock corridor across the Southern
#       Ocean feeding grounds — fixes the unphysical summer pile-up at -65
#       seen in script 09 (model summer centroid -65 vs observed -58..-62).
#
# Calving note: mysticete calves are born on the breeding grounds in winter
# and do not migrate independently; with C a relative density, births are
# folded into mechanism (1) via the breeding envelope rather than modelled
# as a separate pulse.
#
# Run: setwd() to repo root, then source("analysis/10_pde_reaction.R")
# ============================================================================

library(ggplot2)

if (!file.exists("data/water_extended.rds"))
  stop("Working directory is not the repo root. Set it via Files pane -> ",
       "More -> Set As Working Directory. Current getwd(): ", getwd())

# ---- 1. Forcing climatology (as 09) -------------------------------------------
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
for (m in 1:12) {
  col <- sst[, m]
  if (any(is.na(col)))
    sst[, m] <- approx(seq_len(N), col, xout = seq_len(N), rule = 2)$y
}
if (any(is.na(sst))) stop("sst still has NA after imputation")

cur <- read.csv("data/currents_cells.csv")
cur <- cur[cur$cell_centre_lon >= 145 & cur$cell_centre_lon <= 156, ]
cur_clim <- aggregate(v_ms ~ cell_centre_lat + month, data = cur,
                      FUN = mean, na.rm = TRUE)
cur_kmh <- matrix(0, N, 12)
for (i in seq_len(N)) for (m in 1:12) {
  j <- which(cur_clim$cell_centre_lat == grid[i] & cur_clim$month == m)
  if (length(j) == 1) cur_kmh[i, m] <- -cur_clim$v_ms[j] * 3.6
}

# ---- 2. Suitability envelopes (from 02 GAM fits) -------------------------------
S_breeding <- function(T) exp(-(T - 21.32)^2 / (2 * 4.64^2))
S_feeding  <- function(T) 1 / (1 + exp((T - 12) / 1.5))
# life-history phase by month: breeding May-Sep, feeding Nov-Mar, transit else
phase <- c("1" = "feeding", "2" = "feeding", "3" = "feeding",
           "4" = "transit", "5" = "breeding", "6" = "breeding",
           "7" = "breeding", "8" = "breeding", "9" = "breeding",
           "10" = "transit", "11" = "feeding", "12" = "feeding")
S_b <- S_breeding(sst)
S_f <- S_feeding(sst)

# ---- 3. Parameters + switches --------------------------------------------------
track_params <- if (file.exists("data/track_params.rds")) readRDS("data/track_params.rds") else NULL
D0 <- if (is.null(track_params)) 10 else track_params$D0
h  <- 0.4

USE_GROWTH    <- TRUE     # mechanism (1)
USE_MORTALITY <- TRUE     # mechanism (2)
USE_WALL_LOSS <- FALSE     # mechanism (3)

r_max     <- 0.25   # max fractional gain per month at S = 1, C << K
K         <- 1.0    # carrying density, units of C (local peak)
m0        <- 0.05   # mismatch mortality: max fractional loss per month at S = 0
m_wall    <- 0.20   # extra fractional loss per month in wall cells
WALL_CELLS <- 2     # southernmost cells (-63, -65) treated as open boundary

# NOTE from dry-run testing: R(C,T) sets total mass / population trend but
# barely moves the monthly centroid (drift + diffusion own the timing). It is
# the right term for the climate-change question (trend under warming), not a
# fix for the Oct/Nov southbound lag — that lives in the drift curve and the
# missing second spatial dimension.

hours_per_month <- 365 / 12 * 24
dt <- 0.5
n_steps <- round(36 * hours_per_month / dt)

# ---- 4. PDE step + reaction -----------------------------------------------------
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

reaction_step <- function(C, S, dt_h) {
  # rates are per month; dt_h is hours -> scale
  f <- dt_h / hours_per_month
  R <- numeric(N)
  if (USE_GROWTH)    R <- R + r_max * S * pmax(0, 1 - C / K) * C
  if (USE_MORTALITY) R <- R - m0 * (1 - S) * C
  if (USE_WALL_LOSS) R[(N - WALL_CELLS + 1):N] <-
    R[(N - WALL_CELLS + 1):N] - m_wall * C[(N - WALL_CELLS + 1):N]
  pmax(0, C + f * R)     # positivity clamp: losses can't take C below 0
}

drift_obs <- c("1" = -1.39, "2" = -1.06, "3" = -0.88, "4" = -0.87,
               "5" =  1.06, "6" =  4.71, "7" =  1.90, "8" =  1.99,
               "9" =  0.00, "10" = -2.55, "11" = -2.55, "12" = -2.87)

# ---- 5. Time loop -----------------------------------------------------------------
C <- numeric(N); C[(N - 2):N] <- 1
t_months <- 3
snaps <- list()
for (s in seq_len(n_steps)) {
  m   <- t_months %% 12
  mi0 <- floor(m) %% 12 + 1          # month indices; named mi0/mi1, NOT m0 —
  mi1 <- mi0 %% 12 + 1               # m0 is the mortality parameter above and
  w   <- m - floor(m)                # a name collision here silently set m0=3

  drift <- (1 - w) * drift_obs[as.character(mi0)] + w * drift_obs[as.character(mi1)]
  v <- rep(-drift, N) + (1 - w) * cur_kmh[, mi0] + w * cur_kmh[, mi1]
  D <- D0 / (1 + ((1 - w) * chl[, mi0] + w * chl[, mi1]) / h)

  # reaction uses the suitability of the dominant phase this sub-interval
  ph <- phase[as.character(mi0)]
  S <- switch(ph, breeding = (1 - w) * S_b[, mi0] + w * S_b[, mi1],
                  feeding  = (1 - w) * S_f[, mi0] + w * S_f[, mi1],
                  transit  = rep(0, N))

  C <- pde_step(C, v, D, dt)
  C <- reaction_step(C, S, dt)
  t_months <- t_months + dt / hours_per_month

  if (s %% round(hours_per_month / dt) == 0)
    snaps[[length(snaps) + 1]] <- data.frame(month = t_months, lat = grid, density = C)
}
out <- do.call(rbind, snaps)
out$year <- floor(out$month / 12) + 1
out$month_of_year <- round(out$month) %% 12
out$month_of_year[out$month_of_year == 0] <- 12

mass <- tapply(out$density, out$month, sum)
cat(sprintf("Mass: %.3f .. %.3f (initial 3.000; growth/loss expected with R)\n",
            min(mass), max(mass)))
cat(sprintf("Min density: %.5f\n", min(out$density)))

# ---- 6. Validation vs sightings (as 09) ----------------------------------------------
all_s <- readRDS("data-raw/humpback_all.rds")
lat_col <- intersect(c("decimalLatitude", "latitude", "lat"), names(all_s))[1]
lon_col <- intersect(c("decimalLongitude", "longitude", "lon"), names(all_s))[1]
keep <- all_s$year >= 2002 &
  ((all_s[[lon_col]] >= 145 & all_s[[lon_col]] <= 156 & all_s[[lat_col]] >= -45) |
     (all_s[[lon_col]] >= 80  & all_s[[lon_col]] <= 180 & all_s[[lat_col]] < -45))
sight <- all_s[keep, ]
obs_cent <- tapply(sight[[lat_col]], sight$month, mean)

yr3 <- out[out$year == 3, ]
mod_cent <- tapply(seq_len(nrow(yr3)), yr3$month_of_year, function(idx) {
  d <- yr3$density[idx]
  if (sum(d) <= 0) return(NA_real_)
  sum(yr3$lat[idx] * d) / sum(d)
})
cmp <- data.frame(month = 1:12,
                  model = as.numeric(mod_cent[as.character(1:12)]),
                  observed = as.numeric(obs_cent[as.character(1:12)]))
cmp$err <- cmp$model - cmp$observed
cat("\nMonth  model  observed  err (deg lat)\n")
print(round(cmp, 1))
cat(sprintf("\nMean err %.1f, mean |err| %.1f, corr %.2f\n",
            mean(cmp$err, na.rm = TRUE), mean(abs(cmp$err), na.rm = TRUE),
            cor(cmp$model, cmp$observed, use = "complete.obs")))

p1 <- ggplot(yr3, aes(lat, density, group = month_of_year, colour = factor(month_of_year))) +
  geom_line() + scale_colour_viridis_d(name = "month") +
  labs(title = "Year 3: density with R(C,T)", x = "latitude", y = "density (relative)") +
  theme_minimal()
print(p1)
ggsave("analysis/pde_reaction_year3.png", p1, width = 8, height = 5)

p2 <- ggplot(out, aes(month, lat, fill = density)) +
  geom_raster() + scale_fill_viridis_c() +
  labs(title = "Hovmoller: PDE with R(C,T)", x = "months since start", y = "latitude") +
  theme_minimal()
print(p2)
ggsave("analysis/pde_reaction_hovmoller.png", p2, width = 8, height = 5)

p3 <- ggplot(cmp, aes(month)) +
  geom_line(aes(y = model, colour = "model"), linewidth = 1) +
  geom_line(aes(y = observed, colour = "observed"), linewidth = 1) +
  scale_x_continuous(breaks = 1:12) +
  labs(title = "Monthly centroid: model with R(C,T) vs sightings",
       y = "latitude", colour = NULL) +
  theme_minimal()
print(p3)
ggsave("analysis/pde_reaction_centroids.png", p3, width = 8, height = 4)

saveRDS(list(output = out, centroid_comparison = cmp,
             params = list(r_max = r_max, K = K, m0 = m0, m_wall = m_wall,
                           WALL_CELLS = WALL_CELLS)),
        "data/pde_reaction_output.rds")
cat("\nSaved data/pde_reaction_output.rds + three PNGs in analysis/\n")
