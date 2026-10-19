# ============================================================================
# 11_thermal_departure.R
# Script 10, but migration timing is CUED BY TEMPERATURE instead of being a
# fixed calendar clock. This is the model version that can actually respond
# to climate scenarios: as the ocean warms, the thermal cue fires earlier,
# shifting migration timing — an effect the fixed drift curve cannot express.
#
# Mechanism (hypothesis, NOT fitted — see note below):
#   The seasonal drift curve drift(month) from the tags defines the SHAPE of
#   the migration clock. The phase of that clock shifts with the thermal
#   state of the corridor:
#
#     phase_shift (months) = beta * T_anom
#
#   where T_anom = SST averaged over the departure region minus its 2002-2015
#   climatological mean, in °C; beta (months/°C) is the cue sensitivity.
#
#   Southbound departure region: Coral Sea / breeding grounds (lat -25..-15),
#   applied Sep-Jan. Northbound departure region: feeding grounds (lat -65..
#   -55), applied Mar-Aug. The drift evaluated at time t uses the curve at
#   month (t + shift), i.e. positive shift = earlier migration in a warm year.
#
#   beta is a LITERATURE-STYLE sensitivity parameter, not estimable from the
#   E1 tags (only 4 seasons, 2008-2011 — cannot fit an interannual timing-SST
#   relationship). Cetacean literature values for thermal migration cues are
#   order 2-4 weeks per °C; we default beta = 0.75 months/°C (~3 weeks) and
#   the script prints a sensitivity sweep over beta = 0, 0.25, 0.5, 0.75, 1.0
#   plus a +2C scenario run to show the climate response.
#
# Everything else (currents, ARS diffusion, R(C,T), mass-conserving scheme,
# drift shape, validation vs sightings) is identical to script 10.
#
# Run: setwd() to repo root, then source("analysis/11_thermal_departure.R")
# ============================================================================

library(ggplot2)

if (!file.exists("data/cell_month_2002_2015.csv"))
  stop("Working directory is not the repo root. Set it via Files pane -> ",
       "More -> Set As Working Directory. Current getwd(): ", getwd())

# ---- 1. Interannual forcing (as 08) + climatology for anomalies --------------
cm  <- read.csv("data/cell_month_2002_2015.csv")
cur <- read.csv("data/currents_cells.csv")
stopifnot(nrow(cm) > 0, nrow(cur) > 0)

grid <- seq(-29, -65, by = -2)
N <- length(grid)
dx_km <- 2 * 111

chl_med <- median(cm$chlorophyll_mg_per_m3, na.rm = TRUE)

get_env <- function(year, month) {
  rows <- cm[cm$year == year & cm$month == month &
               cm$cell_centre_lat %in% grid &
               cm$cell_centre_lon >= 145 & cm$cell_centre_lon <= 156, ]
  if (nrow(rows) == 0)
    stop("get_env(", year, ", ", month, "): zero rows from cm")
  a_sst <- aggregate(sea_surface_temp_C ~ cell_centre_lat, data = rows,
                     FUN = mean, na.rm = TRUE)
  rows_chl <- rows[!is.na(rows$chlorophyll_mg_per_m3), ]
  a_chl <- if (nrow(rows_chl) > 0)
    aggregate(chlorophyll_mg_per_m3 ~ cell_centre_lat, data = rows_chl,
              FUN = mean, na.rm = TRUE)
  else data.frame(cell_centre_lat = grid, chlorophyll_mg_per_m3 = NA_real_)
  sst <- chl <- rep(NA_real_, N)
  for (i in seq_len(N)) {
    j <- which(a_sst$cell_centre_lat == grid[i])
    if (length(j) == 1) sst[i] <- a_sst$sea_surface_temp_C[j]
    j <- which(a_chl$cell_centre_lat == grid[i])
    if (length(j) == 1) chl[i] <- a_chl$chlorophyll_mg_per_m3[j]
  }
  chl[is.na(chl)] <- chl_med
  if (any(is.na(sst))) {
    ok <- !is.na(sst)
    sst <- approx(which(ok), sst[ok], xout = seq_len(N), rule = 2)$y
  }
  rows <- cur[cur$year == year & cur$month == month &
                cur$cell_centre_lat %in% grid &
                cur$cell_centre_lon >= 145 & cur$cell_centre_lon <= 156, ]
  a <- aggregate(v_ms ~ cell_centre_lat, data = rows, FUN = mean, na.rm = TRUE)
  v <- rep(0, N)
  for (i in seq_len(N)) {
    j <- which(a$cell_centre_lat == grid[i])
    if (length(j) == 1) v[i] <- -a$v_ms[j] * 3.6
  }
  list(sst = sst, chl = chl, vcur = v)
}

# climatological SST per (lat, month) for anomaly reference
clim_sst <- array(NA_real_, c(N, 12))
for (m in 1:12) {
  for (y in 2002:2015) {
    e <- get_env(y, m)
    clim_sst[, m] <- if (y == 2002) e$sst else clim_sst[, m] + e$sst
  }
}
clim_sst <- clim_sst / length(2002:2015)

# ---- 2. Model pieces (as 10) ---------------------------------------------------
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

# drift at fractional (possibly shifted) month; the curve is periodic 1..12
drift_at <- function(m_frac) {
  m_wrapped <- (m_frac - 1) %% 12 + 1        # wrap into [1, 12]
  m0 <- floor(m_wrapped)
  m1 <- m0 %% 12 + 1
  w  <- m_wrapped - m0
  (1 - w) * drift_obs[as.character(m0)] + w * drift_obs[as.character(m1)]
}

S_breeding <- function(T) exp(-(T - 21.32)^2 / (2 * 4.64^2))
S_feeding  <- function(T) 1 / (1 + exp((T - 12) / 1.5))
phase <- c("1" = "feeding", "2" = "feeding", "3" = "feeding",
           "4" = "transit", "5" = "breeding", "6" = "breeding",
           "7" = "breeding", "8" = "breeding", "9" = "breeding",
           "10" = "transit", "11" = "feeding", "12" = "feeding")

r_max <- 0.25; K <- 1.0; m0_mort <- 0.05
USE_WALL_LOSS <- FALSE; m_wall <- 0.20; WALL_CELLS <- 2

reaction_step <- function(C, S, dt_h) {
  f <- dt_h / hours_per_month
  R <- r_max * S * pmax(0, 1 - C / K) * C - m0_mort * (1 - S) * C
  if (USE_WALL_LOSS) R[(N - WALL_CELLS + 1):N] <-
    R[(N - WALL_CELLS + 1):N] - m_wall * C[(N - WALL_CELLS + 1):N]
  pmax(0, C + f * R)
}

hours_per_month <- 365 / 12 * 24
dt <- 1.0

# ---- 3. Thermal cue --------------------------------------------------------------
# departure-region indices
idx_sb <- which(grid >= -31)                 # Coral Sea proxy: grid starts at
# -29, so the ideal -25..-15 breeding-ground cue region is off-grid; the
# northernmost cells (-29, -31) are the closest available.
idx_nb <- which(grid <= -43)                 # feeding-grounds proxy: cm data
# stops at -45, so the true -65..-55 cue region is unavailable; -45/-43 is
# the nearest Southern-Ocean-influenced cell in the forcing table.

# shift in months for a given calendar month + SST field. Positive shift =
# migration phase advanced (warm year -> earlier departure).
MAX_SHIFT <- 0.5    # months. Thermal cues nudge timing; they don't rephase
# the calendar. Uncapped, a +2C scenario at beta=0.75 gives a 1.5-month phase
# jump, which is a huge velocity change in transition months (Apr->Jun drift
# goes -0.87 -> +4.71 km/h) and catapults the wave off the domain (dry-run
# produced non-physical +24 deg centroid shifts in May).
cue_shift <- function(sst_field, month, beta) {
  raw <- if (month %in% c(9, 10, 11, 12, 1)) {     # southbound window
    anom <- mean(sst_field[idx_sb] - clim_sst[idx_sb, month], na.rm = TRUE)
    beta * anom
  } else if (month %in% 3:8) {                      # northbound window
    anom <- mean(sst_field[idx_nb] - clim_sst[idx_nb, month], na.rm = TRUE)
    beta * anom
  } else 0
  max(-MAX_SHIFT, min(MAX_SHIFT, raw))
}

# ---- 4. Hindcast 2002-2015 for a given beta ----------------------------------------
run_model <- function(beta, sst_offset = 0, verbose = FALSE) {
  C <- numeric(N); C[(N - 2):N] <- 1
  cents <- data.frame(year = integer(), month = integer(), mod_lat = numeric())
  for (year in 2002:2015) {
    for (month in 1:12) {
      env <- get_env(year, month)
      sst_f <- env$sst + sst_offset
      shift <- cue_shift(sst_f, month, beta)
      if (!is.finite(shift)) stop("cue_shift returned NA at ", year, "-", month,
                                  " — check idx_sb/idx_nb are non-empty")
      v <- rep(-drift_at(month + shift), N) + env$vcur
      D <- D0 / (1 + env$chl / h)
      ph <- phase[as.character(month)]
      S <- switch(ph, breeding = S_breeding(sst_f),
                      feeding  = S_feeding(sst_f),
                      transit  = rep(0, N))
      for (hh in seq_len(round(hours_per_month / dt))) {
        C <- pde_step(C, v, D, dt)
        C <- reaction_step(C, S, dt)
      }
      if (sum(C) > 0)
        cents[nrow(cents) + 1, ] <- list(year, month, sum(grid * C) / sum(C))
    }
  }
  cents
}

# observed monthly centroids
all_s <- readRDS("data-raw/humpback_all.rds")
lat_col <- intersect(c("decimalLatitude", "latitude", "lat"), names(all_s))[1]
lon_col <- intersect(c("decimalLongitude", "longitude", "lon"), names(all_s))[1]
keep <- all_s$year >= 2002 & all_s$year <= 2015 &
  ((all_s[[lon_col]] >= 145 & all_s[[lon_col]] <= 156 & all_s[[lat_col]] >= -45) |
     (all_s[[lon_col]] >= 80  & all_s[[lon_col]] <= 180 & all_s[[lat_col]] < -45))
sight <- all_s[keep, ]
obs_ym <- aggregate(sight[[lat_col]], list(year = sight$year, month = sight$month),
                    FUN = mean)
names(obs_ym)[3] <- "obs_lat"

# ---- 5. Sensitivity sweep over beta -------------------------------------------------
cat("beta (mo/°C)   mean|err|   corr   Oct|err|  Nov|err|\n")
sweep <- data.frame()
for (b in c(0, 0.25, 0.5, 0.75, 1.0)) {
  cents <- run_model(b)
  cmp <- merge(cents, obs_ym)
  err <- cmp$mod_lat - cmp$obs_lat
  oct_e <- mean(abs(err[cmp$month == 10]))
  nov_e <- mean(abs(err[cmp$month == 11]))
  cat(sprintf("%6.2f %12.1f %8.2f %9.1f %9.1f\n",
              b, mean(abs(err)), cor(cmp$mod_lat, cmp$obs_lat), oct_e, nov_e))
  sweep <- rbind(sweep, data.frame(beta = b, abs_err = mean(abs(err)),
                                   corr = cor(cmp$mod_lat, cmp$obs_lat)))
}

# ---- 6. Climate scenario: +2C uniform warming, best-guess beta -------------------------
BETA <- 0.75
cat(sprintf("\nClimate scenario with beta = %.2f:\n", BETA))
base_run <- run_model(BETA, sst_offset = 0)
warm_run <- run_model(BETA, sst_offset = 2)
cmp <- merge(base_run, warm_run, by = c("year", "month"),
             suffixes = c("_base", "_warm"))
cmp$shift_deg <- cmp$mod_lat_warm - cmp$mod_lat_base
per_month <- aggregate(shift_deg ~ month, data = cmp, FUN = mean)
cat("Mean centroid shift under +2C (degrees lat, by month):\n")
print(round(per_month, 2))
cat(sprintf("Overall mean shift: %.2f deg lat (negative = further south)\n",
            mean(cmp$shift_deg)))

p <- ggplot(per_month, aes(month, shift_deg)) +
  geom_col(fill = "steelblue") +
  scale_x_continuous(breaks = 1:12) +
  labs(title = sprintf("+2C warming: centroid shift by month (beta = %.2f mo/°C)", BETA),
       x = "month", y = "shift (degrees latitude)") +
  theme_minimal()
print(p)
ggsave("analysis/thermal_departure_scenario.png", p, width = 8, height = 4)

saveRDS(list(sweep = sweep, baseline = base_run, warm = warm_run),
        "data/thermal_departure.rds")
cat("\nSaved data/thermal_departure.rds + analysis/thermal_departure_scenario.png\n")
