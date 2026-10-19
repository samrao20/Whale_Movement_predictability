# ============================================================================
# 08_hindcast_validation.R
# Hindcast the PDE 2002-2015 against observed monthly presence by latitude.
#
# Model: same equations as 07 (tag drift + IMOS currents + ARS diffusion),
# but forced with INTERANNUAL forcing per year-month, not climatology:
#   SST/Chl from data/cell_month_2002_2015.csv (monthly, per cell)
#   u/v from data/currents_cells.csv (2002-2015)
# Domain: -29..-45 (currents coverage). Initial condition: climatological
# March density from the year-3 spin-up of the 07 model.
#
# Validation: for each month 2002-2015, compare model density by latitude to
# observed presence (present == 1) by latitude, summarised by:
#   * centroid latitude error (model vs observed)
#   * skill = correlation between model density profile and observed presence
#     counts across latitudes
# plus Hovmoller panels side by side.
#
# Caveat: presence is effort-biased (whale-watch at -27, -33, -38) — the
# centroid comparison is about PHENOLOGY (timing of traverse), not absolute
# latitude. Report both honestly.
#
# Run: setwd() to repo root, then source("analysis/08_hindcast_validation.R")
# ============================================================================

library(ggplot2)

if (!file.exists("data/cell_month_2002_2015.csv"))
  stop("Working directory is not the repo root.\n",
       "  In RStudio: Files pane -> navigate into Whale_Movement_predictability ->\n",
       "  More -> Set As Working Directory, then re-source this script.\n",
       "  Current getwd(): ", getwd())

cm  <- read.csv("data/cell_month_2002_2015.csv")
cur <- read.csv("data/currents_cells.csv")
stopifnot(nrow(cm) > 0, nrow(cur) > 0)

grid <- seq(-29, -45, by = -2)
N <- length(grid)
dx_km <- 2 * 111

# index forcing into [lat, month] lookup matrices per year
get_env <- function(year, month) {
  rows <- cm[cm$year == year & cm$month == month &
               cm$cell_centre_lat %in% grid &
               cm$cell_centre_lon >= 145 & cm$cell_centre_lon <= 156, ]
  if (nrow(rows) == 0)
    stop("get_env(", year, ", ", month, "): zero rows from cm.\n",
         "  grid = ", paste(grid, collapse = " "), "\n",
         "  cm lats = ", paste(sort(unique(cm$cell_centre_lat)), collapse = " "),
         "\n  (stale `grid` from script 05 [-29..-55] is the usual cause — ",
         "rm(list=ls()) then re-source this script)")
  # aggregate SST and Chl separately: base aggregate() applies na.omit to the
  # whole cbind() response BEFORE na.rm, so a year-month where Chl is all-NA
  # would drop every row ("no rows to aggregate") even though SST is fine.
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
  chl[is.na(chl)] <- median(cm$chlorophyll_mg_per_m3, na.rm = TRUE)
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
    if (length(j) == 1) v[i] <- -a$v_ms[j] * 3.6      # + poleward km/h
  }
  list(chl = chl, vcur = v)
}

# ---- parameters + PDE step (as 07, with mass-conserving boundaries) --------
track_params <- if (file.exists("data/track_params.rds")) readRDS("data/track_params.rds") else NULL
D0 <- if (is.null(track_params)) 10 else track_params$D0
h  <- 0.4

pde_step <- function(C, v, D, dt) {
  vf <- 0.5 * (v[-N] + v[-1])
  up <- ifelse(vf > 0, C[-N], C[-1])
  adv_flux <- vf * up
  Df <- 0.5 * (D[-N] + D[-1])
  dif_flux <- -Df * (C[-1] - C[-N]) / dx_km
  net <- adv_flux + dif_flux
  dC <- numeric(N)
  dC[1] <- -net[1] / dx_km
  dC[N] <-  net[N - 1] / dx_km
  dC[2:(N - 1)] <- -(net[2:(N - 1)] - net[1:(N - 2)]) / dx_km
  C + dt * dC
}

drift_obs <- c("1" = -1.39, "2" = -1.06, "3" = -0.88, "4" = -0.87,
               "5" =  1.06, "6" =  4.71, "7" =  1.90, "8" =  1.99,
               "9" =  0.00, "10" = -2.55, "11" = -2.55, "12" = -2.87)

# ---- run 2002-2015, monthly forcing ------------------------------------------
dt <- 1.0                      # hours (drift <= 7.2 km/h => CFL ok)
hours_per_month <- 365 / 12 * 24
C <- numeric(N); C[(N-2):N] <- 1        # March-like start
snaps <- list()

for (year in 2002:2015) {
  for (month in 1:12) {
    env <- get_env(year, month)
    v <- rep(-drift_obs[as.character(month)], N) + env$vcur
    D <- D0 / (1 + env$chl / h)
    for (hh in seq_len(round(hours_per_month / dt))) C <- pde_step(C, v, D, dt)
    snaps[[length(snaps) + 1]] <- data.frame(year = year, month = month,
                                             lat = grid, density = C)
  }
}
mod <- do.call(rbind, snaps)
cat(sprintf("Mass: %.3f (initial 3.000); min density %.5f\n",
            sum(mod$density[mod$year == 2015 & mod$month == 12]), min(mod$density)))

# ---- observed presence by latitude ----------------------------------------------
obs <- aggregate(present ~ year + month + cell_centre_lat,
                 data = cm[cm$cell_centre_lat %in% grid, ], FUN = sum)
obs <- obs[obs$present > 0, ]

cent <- do.call(rbind, lapply(split(obs, list(obs$year, obs$month)), function(d) {
  data.frame(year = d$year[1], month = d$month[1],
             obs_lat = sum(d$cell_centre_lat * d$present) / sum(d$present),
             n = sum(d$present))
}))
modc <- do.call(rbind, lapply(split(mod, list(mod$year, mod$month)), function(d) {
  data.frame(year = d$year[1], month = d$month[1],
             mod_lat = sum(d$lat * d$density) / sum(d$density))
}))
cmp <- merge(cent, modc)
if (nrow(cmp) == 0)
  stop("No overlapping year-months between model output and sightings.")
cmp$err <- cmp$mod_lat - cmp$obs_lat
cat(sprintf("\nCentroid error (model - observed, degrees lat): mean %.1f, |mean| err %.1f, corr %.2f\n",
            mean(cmp$err), mean(abs(cmp$err)), cor(cmp$mod_lat, cmp$obs_lat)))

p1 <- ggplot(cmp, aes(obs_lat, mod_lat, colour = factor(month))) +
  geom_point() + geom_abline(slope = 1, intercept = 0, linetype = 2) +
  coord_equal() +
  scale_colour_viridis_d(name = "month") +
  labs(title = "Hindcast 2002-2015: model vs observed monthly centroid",
       x = "observed centroid latitude", y = "model centroid latitude") +
  theme_minimal()
print(p1)
ggsave("analysis/hindcast_centroid_scatter.png", p1, width = 7, height = 6)

mod$what <- "model"; mod$val <- mod$density / max(mod$density)
obs2 <- aggregate(present ~ year + month + cell_centre_lat, data = cm[cm$cell_centre_lat %in% grid, ], FUN = sum)
names(obs2)[3] <- "lat"; obs2$what <- "observed presence"; obs2$val <- obs2$present / max(obs2$present)
both <- rbind(mod[, c("year", "month", "lat", "what", "val")],
              obs2[, c("year", "month", "lat", "what", "val")])
both$t <- (both$year - 2002) * 12 + both$month
p2 <- ggplot(both, aes(t, lat, fill = val)) +
  geom_raster() + facet_wrap(~what, ncol = 1) +
  scale_fill_viridis_c(name = "scaled") +
  labs(title = "Hindcast vs sightings 2002-2015 (each scaled to its own max)",
       x = "months since Jan 2002", y = "latitude") +
  theme_minimal()
print(p2)
ggsave("analysis/hindcast_hovmoller.png", p2, width = 9, height = 7)

saveRDS(list(model = mod, comparison = cmp), "data/hindcast_2002_2015.rds")
cat("\nSaved data/hindcast_2002_2015.rds + two PNGs in analysis/\n")
