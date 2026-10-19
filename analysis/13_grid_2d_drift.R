# ============================================================================
# 13_grid_2d_drift.R
# Fit the 2-D seasonal drift field from the E1 tag tracks — the advection
# field for the 2-D PDE.
#
# Unlike the 1-D script (06), which used the pre-computed meridional step
# file, here we recompute steps directly from the raw track fixes
# (data-raw/humpback_tracks_clean.csv) so we get BOTH components:
#   v_mer = dlat * 111 / dt           (km/h, + = north)
#   v_zon = dlon * 111 * cos(lat) / dt (km/h, + = east)
# Fixes closer than 1h or further than 48h apart are dropped; speeds > 15 km/h
# are dropped as Argos noise.
#
# Steps are binned on a 4-degree grid x month; bin medians are the drift
# field, smoothed with a 3x3 box mean to fill single-bin holes. Bins with no
# data after smoothing fall back to the corridor-wide monthly meridional
# drift (script 06 values) and zero zonal drift.
#
# Known limitation (from 1-D testing): within a bin, transiting and resident
# whales mix, so medians understate peak directed speeds. Same caveat as the
# 1-D drift curve — this is the mean seasonal flow, not individual velocity.
#
# Outputs: data/drift_field_2d.rds + seasonal quiver maps in analysis/.
# Run: setwd() to repo root, then source("analysis/13_grid_2d_drift.R")
# ============================================================================

library(ggplot2)

if (!file.exists("data-raw/humpback_tracks_clean.csv"))
  stop("Working directory is not the repo root. Set it via Files pane -> ",
       "More -> Set As Working Directory. Current getwd(): ", getwd())

tr <- read.csv("data-raw/humpback_tracks_clean.csv")
tr$timestamp <- as.POSIXct(tr$timestamp, tz = "UTC")
# read.csv imports the True/False 'visible' column as character, so
# tr$visible == TRUE would silently drop every row — compare as string
tr <- tr[as.character(tr$visible) == "TRUE" |
         as.character(tr$visible) == "True", ]
tr <- tr[order(tr$tag, tr$timestamp), ]

# ---- 1. Steps -------------------------------------------------------------------
dt_h  <- as.numeric(diff(tr$timestamp), units = "hours")
same  <- diff(tr$tag) == 0
dlat  <- diff(tr$lat)
dlon  <- diff(tr$lon)
lat0  <- tr$lat[-nrow(tr)]
steps <- data.frame(
  tag   = tr$tag[-nrow(tr)][same],
  month = as.integer(format(tr$timestamp[-nrow(tr)][same], "%m")),
  lat   = lat0[same],
  lon   = tr$lon[-nrow(tr)][same],
  v_mer = (dlat[same] * 111) / dt_h[same],
  v_zon = (dlon[same] * 111 * cos(lat0[same] * pi / 180)) / dt_h[same],
  dt_h  = dt_h[same]
)
steps <- steps[steps$dt_h >= 1 & steps$dt_h <= 48 &
               abs(steps$v_mer) <= 15 & abs(steps$v_zon) <= 15, ]
cat(sprintf("Clean steps: %d from %d fixes\n", nrow(steps), nrow(tr)))

# ---- 2. Bin 4x4 deg x month --------------------------------------------------------
BIN <- 4
domain <- list(lat = c(-68, -12), lon = c(140, 180))
steps$latb <- floor(steps$lat / BIN) * BIN + BIN / 2
steps$lonb <- floor(steps$lon / BIN) * BIN + BIN / 2
steps <- steps[steps$latb >= domain$lat[1] & steps$latb <= domain$lat[2] &
               steps$lonb >= domain$lon[1] & steps$lonb <= domain$lon[2], ]

lats <- seq(domain$lat[1] + BIN/2, domain$lat[2] - BIN/2, by = BIN)
lons <- seq(domain$lon[1] + BIN/2, domain$lon[2] - BIN/2, by = BIN)
nla <- length(lats); nlo <- length(lons)

v_mer <- array(NA_real_, c(nla, nlo, 12))
v_zon <- array(NA_real_, c(nla, nlo, 12))
nbin  <- array(0L,       c(nla, nlo, 12))
for (i in seq_len(nla)) for (j in seq_len(nlo)) for (m in 1:12) {
  s <- steps[steps$latb == lats[i] & steps$lonb == lons[j] & steps$month == m, ]
  nbin[i, j, m] <- nrow(s)
  if (nrow(s) >= 10) {                       # min steps for a stable median
    v_mer[i, j, m] <- median(s$v_mer)
    v_zon[i, j, m] <- median(s$v_zon)
  }
}
cat(sprintf("Bins with >=10 steps: %d of %d\n", sum(nbin >= 10), nla * nlo * 12))

# ---- 3. Smooth (3x3 box, preserves NA) then fallback ----------------------------------
smooth3 <- function(A) {
  S <- array(NA_real_, dim(A))
  for (m in 1:12) for (i in seq_len(nla)) for (j in seq_len(nlo)) {
    ii <- max(1, i - 1):min(nla, i + 1)
    jj <- max(1, j - 1):min(nlo, j + 1)
    vals <- A[ii, jj, m]
    if (any(!is.na(vals))) S[i, j, m] <- mean(vals, na.rm = TRUE)
  }
  S
}
v_mer_s <- smooth3(v_mer)
v_zon_s <- smooth3(v_zon)

drift_corridor <- c("1" = -1.39, "2" = -1.06, "3" = -0.88, "4" = -0.87,
                    "5" =  1.06, "6" =  4.71, "7" =  1.90, "8" =  1.99,
                    "9" =  0.00, "10" = -2.55, "11" = -2.55, "12" = -2.87)
for (m in 1:12) {
  v_mer_s[, , m][is.na(v_mer_s[, , m])] <- drift_corridor[as.character(m)]
  v_zon_s[, , m][is.na(v_zon_s[, , m])] <- 0
}

saveRDS(list(v_mer = v_mer_s, v_zon = v_zon_s, lats = lats, lons = lons,
             bin_deg = BIN, n_steps = nbin),
        "data/drift_field_2d.rds")

# ---- 4. Seasonal quiver maps ---------------------------------------------------------
season_months <- list("Northbound (Jun)" = 6, "Breeding (Jul)" = 7,
                      "Southbound (Nov)" = 11, "Feeding (Jan)" = 1)
plt <- do.call(rbind, lapply(names(season_months), function(nm) {
  m <- season_months[[nm]]
  expand.grid(lat = lats, lon = lons, season = nm, month = m,
              stringsAsFactors = FALSE)
}))
plt$v_mer <- mapply(function(la, lo, m)
  v_mer_s[which(lats == la), which(lons == lo), m], plt$lat, plt$lon, plt$month)
plt$v_zon <- mapply(function(la, lo, m)
  v_zon_s[which(lats == la), which(lons == lo), m], plt$lat, plt$lon, plt$month)
# flag data-supported bins (direct or 3x3-smoothed neighbours) vs pure fallback
plt$n_steps <- mapply(function(la, lo, m)
  nbin[which(lats == la), which(lons == lo), m], plt$lat, plt$lon, plt$month)
plt$src <- ifelse(plt$n_steps >= 10, "tag data (n>=10)", "fallback (corridor drift)")
plt$season <- factor(plt$season, levels = names(season_months))
sc <- 0.6   # arrow scale: 1 km/h ~ sc degrees of lat
plt$lat_e <- plt$lat + plt$v_mer * sc
plt$lon_e <- plt$lon + plt$v_zon * sc

p <- ggplot(plt, aes(lon, lat, colour = src)) +
  geom_segment(aes(xend = lon_e, yend = lat_e),
               arrow = arrow(length = unit(0.12, "cm"))) +
  scale_colour_manual(values = c("tag data (n>=10)" = "firebrick",
                                 "fallback (corridor drift)" = "grey70"),
                      name = NULL) +
  facet_wrap(~season, ncol = 2) +
  coord_quickmap(xlim = domain$lon, ylim = domain$lat) +
  labs(title = "Tag-derived seasonal drift field (km/h, 4° bins, smoothed)",
       subtitle = "Red = fitted from tag steps; grey = corridor-wide fallback (95% of bin-months have no tags)",
       x = "longitude", y = "latitude") +
  theme_minimal()
print(p)
ggsave("analysis/drift_field_2d_quivers.png", p, width = 10, height = 9)
cat("\nSaved data/drift_field_2d.rds + analysis/drift_field_2d_quivers.png\n")
