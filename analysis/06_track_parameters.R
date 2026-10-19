# ============================================================================
# 06_track_parameters.R
# Derive PDE movement parameters (v_max, D0) from the Andrews-Goff et al. 2023
# E1 humpback satellite tracks (Movebank doi:10.5441/001/1.294).
#
# Input : data-raw/humpback_track_steps.csv  (built from the modelled/SSM
#         locations; one row per step with dt_h, dist_km, dlat, month)
# Output: printed parameter estimates + data/track_params.rds
#
# What we estimate:
#   v_max      = 90th percentile of |meridional speed| on directed steps
#                (robust to Argos noise; tags show bursts up to ~8 km/h)
#   v_mig(m)   = median meridional speed of directed steps, by month and
#                phase (northbound May-Aug, southbound Oct-Feb) -> the
#                seasonal drift curve the PDE advection should reproduce
#   D0         = effective along-coast diffusivity on resident (non-directed)
#                steps, D = <dx^2>/(4 dt), expressed in km^2/h
#
# Run: setwd() to the repo root, then source("analysis/06_track_parameters.R")
# ============================================================================

steps <- read.csv("data-raw/humpback_track_steps.csv")
steps$merid_speed <- steps$dlat * 111 / steps$dt_h
steps$speed       <- steps$dist_km / steps$dt_h

directed <- steps[abs(steps$merid_speed) > 0.3, ]
directed$phase <- ifelse(directed$merid_speed > 0, "northbound", "southbound")

# ---- v_max: robust upper quantile of directed meridional speed -------------
v_max <- as.numeric(quantile(abs(directed$merid_speed), 0.90, names = FALSE))
cat(sprintf("v_max (90th pct |meridional speed|, directed steps): %.2f km/h\n", v_max))
cat(sprintf("  (95th pct = %.2f, 99th pct = %.2f km/h)\n",
            quantile(abs(directed$merid_speed), 0.95),
            quantile(abs(directed$merid_speed), 0.99)))

# ---- seasonal drift curve --------------------------------------------------
mig <- aggregate(merid_speed ~ month + phase, data = directed,
                 FUN = function(x) c(n = length(x), median = median(x)))
mig <- do.call(data.frame, mig)
names(mig) <- c("month", "phase", "n", "median_merid_kmh")
mig <- mig[order(mig$month, mig$phase), ]
cat("\nSeasonal drift (median meridional speed by month & phase):\n")
print(mig)

# ---- D0: effective diffusivity from resident steps -------------------------
resident <- steps[abs(steps$merid_speed) <= 0.3 &
                  steps$dt_h >= 12 & steps$dt_h <= 36, ]
D_vals <- resident$dist_km^2 / (4 * resident$dt_h)
D0 <- median(D_vals)
cat(sprintf("\nD0 (median effective diffusivity, resident steps): %.0f km^2/h\n", D0))
cat(sprintf("  (IQR %.0f .. %.0f km^2/h)\n",
            quantile(D_vals, 0.25), quantile(D_vals, 0.75)))

params <- list(v_max = v_max, D0 = D0, seasonal_drift = mig)
saveRDS(params, "data/track_params.rds")
cat("\nSaved data/track_params.rds\n")
