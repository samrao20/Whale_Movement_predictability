# 02_fit_suitability.R
# Fit seasonal suitability models (presence-background GAMs) and extract the
# thermal envelope parameters (T_opt, sigma_T) needed for the PDE's S(T, z).
#
# Inputs : data/sightings_with_env.rds, data/background_with_env.rds
# Outputs: data/suitability_models.rds, data/thermal_envelopes.csv
#
# Run after: analysis/01_join_sightings_water.R

library(dplyr)
library(mgcv)

sightings  <- readRDS("data/sightings_with_env.rds")
background <- readRDS("data/background_with_env.rds")

# --- Build presence/background dataset ---------------------------------------
pres <- sightings |>
  filter(!is.na(sst)) |>
  transmute(pres = 1, year, month, lat = latitude, lon = longitude,
            sst, chl)

bg <- background |>
  filter(!is.na(sst)) |>
  transmute(pres = 0, year, month, lat = cell_lat, lon = cell_lon,
            sst, chl)

all_dat <- bind_rows(pres, bg) |>
  mutate(season = case_when(
    month %in% 6:9            ~ "breeding",   # Jun-Sep: QLD calving grounds
    month %in% c(11, 12, 1, 2) ~ "feeding",   # Nov-Feb: southern feeding grounds
    TRUE                      ~ "transit"
  ))

# Presence-background GAMs are sensitive to the background ratio; use case-
# control weights (presences weight 1, background scaled to sum ~ n presences)
fit_season <- function(dat, label) {
  d <- dat |> filter(season == label)

  pres_d <- d |> filter(pres == 1, !is.na(sst))
  n_pres <- nrow(pres_d)

  # 1. Cap background to the presence SST range (+/- 1 C) so the model can't
  #    extrapolate suitability into environments no whale ever experienced.
  sst_lo <- min(pres_d$sst) - 1
  sst_hi <- max(pres_d$sst) + 1
  # 2. Restrict background to the latitudinal band where presences occur
  #    (+/- 2 deg) -- background must be *reachable* habitat.
  lat_lo <- min(pres_d$lat) - 2
  lat_hi <- max(pres_d$lat) + 2
  bg_d <- d |>
    filter(pres == 0, !is.na(sst),
           sst >= sst_lo, sst <= sst_hi,
           lat >= lat_lo, lat <= lat_hi)

  d <- bind_rows(pres_d, bg_d)

  # 3. Chl only where coverage is adequate (breeding season is ~47% missing)
  use_chl <- label != "breeding" &&
             mean(is.na(d$chl)) < 0.2 &&
             length(unique(na.omit(d$chl))) >= 5

  smooth_terms <- c("s(sst, k = 5)")
  if (use_chl) smooth_terms <- c(smooth_terms, "s(chl, k = 4)")
  smooth_terms <- c(smooth_terms, "s(lat, k = 6)")
  fml <- as.formula(paste("pres ~", paste(smooth_terms, collapse = " + ")))

  d_fit <- if (use_chl) filter(d, !is.na(chl)) else d

  cat(label, "| presences:", n_pres,
      "| background kept:", nrow(bg_d),
      "| SST window: [", round(sst_lo, 1), ",", round(sst_hi, 1), "]",
      "| formula:", deparse(fml), "\n")

  d_fit <- d_fit |>
    mutate(w = if_else(pres == 1, 1, n_pres / sum(pres == 0)))

  m <- gam(fml, family = binomial, data = d_fit, weights = w, method = "REML")

  # --- Extract thermal envelope: T_opt and sigma_T from the SST curve -------
  # Evaluate only within the observed presence SST range (no extrapolation).
  grid <- data.frame(
    sst = seq(min(pres_d$sst), max(pres_d$sst), length.out = 400),
    lat = median(d_fit$lat)
  )
  if (use_chl) grid$chl <- median(d_fit$chl, na.rm = TRUE)
  fit <- predict(m, newdata = grid, type = "response")
  t_opt <- grid$sst[which.max(fit)]

  # sigma_T: width of the curve at half-max (FWHM/2.355, Gaussian analogy)
  half <- max(fit) / 2
  above <- grid$sst[fit >= half]
  fwhm <- max(above) - min(above)
  sigma_T <- fwhm / 2.355

  list(model = m,
       envelope = data.frame(season = label,
                             n_presences = n_pres,
                             T_opt = round(t_opt, 2),
                             sigma_T = round(sigma_T, 2)),
       curve = data.frame(season = label, sst = grid$sst, suitability = fit))
}

breeding <- fit_season(all_dat, "breeding")
feeding  <- fit_season(all_dat, "feeding")

# --- Save --------------------------------------------------------------------
envelopes <- bind_rows(breeding$envelope, feeding$envelope)
curves    <- bind_rows(breeding$curve,    feeding$curve)

saveRDS(list(breeding = breeding$model, feeding = feeding$model),
        "data/suitability_models.rds")
write.csv(envelopes, "data/thermal_envelopes.csv", row.names = FALSE)
write.csv(curves,    "data/suitability_curves.csv", row.names = FALSE)

print(envelopes)

# --- Plot the two thermal envelopes ------------------------------------------
library(ggplot2)
p <- ggplot(curves, aes(sst, suitability, colour = season)) +
  geom_line(linewidth = 1.2) +
  geom_vline(data = envelopes, aes(xintercept = T_opt, colour = season),
             linetype = "dashed") +
  labs(x = "Sea surface temperature (°C)", y = "Relative suitability",
       title = "Humpback thermal envelopes, East Australia",
       subtitle = "Dashed lines = T_opt per season") +
  theme_minimal()
ggsave("analysis/thermal_envelopes.png", p, width = 7, height = 4.5)
print(p)
