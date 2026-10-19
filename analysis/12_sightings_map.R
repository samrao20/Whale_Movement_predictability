# ============================================================================
# 12_sightings_map.R
# Map of humpback sightings off eastern Australia, colour-coded by density,
# faceted by migration season.
#
# Data: data-raw/humpback_all.rds (OBIS/GBIF, presence-only — density here
# means SIGHTING density, which mixes whale density with observer effort;
# whale-watch hotspots like Hervey Bay will glow even in faceted views).
#
# Seasons follow the migration calendar used in the PDE:
#   Northbound  May-Aug,  Breeding  Sep (transit/calving),
#   Southbound  Oct-Dec,  Feeding   Jan-Apr (Southern Ocean)
#
# Run: setwd() to repo root, then source("analysis/12_sightings_map.R")
# Packages: ggplot2 + maps (install.packages("maps") once).
# ============================================================================

library(ggplot2)
library(maps)

if (!file.exists("data-raw/humpback_all.rds"))
  stop("Working directory is not the repo root. Set it via Files pane -> ",
       "More -> Set As Working Directory. Current getwd(): ", getwd())

all_s <- readRDS("data-raw/humpback_all.rds")
lat_col <- intersect(c("decimalLatitude", "latitude", "lat"), names(all_s))[1]
lon_col <- intersect(c("decimalLongitude", "longitude", "lon"), names(all_s))[1]

d <- data.frame(lat = all_s[[lat_col]], lon = all_s[[lon_col]],
                month = all_s$month, year = all_s$year)
d <- d[d$year >= 2002 & !is.na(d$lat) & !is.na(d$lon), ]
d <- d[d$lon >= 140 & d$lon <= 180 & d$lat >= -70 & d$lat <= -10, ]

season <- function(m) {
  ifelse(m %in% 5:8,  "Northbound (May-Aug)",
  ifelse(m == 9,      "Calving/transit (Sep)",
  ifelse(m %in% 10:12,"Southbound (Oct-Dec)",
                      "Feeding (Jan-Apr)")))
}
d$season <- factor(season(d$month),
                   levels = c("Northbound (May-Aug)", "Calving/transit (Sep)",
                              "Southbound (Oct-Dec)", "Feeding (Jan-Apr)"))
cat("Records by season:\n"); print(table(d$season))

# coastline
coast <- map_data("world2", region = c("Australia", "New Zealand", "Antarctica"))
coast$long <- ifelse(coast$long > 180, coast$long - 360, coast$long)
coast <- coast[coast$long >= 140 & coast$long <= 190, ]

# 2-D binned density per season (1-degree bins)
p <- ggplot(d, aes(lon, lat)) +
  geom_bin2d(binwidth = c(1, 1)) +
  geom_path(data = coast, aes(long, lat, group = group),
            colour = "grey30", linewidth = 0.4, inherit.aes = FALSE) +
  scale_fill_viridis_c(name = "sightings\nper 1° cell", trans = "log10") +
  coord_quickmap(xlim = c(140, 180), ylim = c(-70, -10)) +
  facet_wrap(~season, ncol = 2) +
  labs(title = "Humpback sightings by migration season (2002+)",
       subtitle = "Sighting density = whales + observer effort (presence-only data)",
       x = "longitude", y = "latitude") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())
print(p)
ggsave("analysis/sightings_map_seasons.png", p, width = 10, height = 9)
cat("\nSaved analysis/sightings_map_seasons.png\n")
