rm(list = ls())

# -------------------------
# LIBRARIES
# -------------------------
library(sf)
library(ggplot2)
library(dplyr)
library(readr)
library(purrr)
library(ggspatial)

mycolor <- "#758bd1"

# -------------------------
# 1. LOAD + CLEAN SHAPEFILE
# -------------------------
regions_sf <- st_read("Scrapping/geo_data/pub_commcnc.shp") %>%
  st_transform(27700) %>%
  st_make_valid()

# Dissolve into local authorities (structure only, not for detailed drawing)
regions_clean <- regions_sf %>%
  group_by(local_auth) %>%
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  st_make_valid()

# -------------------------
# SCOTLAND CLEAN SURFACE (KEY STEP)
# -------------------------
scotland_poly <- st_union(regions_clean) %>%
  st_make_valid()

# Clean coastline (removes micro artifacts)
scotland_poly <- st_buffer(scotland_poly, 0)

# Optional but very effective smoothing
scotland_poly <- st_simplify(
  scotland_poly,
  dTolerance = 300,
  preserveTopology = TRUE
)

# SINGLE CLEAN OUTLINE (NO RESIDUAL INTERNAL LINES POSSIBLE)
scotland_outline <- st_boundary(scotland_poly)

# -------------------------
# 2. INFRASTRUCTURE DATA
# -------------------------
base_path <- "Scrapping/data/result/"
local_auths <- list.dirs(base_path, recursive = FALSE, full.names = FALSE)

read_infra <- function(la) {
  path <- file.path(base_path, la, "charging_infrastructure.csv")
  
  if (!file.exists(path)) return(NULL)
  
  read_csv(
    path,
    col_types = cols(
      `CP ID` = col_character(),
      Latitude = col_double(),
      Longitude = col_double()
    ),
    show_col_types = FALSE
  ) %>%
    mutate(local_auth = la)
}

infra_all <- map_dfr(local_auths, read_infra)

infra_sf <- st_as_sf(
  infra_all,
  coords = c("Longitude", "Latitude"),
  crs = 4326
) %>%
  st_transform(27700)

# -------------------------
# 3. PLOT
# -------------------------
g_scotland <- ggplot() +
  
  # -------------------------
# CARTO BASE
# -------------------------
annotation_map_tile(
  type = "cartolight",
  zoom = 8
) +
  
  # -------------------------
# LIGHT STRUCTURE (local authorities ONLY as fill, no borders)
# -------------------------
geom_sf(
  data = regions_clean,
  fill = NA,
  color = "grey85",
  linewidth = 0
) +
  
  # -------------------------
# CLEAN SCOTLAND OUTLINE (ONLY BORDER THAT MATTERS)
# -------------------------
geom_sf(
  data = scotland_outline,
  fill = NA,
  color = "black",
  linewidth = 0.6,
  lineend = "round"
) +
  
  # -------------------------
# INFRASTRUCTURE (density effect)
# -------------------------
geom_sf(
  data = infra_sf,
  shape = 16,
  color = mycolor,
  size = 2.5,
  alpha = 0.05
) +
  geom_sf(
    data = infra_sf,
    shape = 16,
    color = mycolor,
    size = 1.2,
    alpha = 0.15
  ) +
  geom_sf(
    data = infra_sf,
    shape = 16,
    color = mycolor,
    size = 0.7,
    alpha = 1
  ) +
  
  coord_sf(expand = FALSE) +
  theme_void()

# -------------------------
# 4. SAVE
# -------------------------
ggsave(
  filename = "./Figures_AOAS/Scotland.pdf",
  plot = g_scotland,
  width = 6,
  height = 6
)
