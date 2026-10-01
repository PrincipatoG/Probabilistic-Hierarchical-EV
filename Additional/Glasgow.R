rm(list = ls())
source("R/forecast_function.R")
source("R/reconciliation_function.R")

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(purrr)
library(ggplot2)
library(mgcv)
library(MASS)
library(fastmatrix)
library(data.table)

windows <- generate_rolling_windows()
df_address <- read_csv("Data/dataset_address_main.csv")

df_long <- readRDS("Output_ponctual/Seed_1/reconciled_forecasts_Combination_Combination_Combination.RDS")
H <- as.matrix(readRDS("Data/structural.RDS")[,-1])
stations_glasgow <- which(H[,"Region_Glasgow City"] == 1)
ids_glasgow  <- as.numeric(sub("Station_", "", colnames(H[,33+stations_glasgow])))

res_glasgow <- df_address %>%
  dplyr::filter(Station.ID %in% ids_glasgow) %>%
  dplyr::group_by(Station.ID) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::select(Station.ID, Region, Longitude, Latitude, Address)

df_long_glasgow <- df_long %>% filter(node %in% stations_glasgow)

Y <- df_long_glasgow %>%
  dplyr::select(time, y) %>%
  dplyr::group_by(time) %>%
  dplyr::mutate(id = dplyr::row_number()) %>%
  tidyr::pivot_wider(
    names_from = id,
    values_from = y
  ) %>%
  dplyr::ungroup()
Y <- as.matrix(Y[,-1])
colnames(Y) <- colnames(H[,33+stations_glasgow])

M <- is.na(Y) * 1
n_unique_rows <- nrow(unique(M))

# --- Figures --- 
library(sf)
library(ggplot2)
library(dplyr)
library(readr)
library(purrr)
library(ggspatial)
mycolor <- "#758bd1"
# -------------------------
# 1. LOAD SHAPEFILE
# -------------------------
regions_sf <- st_read("Scrapping/geo_data/pub_commcnc.shp") %>%
    st_transform(27700)

# Dissolve local authorities (clean administrative level)
regions_clean <- regions_sf %>%
  group_by(local_auth) %>%
  summarise(geometry = st_union(geometry), .groups = "drop")

# -------------------------
# Glasgow polygon
# -------------------------
glasgow_sf <- regions_clean %>%
  filter(local_auth == "Glasgow City")

# -------------------------
# 2. LOAD INFRASTRUCTURE (ALL AREAS)
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
infra_glasgow <- infra_all %>% filter(`Local Authority` == "Glasgow City")

# infra_glasgow <- st_as_sf(
#   infra_glasgow,
#   coords = c("Longitude", "Latitude"),
#   crs = 4326
# ) %>%
#   st_transform(27700)

infra_glasgow <- st_as_sf(
  res_glasgow,
  coords = c("Longitude", "Latitude"),
  crs = 4326
) %>%
  st_transform(27700)

glasgow_quarters <- regions_sf %>%
  filter(local_auth == "Glasgow City") %>%
  st_transform(27700)

glasgow_bbox <- st_bbox(glasgow_sf)



# --- Plot ----
options(ggspatial.raster = "cairo")
g2 <- ggplot() +
  
  # -------------------------
# fond OSM
# -------------------------
annotation_map_tile(
  type = "cartolight",
  zoom = 13
) +
  # -------------------------
# quartiers
# -------------------------
geom_sf(
  data = glasgow_quarters,
  fill = NA,
  color = "grey60",
  linewidth = 0.15
) +
  
  # -------------------------
# boundary Glasgow
# -------------------------
geom_sf(
  data = glasgow_sf,
  fill = NA,
  color = "black",
  linewidth = 0.4
) +
  
  # -------------------------
# stations VE (heat effect)
# -------------------------
geom_sf(
  data = infra_glasgow,
  shape = 16,
  color = mycolor,
  size = 7.5,
  alpha = 0.15
) +
  
  geom_sf(
    data = infra_glasgow,
    shape = 16,
    color = mycolor,
    size = 5,
    alpha = 0.4
  ) +
  
  geom_sf(
    data = infra_glasgow,
    shape = 16,
    color = mycolor,
    size = 2,
    alpha = 1
  ) +
  
  coord_sf(expand = FALSE) +
  theme_void()
ggsave(
  filename = "./Figures_AOAS/Glasgow.pdf",
  plot = g2,
  width = 6,
  height = 6
)

#####################################################
################ Simulation #########################
#####################################################

common_nodes <- infra_glasgow

glasgow_centroid <- st_centroid(glasgow_sf)
dx <- 300   # est (+)
dy <- 300   # nord (+)
glasgow_centroid_shifted <- st_sfc(
  st_geometry(glasgow_centroid)[[1]] + c(dx, dy),
  crs = st_crs(glasgow_centroid)
)

glasgow_quarters$dist_to_center <- st_distance(
  st_centroid(glasgow_quarters),
  glasgow_centroid_shifted
)
central_quarters <- glasgow_quarters %>%
  arrange(dist_to_center) %>%
  slice(1:4)
zoom_area <- st_union(central_quarters)

nodes_zoom <- common_nodes %>%
  filter(st_within(., zoom_area, sparse = FALSE))

n_nodes <- nrow(nodes_zoom)

cols_keep <- paste0("Station_", nodes_zoom$Station.ID)

Y_zoom <- Y[, cols_keep, drop = FALSE]
M_zoom <- is.na(Y_zoom) * 1
n_unique_rows <- nrow(unique(M_zoom))

M_A <- M_zoom[1,]
M_B <- M_zoom[125,]

nodes_zoom$state_A <- ifelse(M_A == 0, "active", "inactive")

nodes_zoom$state_B <- ifelse(M_B == 0, "active", "inactive")

nodes_A <- nodes_zoom %>%
  mutate(
    active = state_A == "active"
  )
nodes_B <- nodes_zoom %>%
  mutate(
    active = state_B == "active"
  )
nodes_A_active   <- nodes_A %>% filter(active)
nodes_A_inactive <- nodes_A %>% filter(!active)

nodes_B_active   <- nodes_B %>% filter(active)
nodes_B_inactive <- nodes_B %>% filter(!active)
infra_zoom <- infra_glasgow[
  st_within(infra_glasgow, zoom_area, sparse = FALSE),
]
quarters_zoom <- central_quarters
plot_sim <- function(version = c("A","B")) {
  
  version <- match.arg(version)
  
  if (version == "A") {
    active <- nodes_A_active
    inactive <- nodes_A_inactive
  } else {
    active <- nodes_B_active
    inactive <- nodes_B_inactive
  }
  
  ggplot() +
    
    # fond OSM
    annotation_map_tile(
      type = "cartolight",
      zoom = 15
    ) +
    
    # quartiers
    geom_sf(
      data = central_quarters,
      fill = NA,
      color = "grey20",
      linewidth = 0.4
    ) +
    
    # boundary zoom
    geom_sf(
      data = zoom_area,
      fill = NA,
      color = "grey20",
      linewidth = 0.4
    ) +
    
    geom_sf(
      data = inactive,
      shape = 0,
      color = mycolor,
      size = 3,
      alpha = 0.5
    ) +
    
    geom_sf(
      data = active,
      shape = 15,
      color = mycolor,
      fill = mycolor,
      size = 2.8,
      alpha = 0.9
    ) +
    
    coord_sf(expand = FALSE) +
    theme_void()
}
g_A <- plot_sim("A")
g_B <- plot_sim("B")

# -------------------------
# SAVE
# -------------------------
ggsave("./Figures_AOAS/Glasgow_sim_A.pdf", g_A, width = 5, height = 5)
ggsave("./Figures_AOAS/Glasgow_sim_B.pdf", g_B, width = 5, height = 5)

# --- Mask Conditional Study ---
Y <- df_long_glasgow %>%
  dplyr::select(time, y) %>%
  dplyr::group_by(time) %>%
  dplyr::mutate(id = dplyr::row_number()) %>%
  tidyr::pivot_wider(
    names_from = id,
    values_from = y
  ) %>%
  dplyr::ungroup()

Y <- as.matrix(Y[,-1])
colnames(Y) <- colnames(H[,33+stations_glasgow])

# --- Zoom on the center of Glasgow (if necessary) ---

Y_zoom <- Y[, cols_keep, drop = FALSE]
M_zoom <- is.na(Y_zoom) * 1

mask_size <- rowSums(M_zoom)
ignore    <- which((mask_size > 10))
to_keep <- setdiff(seq(nrow(Y_zoom)), ignore)
Glasgow_regional_node <- df_long %>% filter(node == 16, 
                                            time %in% to_keep)

# Consider 0's as Na's
Y_zoom[is.na(as.matrix(Y_zoom))] <- 0

cols_keep <- paste0("Station_", nodes_zoom$Station.ID)

make_matrix_zoom <- function(df, value_col, stations_glasgow, nodes_zoom, H) {
  
  pred <- df %>%
    dplyr::select(time, {{ value_col }}) %>%
    dplyr::group_by(time) %>%
    dplyr::mutate(id = dplyr::row_number()) %>%
    tidyr::pivot_wider(
      names_from = id,
      values_from = {{ value_col }}
    ) %>%
    dplyr::ungroup()
  
  pred <- as.matrix(pred[, -1])
  
  colnames(pred) <- colnames(H[, 33 + stations_glasgow])
  
  cols_keep <- paste0("Station_", nodes_zoom$Station.ID)
  
  pred_zoom <- pred[, cols_keep, drop = FALSE]
  
  return(pred_zoom)
}

pred_direct   <- make_matrix_zoom(df_long_glasgow, Direct,  stations_glasgow, nodes_zoom, H)

pred_OLS      <- make_matrix_zoom(df_long_glasgow, OLS,     stations_glasgow, nodes_zoom, H)

pred_OLS_ref  <- make_matrix_zoom(df_long_glasgow, OLS_ref, stations_glasgow, nodes_zoom, H)

# Compute the absolute residuals
res_direct  <- abs(Y_zoom - pred_direct)[to_keep,]
res_OLS     <- abs(Y_zoom - pred_OLS)[to_keep,]
res_OLS_ref <- abs(Y_zoom - pred_OLS_ref)[to_keep,]
M_zoom <- M_zoom[to_keep,]
mask_size <- rowSums(M_zoom)

# noms des noeuds
node_names <- colnames(Y_zoom)

build_df <- function(res_mat, method_name) {
  
  map_dfr(seq_len(ncol(res_mat)), function(target_node) {
    
    map_dfr(seq_len(ncol(M_zoom)), function(masked_node) {
      
      idx <- M_zoom[, masked_node] == 1
      
      tibble(
        target_node = node_names[target_node],
        masked_node = node_names[masked_node],
        residual = res_mat[idx, target_node],
        method = method_name
      )
    })
  })
}

df_plot <- bind_rows(
  build_df(res_direct,  "Direct"),
  build_df(res_OLS,     "OLS"),
  build_df(res_OLS_ref, "OLS_ref")
) %>%
  filter(target_node != masked_node)

df_plot <- df_plot %>%
  mutate(
    target_node = str_replace(target_node, "Station_", ""),
    masked_node = str_replace(masked_node, "Station_", "")
  )
df_plot <- df_plot %>%
  mutate(
    target_node = paste0("Node ", target_node),
    masked_node = paste0("Node ", masked_node)
  )
plot_target_node <- function(t, df_plot) {
  
  ggplot(
    dplyr::filter(df_plot, target_node == t),
    aes(
      x = masked_node,
      y = residual,
      fill = method
    )
  ) +
    geom_boxplot(outlier.alpha = 0.15, width = 0.6) +
    theme_classic() +
    theme(
      legend.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "top"
    )
}
targets <- unique(df_plot$target_node)
walk(targets, function(t) {
  
  p <- plot_target_node(t, df_plot)
  
  ggsave(
    filename = file.path(
      "Figures/residuals",
      paste0("mask_effect_target_", t, ".pdf")
    ),
    plot = p,
    width = 7,
    height = 4.5
  )
})

t <- targets[1]
nrow(dplyr::filter(df_plot, target_node == t))

mask_size <- rowSums(M[to_keep,])

df_plot <- Glasgow_regional_node %>%
  mutate(mask_size = mask_size) %>%
  pivot_longer(
    cols = c(Direct, OLS, OLS_ref),
    names_to = "method",
    values_to = "prediction"
  ) %>%
  mutate(
    abs_residual = abs(y - prediction)) %>% 
  filter(mask_size < 60)


df_plot <- df_plot %>%
  mutate(mask_bin = cut(
    mask_size,
    breaks = seq(
      floor(min(mask_size)),
      ceiling(max(mask_size)),
      by = 2
    ),
    include.lowest = TRUE,
    right = FALSE
  ))
ggplot(df_plot,
       aes(x = mask_bin,
           y = abs_residual,
           fill = method)) +
  geom_boxplot(outlier.alpha = 0.15, width = 0.7) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

breaks_main <- seq(18, 46, by = 4)
breaks_all <- c(breaks_main, 62, 64)

df_plot <- Glasgow_regional_node %>%
  mutate(mask_size = mask_size) %>%
  pivot_longer(
    cols = c(Direct, OLS, OLS_ref),
    names_to = "method",
    values_to = "prediction"
  ) %>%
  mutate(
    abs_residual = abs(y - prediction),
    mask_bin = cut(
      mask_size,
      breaks = c(breaks_main, 62, 64.001),
      include.lowest = TRUE,
      right = FALSE,
      ordered_result = TRUE
    )
  )
df_plot <- df_plot %>%
  mutate(mask_bin = as.character(mask_bin)) %>%
  mutate(
    mask_bin = case_when(
      mask_size >= 62 & mask_size < 64 ~ "[62,64]",
      TRUE ~ mask_bin
    )
  )
df_plot <- df_plot %>%
  mutate(
    mask_bin = str_replace(mask_bin, "\\[62,64\\]", "[44,64)")
  )
df_plot <- df_plot %>%
  mutate(
    mask_bin = str_replace(mask_bin, "\\[46,62\\)", "[44,64)")
  )
# df_plot <- df_plot %>%
#   mutate(
#     mask_bin = str_replace(mask_bin, "\\[44,46\\)", "[44,64)")
#   )
g <- ggplot(df_plot,
       aes(x = mask_bin,
           y = abs_residual,
           fill = method)) +
  geom_boxplot(outlier.alpha = 0.15,, width = 0.7) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    legend.title = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
g
ggsave("Figures/residuals/regions.pdf",
       plot = g,
       width = 7,
       height = 4.5)
