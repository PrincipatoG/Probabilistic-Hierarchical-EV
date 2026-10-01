rm(list = ls())

source("R/forecast_function.R")
source("R/reconciliation_function.R")
source("R/metric_function.R")

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(purrr)
library(ggplot2)
library(mgcv)
library(MASS)
library(fastmatrix)
library(scales)

# =====================================================
# PARAMETERS
# =====================================================

results_root <- "Output_ponctual"
figures_dir <- "Figures"

dir.create(figures_dir, showWarnings = FALSE)

# =====================================================
# DETECT ALL SEED FOLDERS
# =====================================================

seed_dirs <- list.dirs(
  results_root,
  recursive = FALSE,
  full.names = TRUE
)

seed_dirs <- seed_dirs[
  grepl("Seed_", basename(seed_dirs))
]

cat("Detected seed folders:\n")
print(basename(seed_dirs))

test <- readRDS(file = "Output_ponctual/Seed_1/results_stations.RDS")
forecast_cols <- c(
  "LOCAL_GAM",
  "LOCAL_RF",
  "LOCAL_XBG",
  "LOCAL_ARIMA",
  "GLOBAL_RF",
  "GLOBAL_XGB",
  "GLOBAL_GAM",
  "tabICL",
  "Combi"
)

# =====================================================
# LOAD ALL RESULTS ACROSS SEEDS
# =====================================================

prev_nat_list <- list()
prev_reg_list <- list()
prev_sta_list <- list()
prev_recon_list <- list()

for (seed_dir in seed_dirs) {
  
  seed_id <- str_extract(basename(seed_dir), "\\d+")
  
  result_files <- list.files(
    seed_dir,
    pattern = "\\.RDS$",
    full.names = TRUE
  )
  
  result_files <- result_files[
    !grepl("global", basename(result_files), ignore.case = TRUE)
  ]
  result_files <- result_files[
    !grepl("reconciled", basename(result_files), ignore.case = TRUE)
  ]

  cat("\nSeed", seed_id, ":\n")
  print(basename(result_files))
  
  for (f in result_files) {

    dat <- readRDS(f)
    
    # Columns to keep
    cols_keep <- c(
      "window_id",
      "Date",
      "type",
      forecast_cols,
      "Consumed_kWh",
      "y", 
      "node", 
      "time", 
      "Direct",
      "OLS",
      "OLS_refined",
      "node_type", 
      "Region",
      "Station.ID"
    )
    # Keep only columns that are present
    cols_keep <- intersect(cols_keep, names(dat))
    
    dat_reduced <- dat %>%
      dplyr::select(all_of(cols_keep)) %>%
      mutate(seed = as.integer(seed_id))
    
    file_name <- basename(f)
    
    if (grepl("scotland", file_name, ignore.case = TRUE)) {
      
      prev_nat_list[[length(prev_nat_list) + 1]] <- dat_reduced
      
    } else if (grepl("region", file_name, ignore.case = TRUE)) {
      
      prev_reg_list[[length(prev_reg_list) + 1]] <- dat_reduced
      
    } else if (grepl("station", file_name, ignore.case = TRUE)) {
      
      prev_sta_list[[length(prev_sta_list) + 1]] <- dat_reduced
      
    }
}
}

for (seed_dir in seed_dirs) {
  seed_id <- str_extract(basename(seed_dir), "\\d+")
  
  result_file <- paste0(seed_dir, "/reconciled_forecasts_Combination_Combination_Combination.RDS")
  
  dat <- readRDS(result_file)
  
  # Columns to keep
  cols_keep <- c(
    "window_id",
    "Date",
    "type",
    forecast_cols,
    "Consumed_kWh",
    "y", 
    "node", 
    "time", 
    "Direct",
    "OLS",
    "OLS_refined",
    "node_type", 
    "Region",
    "Station.ID"
  )
  # Keep only columns that are present
  cols_keep <- intersect(cols_keep, names(dat))
  
  dat_reduced <- dat %>%
    dplyr::select(all_of(cols_keep)) %>%
    mutate(seed = as.integer(seed_id))
  
  file_name <- basename(f)
  prev_recon_list[[length(prev_recon_list) + 1]] <- dat_reduced
}  


prev_nat <- bind_rows(prev_nat_list) %>% filter(type == "test")
prev_reg <- bind_rows(prev_reg_list) %>% filter(type == "test")
prev_sta <- bind_rows(prev_sta_list) %>% filter(type == "test")
prev_recon <- bind_rows(prev_recon_list)

time_index <- prev_recon %>%
  distinct(time) %>%
  mutate(Date = sort(unique(prev_nat$Date))) %>%
  dplyr::select(Date, time)

prev_nat2 <- prev_nat %>%
  left_join(time_index, by = "Date") %>%
  mutate(
    node = 1,
    node_type = "national"
  ) %>%
  rename(y = Consumed_kWh, 
         XGB = LOCAL_ARIMA) %>%        # Adapt if necessary
  dplyr::select(-Date) %>%
  rename_with(~str_remove(.x, "^(LOCAL_|GLOBAL_)")) %>%
  dplyr::select(time, seed, node, GAM, XGB, tabICL, RF)

region_key <- tibble(
  Region = sort(unique(prev_reg$Region)),
  node = 2:(length(unique(prev_reg$Region)) + 1)
)

prev_reg2 <- prev_reg %>%
  left_join(region_key, by = "Region") %>%
  left_join(time_index, by = "Date") %>%
  mutate(node_type = "regional") %>%
  rename(y = Consumed_kWh) %>%        # Adapt if necessary
  dplyr::select(-Date) %>%
  rename_with(~str_remove(.x, "^(LOCAL_|GLOBAL_)")) %>%
  dplyr::select(time, seed, node, GAM, XGB, tabICL, RF)

prev_sta2 <- prev_sta %>%
  left_join(time_index, by = "Date") %>%
  mutate(
    node = Station.ID + 33,
    node_type = "station"
  ) %>%
  rename(y = Consumed_kWh) %>%        # Adapt if necessary
  dplyr::select(-Date) %>%
  rename_with(~str_remove(.x, "^(LOCAL_|GLOBAL_)")) %>%
  dplyr::select(time, seed, node, GAM, XGB, tabICL, RF)

prev_vect <- bind_rows(prev_nat2, prev_reg2, prev_sta2)
prev_recon2 <-
  prev_recon %>%
  left_join(prev_vect,
            by = c("time", "seed", "node")) %>%
  rename(XGB_ARIMA = XGB)

valid_nodes <- prev_recon2 %>%
  group_by(seed, node) %>%
  summarise(n = sum(!is.na(y)), .groups = "drop") %>%
  filter(n >= 20)

dat <- prev_recon2 %>%
  semi_join(valid_nodes, by = c("seed", "node")) %>%
  rename(rOLS = OLS_refined)

# methods <- c("Direct", "OLS", "rOLS",
#              "GAM", "XGB_ARIMA", "tabICL", "RF")

saveRDS(dat, "Output_ponctual/aggregated_forecast.RDS")

dat <- readRDS("Output_ponctual/aggregated_forecast.RDS")
methods <- c("Direct", "GAM", "XGB", "tabICL", "RF", "Combi")#, "OLS", "rOLS")

cols <- c(
  RF          = "grey75",
  XGB         = "grey75",
  ARIMA       = "grey75",
  GAM         = "grey75",
  tabICL      = "grey75",
  Combi = "grey30"
)
valid_nodes <- dat %>%
  group_by(seed, node) %>%
  summarise(n = sum(!is.na(y)), .groups = "drop") %>%
  filter(n >= 20)
dat <- dat %>%
  semi_join(valid_nodes, by = c("seed", "node")) 

dat_long <- dat %>%
  pivot_longer(
    cols = c(Direct, GAM, XGB_ARIMA, tabICL, RF),#, OLS, rOLS),
    names_to = "method",
    values_to = "prediction"
  ) %>%
  mutate(
    method = case_when(
      method == "XGB_ARIMA" & node_type == "national" ~ "ARIMA",
      method == "XGB_ARIMA"                           ~ "XGB",
      method == "Direct" ~ "Combi",
      TRUE                                            ~ method
    )
  )
dat_metrics <- dat_long %>%
  filter(!is.na(y), !is.na(prediction)) %>%
  group_by(node_type, node, method, seed) %>%
  filter() %>%
  summarise(
    MAE   = mean(abs(y - prediction)),
    MSE   = mean((y - prediction)^2),
    RMSE  = sqrt(MSE),
    NMAE  = mean(abs(y - prediction)) / mean(abs(y)),
    NRMSE = RMSE / mean(abs(y)),
    MAPE  = mean(abs(y - prediction) / abs(y)),
    n_obs = sum(!is.na(y)), 
    .groups = "drop"
  )
aggregate_metric <- function(data, metric_name) {
  
  data %>%
    group_by(node_type, node, method, seed) %>%
    summarise(
      # value = mean(.data[[metric_name]]),
      value = weighted.mean(.data[[metric_name]], w = n_obs, na.rm = TRUE),
      .groups = "drop"
    )
}


base_theme <- theme_classic(base_size = 18) +
  theme(
    legend.position = "none",
    
    axis.title = element_blank(),
    
    axis.text.x = element_text(
      size = 18,
      angle = 35,
      hjust = 0.55,
      vjust = 0.8
    ),
    axis.text.y = element_text(
      size = 20,
      angle = 90,
      hjust = 0.5),
    
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.4),
    
    plot.title = element_text(
      face = "bold",
      size = 20,
      hjust = 0.5
    )
  )

make_plot_metric <- function(df, title, x_levels, ylab) {
  
  ggplot(df, aes(method, value, fill = method)) +
    geom_violin(
      width = 0.9,
      trim = TRUE,
      colour = NA,
      alpha = 0.6
    ) +
    geom_boxplot(
      aes(colour = method),
      width = 0.18,
      outlier.shape = 16,
      outlier.size = 0.6,
      outlier.alpha = 0.25,
      outlier.color = "grey75",
      linewidth = 0.3
    ) +
    scale_colour_manual(
      values = c(
        RF = "grey10",
        XGB = "grey10",
        GAM = "grey10",
        tabICL = "grey10",
        Combi = "grey80"
      ) 
      ) + 
    scale_fill_manual(values = cols) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(
      n.breaks = 4
    ) +
    labs(
      title = title,
      y = ylab,
      x = NULL
    ) +
    base_theme
}

plot_metric <- function(metric_name, ylab) {
  
  df_metric <- aggregate_metric(dat_metrics, metric_name) %>%
    filter(method %in% methods)
  
  make_panel <- function(type, x_levels) {
    
    p <- df_metric %>%
      filter(node_type == type) %>%
      mutate(method = factor(method, levels = x_levels)) %>%
      make_plot_metric(
        title = paste0(toupper(substr(type, 1, 1)), substr(type, 2, nchar(type))),
        x_levels = x_levels,
        ylab = ylab
      )
    
    if (type == "station") {
      p <- p +
        coord_cartesian(
          ylim = c(
            0,
            quantile(
              df_metric %>% 
                filter(node_type == "station") %>% 
                pull(value),
              0.99,
              na.rm = TRUE
            )
          )
        )
    }
    
    p
  }
  
  list(
    national = make_panel("national",
                          c("RF", "GAM", "tabICL", "Combi")),
    regional = make_panel("regional",
                          c("RF", "XGB", "GAM", "tabICL", "Combi")),
    station  = make_panel("station",
                          c("RF", "XGB", "GAM", "tabICL", "Combi"))
  )
}

## ------------------------------------------------------------------
## Generate figures
## ------------------------------------------------------------------

figs <- plot_metric("RMSE", "Mean RMSE")

figs$national
figs$regional
figs$station

ggsave(
  "Figures/rmse_national.pdf",
  figs$national,
  width = 3,
  height = 6.5,
  device = cairo_pdf
)

ggsave(
  "Figures/rmse_regional.pdf",
  figs$regional,
  width = 3,
  height = 6.5,
  device = cairo_pdf
)

ggsave(
  "Figures/rmse_station.pdf",
  figs$station,
  width = 3,
  height = 6.5,
  device = cairo_pdf
)


cols <- c(
  Combi = "grey30",
  OLS         = "grey80",
  rOLS        = "grey15"
)
methods <- c("Direct", "Combi", "OLS", "rOLS")

valid_nodes <- dat %>%
  group_by(seed, node) %>%
  summarise(n = sum(!is.na(y)), .groups = "drop") %>%
  filter(n >= 20)
dat <- dat %>%
  semi_join(valid_nodes, by = c("seed", "node")) 

dat_long <- dat %>%
  pivot_longer(
    cols = c(Direct, OLS, rOLS),
    names_to = "method",
    values_to = "prediction"
  ) %>%
  mutate(
    method = case_when(
      method == "XGB_ARIMA" & node_type == "national" ~ "ARIMA",
      method == "XGB_ARIMA"                           ~ "XGB",
      method == "Direct" ~ "Combi",
      TRUE                                            ~ method
    )
  )
dat_metrics <- dat_long %>%
  filter(!is.na(y), !is.na(prediction)) %>%
  group_by(node_type, node, method, seed) %>%
  filter() %>%
  summarise(
    MAE   = mean(abs(y - prediction)),
    MSE   = mean((y - prediction)^2),
    RMSE  = sqrt(MSE),
    NMAE  = mean(abs(y - prediction)) / mean(abs(y)),
    NRMSE = RMSE / mean(abs(y)),
    MAPE  = mean(abs(y - prediction) / abs(y)),
    Mean_y = mean(abs(y)),
    n_obs = n(),
    .groups = "drop"
  )
dat_metrics2 <- dat_metrics %>%
  group_by(method, seed) %>%
  summarise(NnRMSE = sqrt(mean(MSE)) / mean(Mean_y) ,
  nRMSE = sqrt(mean(MSE)),
  .groups = "drop")
aggregate_metric <- function(data, metric_name) {
  
  data %>%
    group_by(node_type, node, method, seed) %>%
    summarise(
      value = mean(.data[[metric_name]]),
      .groups = "drop"
    )
}

make_plot_metric <- function(df, title, x_levels, ylab) {
  
  ggplot(df, aes(method, value, fill = method)) +
    geom_violin(
      width = 0.9,
      trim = TRUE,
      colour = NA,
      alpha = 0.6
    ) +
    geom_boxplot(
      aes(colour = method),
      width = 0.18,
      outlier.shape = 16,
      outlier.size = 0.6,
      outlier.alpha = 0.25,
      outlier.color = "grey40",
      linewidth = 0.3
    ) +
    scale_fill_manual(values = cols) +
    scale_colour_manual(
      values = c(
        OLS = "grey20",
        Combi = "grey80",
        rOLS = "grey90"
      )
    ) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(
      n.breaks = 3.25
    ) +
    labs(
      title = title,
      y = ylab,
      x = NULL
    ) +
    base_theme
}

plot_metric <- function(metric_name, ylab) {
  
  df_metric <- aggregate_metric(dat_metrics, metric_name) %>%
    filter(method %in% methods)
  
  df_metric2 <- dat_metrics2
  
  make_panel <- function(type, x_levels) {
    
    if (type == "global") {
      print("ici")
      p <- df_metric2 %>%
        mutate(
          value = nRMSE,
          method = factor(method, levels = x_levels)
        ) %>%
        make_plot_metric(
          title = "Global",
          x_levels = x_levels,
          ylab = ylab
        )
      
      p
      
    } else {
      
      p <- df_metric %>%
        filter(node_type == type) %>%
        mutate(method = factor(method, levels = x_levels)) %>%
        make_plot_metric(
          title = paste0(toupper(substr(type, 1, 1)), substr(type, 2, nchar(type))),
          x_levels = x_levels,
          ylab = ylab
        )
      
      if (type == "station") {
        ymax <- df_metric %>%
          filter(node_type == "station") %>%
          group_by(method) %>%
          summarise(
            q95 = quantile(value, 0.99, na.rm = TRUE),
            .groups = "drop"
          ) %>%
          pull(q95) %>%
          max()
        
        p <- p +
          coord_cartesian(
            ylim = c(0, ymax)
          )
      }
      
      p
      
    }
  }  
  list(
    global   = make_panel("global",   c("Combi", "OLS", "rOLS")),
    national = make_panel("national", c("Combi", "OLS", "rOLS")),
    regional = make_panel("regional", c("Combi", "OLS", "rOLS")),
    station  = make_panel("station",  c("Combi", "OLS", "rOLS"))
  )
}

## ------------------------------------------------------------------
## Generate figures
## ------------------------------------------------------------------

figs <- plot_metric("RMSE", "Mean RMSE")

## Display
figs$global
figs$national
figs$regional
figs$station

## Save
ggsave("Figures/recon_rmse_global.pdf",
       figs$global,
       width = 2.25, height = 6.5, device = cairo_pdf)

ggsave("Figures/recon_rmse_national.pdf",
       figs$national,
       width = 2.25, height = 6.5, device = cairo_pdf)

ggsave("Figures/recon_rmse_regional.pdf",
       figs$regional,
       width = 2.25, height = 6.5, device = cairo_pdf)

ggsave("Figures/recon_rmse_station.pdf",
       figs$station,
       width = 2.25, height = 6.5, device = cairo_pdf)
