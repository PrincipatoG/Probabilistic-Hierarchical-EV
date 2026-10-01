#!/usr/bin/env Rscript

rm(list = ls())

source("R/forecast_function.R")
source("R/reconciliation_function.R")
source("R/metric_function.R")

library(argparser)

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(purrr)
library(ggplot2)
library(mgcv)
library(MASS)
library(fastmatrix)
library(Matrix)
library(RSpectra)

# -----------------------------
# ARGUMENT PARSER
# -----------------------------

p <- arg_parser("Hierarchical forecast reconciliation")

p <- add_argument(p, 
                  "--seed", 
                  help="Seed for generating windows", 
                  default=1, 
                  type="integer")

p <- add_argument(
  p,
  "--alpha",
  help = "Alpha level",
  type = "numeric",
  default = 0.1
)

p <- add_argument(
  p,
  "--national_model",
  help = "Model used for national forecasts",
  default = "Combination"
)

p <- add_argument(
  p,
  "--regional_model",
  help = "Model used for regional forecasts",
  default = "Combination"
)

p <- add_argument(
  p,
  "--station_model",
  help = "Model used for station forecasts",
  default = "Combination"
)

p <- add_argument(
  p,
  "--run_gather",
  help = "Run the forecast gathering step",
  type = "logical",
  default = FALSE
)

p <- add_argument(
  p,
  "--results_national",
  help = "Path to national forecasts",
  default = "Output_ponctual/Seed_1/results_scotland.RDS"
)

p <- add_argument(
  p,
  "--results_regional",
  help = "Path to regional forecasts",
  default = "Output_ponctual/Seed_1/results_regions.RDS"
)

p <- add_argument(
  p,
  "--results_station",
  help = "Path to station forecasts",
  default = "Output_ponctual/Seed_1/results_stations.RDS"
)

p <- add_argument(
  p,
  "--structural_matrix",
  help = "Path to structural hierarchy matrix",
  default = "Data/structural.RDS"
)

p <- add_argument(
  p,
  "--output",
  help = "Output path for reconciled forecasts",
  default = NA
)

argv <- parse_args(p)

# -----------------------------
# INITIALISATION
# -----------------------------
seed_windows  <- argv$seed

alpha <- argv$alpha

model_names <- list(
  argv$national_model,
  argv$regional_model,
  argv$station_model
)

run_gather <- argv$run_gather

results_national_path <- paste0("Output_ponctual/Seed_", seed_windows, "/results_scotland.RDS")
results_regional_path <- paste0("Output_ponctual/Seed_", seed_windows, "/results_regions.RDS")
results_station_path  <- paste0("Output_ponctual/Seed_", seed_windows, "/results_stations.RDS")

structural_matrix_path <- argv$structural_matrix

model_tag <- paste(
  argv$national_model,
  argv$regional_model,
  argv$station_model,
  sep = "_"
)

output_path <- ifelse(
  is.na(argv$output),
  paste0(
    "Output_ponctual/Seed_", seed_windows,"/reconciled_forecasts_",
    model_tag,
    ".RDS"
  ),
  argv$output
)

windows <- generate_rolling_windows(seed_windows)

# -----------------------------
# OPTIONAL: GATHER FORECASTS
# -----------------------------
if (run_gather) {
  
  # --- Gather all station forecasts ---
  res_sta_global <- lapply(seq_along(windows), function(i) {
    
    readRDS(
      paste0(
        "Output_ponctual/Seed_", seed_windows,"/results_stations_type_global_lags_TRUE_1398/global_win_",
        i,
        ".RDS"
      )
    )
  })
  
  res_sta <- map2_dfr(
    res_sta_global,
    seq_along(res_sta_global),
    ~ mutate(.x, window_id = .y)
  )
  
  station_foundation_file <- paste0("Output_ponctual/Seed_", seed_windows,"/results_stations_tabICL.csv")
  
  if (file.exists(station_foundation_file)) {
    
    res_sta_foundation <- read.csv(
      station_foundation_file,
      sep = ";"
    )
    
    res_sta$tabICL <- res_sta_foundation$pred_tabICL
    
  } else {
    
    warning(
      "Foundation file not found: ",
      station_foundation_file
    )
  }
  
  forecast_cols <- c(
    "LOCAL_GAM",
    "GLOBAL_RF",
    "GLOBAL_XGB",
    "GLOBAL_GAM",
    "tabICL"
  )
  
  available_cols <- intersect(
    forecast_cols,
    names(res_sta)   # or res_reg / res_nat
  )
  
  res_sta$Combination <- rowMeans(
    res_sta[, available_cols, drop = FALSE],
    na.rm = TRUE
  )
  saveRDS(res_sta, results_station_path)
  
  
  # --- Gather all regional forecasts ---
  res_reg_global <- readRDS(
    paste0("Output_ponctual/Seed_", seed_windows,"/results_regions_global_transform_log_lags_TRUE.RDS")
  )[["global"]]
  
  res_reg <- res_reg_global
  
  regional_foundation_file <- paste0("Output_ponctual/Seed_", seed_windows,"/results_regions_tabICL.csv")
  
  if (file.exists(regional_foundation_file)) {
    
    res_reg_foundation <- read.csv(
      regional_foundation_file,
      sep = ","
    )
    
    res_reg$tabICL <- res_reg_foundation$pred_tabICL
    
  } else {
    
    warning(
      "Foundation file not found: ",
      regional_foundation_file
    )
  }
  forecast_cols <- c(
    "LOCAL_GAM",
    "GLOBAL_RF",
    "GLOBAL_XGB",
    "GLOBAL_GAM",
    "tabICL",
    "LOCAL_RF"
  )
  
  available_cols <- intersect(
    forecast_cols,
    names(res_reg)   # or res_reg / res_nat
  )
  
  res_reg$Combination <- rowMeans(
    res_reg[, available_cols, drop = FALSE],
    na.rm = TRUE
  )
  saveRDS(res_reg, results_regional_path)
  
  
  # --- Gather all national forecasts ---
  res_nat_local <- readRDS(results_national_path)

  res_nat <- res_nat_local
  
  national_foundation_file <- paste0("Output_ponctual/Seed_", seed_windows,"/results_scotland_tabICL.csv")
  
  if (file.exists(national_foundation_file)) {
    
    res_nat_foundation <- read.csv(
      national_foundation_file,
      sep = ","
    )
    
    res_nat$tabICL <- res_nat_foundation$pred_tabICL
    
  } else {
    
    warning(
      "Foundation file not found: ",
      national_foundation_file
    )
  }
  
  forecast_cols <- c(
    "LOCAL_GAM",
    "LOCAL_RF",
    "GLOBAL_RF",
    "GLOBAL_XGB",
    "GLOBAL_GAM",
    "tabICL"
  )
  
  available_cols <- intersect(
    forecast_cols,
    names(res_nat)   # or res_reg / res_nat
  )
  
  res_nat$Combination <- rowMeans(
    res_nat[, available_cols, drop = FALSE],
    na.rm = TRUE
  )
  saveRDS(res_nat, results_national_path)
}

# -----------------------------
# LOAD DATA
# -----------------------------

res_nat <- readRDS(results_national_path)
res_reg <- readRDS(results_regional_path)
res_sta <- readRDS(results_station_path)

H <- readRDS(structural_matrix_path)

region_levels  <- sub("Region_", "", grep("^Region_", names(H), value = TRUE))
station_levels <- sub("Station_", "", grep("^Station_", names(H), value = TRUE))

# -----------------------------
# RECONCILIATION
# -----------------------------

metrics_list <- list()

P_OLS <- refined_OLS_projection(as.matrix(H[, -1]))

df_predictions <- list()

i <- 0

for (window in windows) {
  
  i <- i + 1
  
  print(i)
  
  test_period <- window %>%
    filter(type == "test") %>%
    pull(Date)
  
  train_period <- window %>%
    filter(type == "train") %>%
    pull(Date)
  
  # --- Gather targets and forecasts ---
  
  vectors <- build_vectors(
    period = test_period,
    res_nat = res_nat %>% filter(window_id == i),
    res_reg = res_reg %>% filter(window_id == i),
    res_sta = res_sta %>% filter(window_id == i),
    model_names = model_names
  )
  
  Y <- do.call(rbind, lapply(vectors, `[[`, "y"))
  
  Yhat <- do.call(rbind, lapply(vectors, `[[`, "yhat"))
  
  # --- Build masks for missing values ---
  
  M <- is.na(Y) * 1
  
  masks <- split(M, seq_len(nrow(M)))
  
  # --- Compute refined structural matrices ---
  
  H_t_list <- lapply(masks, function(M_t) {
    
    mask_structural(H, M_t)
  })
  
  # --- Compute refined projection matrices ---
  
  P_OLS_refined <- lapply(H_t_list, function(H_t) {
    
    refined_OLS_projection(H_t)
  })
  
  nT <- nrow(Y)
  
  # --- Temporal loop ---
  
  for (t in seq_len(nT)) {
    
    y    <- Y[t, ]
    yhat <- Yhat[t, ]
    
    y_ols <- as.vector(P_OLS %*% yhat)
    
    y_ols_refined <- as.vector(
      P_OLS_refined[[t]] %*% yhat
    )
    
    df_t <- list(
      y = y,
      yhat = yhat,
      OLS = y_ols,
      OLS_refined = y_ols_refined,
      window = i
    )
    
    df_predictions[[length(df_predictions) + 1]] <- df_t
  }
}

# -----------------------------
# FORMAT OUTPUT
# -----------------------------

df_long <- map_dfr(seq_along(df_predictions), function(k) {
  
  obj <- df_predictions[[k]]
  
  n_nodes <- length(obj$y)
  
  data.frame(
    time = k,
    window = obj$window,
    node = seq_len(n_nodes),
    y = obj$y,
    Direct = obj$yhat,
    OLS = obj$OLS,
    OLS_refined = obj$OLS_refined
  )
})

df_long <- df_long %>%
  mutate(
    node_type = case_when(
      node == 1 ~ "national",
      node %in% 2:33 ~ "regional",
      TRUE ~ "station"
    )
  )

# -----------------------------
# EXPORT RESULTS
# -----------------------------

saveRDS(df_long, output_path)
