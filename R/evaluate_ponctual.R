source("R/forecast_function.R")
source("R/reconciliation_function.R")
source("R/metric_function.R")
source("R/conformal_methods.R")

library(argparser)
library(dplyr)
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
# METRICS
# -----------------------------
df_long <- readRDS("Output_ponctual/Seed_3/reconciled_forecasts_Combination_Combination_Combination.RDS")

tab_global <- data.frame(
  method = c("Direct", "OLS", "OLS_refined"),
  mse = c(
    mse_global(df_long, "Direct"),
    mse_global(df_long, "OLS"),
    mse_global(df_long, "OLS_refined")
  )
)

tab_type <- bind_rows(
  lapply(c("Direct", "OLS", "OLS_refined"), function(m) {
    df <- mse_by_type(df_long, m)
    df$method <- m
    df
  })
)

tab_type_pct <- tab_type %>%
  group_by(node_type) %>%
  mutate(
    mse_direct = mse[method == "Direct"],
    pct_diff = 100 * (mse - mse_direct) / mse_direct
  ) %>%
  ungroup()

tab_node <- bind_rows(
  lapply(c("Direct", "OLS", "OLS_refined"), function(m) {
    df <- mse_by_node(df_long, m)
    df$method <- m
    df
  })
  ) %>% mutate(
    node_type = case_when(
    node == "1" ~ "national",
    # node %in% regional_nodes ~ "regional",
    node %in% 1:(33) ~ "regional", 
    TRUE ~ "station"
    )
  )
tab_level <- tab_node %>% group_by(node_type)
  

methods <- c("Direct", "OLS", "OLS_refined")

metrics_window <- bind_rows(lapply(methods, function(m) {
  
  df_long %>%
    group_by(window) %>%
    summarise(
      Gmse = mse_global(cur_data(), m),
      .groups = "drop"
    ) %>%
    mutate(method = m)
}))
metrics_level_window <- bind_rows(lapply(methods, function(m) {
  
  df_long %>%
    group_by(window, node_type) %>%
    summarise(
      Gmse = mse_global(cur_data(), m),
      .groups = "drop"
    ) %>%
    mutate(method = m)
}))

# -----------------------------
# FIGURES
# -----------------------------

ggplot(metrics_window, aes(x = method, y = Gmse)) +
  
  geom_boxplot(outlier.shape = NA, alpha = 0.2) +
  
  geom_jitter(aes(color = window),
              width = 0.15,
              size = 2,
              alpha = 0.8) +
  
  stat_summary(fun = mean,
               geom = "point",
               color = "black",
               size = 3) +
  
  scale_color_gradientn(colors = rainbow(length(unique(metrics_window$window)))) +
  
  theme_minimal() +
  
  labs(
    x = "",
    y = "mse",
    color = "Window"
  )
plot_mse_level <- function(df, level_name) {
  
  ggplot(df, aes(x = method, y = Gmse)) +
    
    geom_boxplot(outlier.shape = NA, alpha = 0.2) +
    
    geom_jitter(aes(color = window),
                width = 0.15,
                size = 2,
                alpha = 0.8) +
    
    stat_summary(fun = mean,
                 geom = "point",
                 color = "black",
                 size = 3) +
    
    scale_color_gradientn(colors = rainbow(length(unique(df$window)))) +
    
    theme_minimal() +
    
    labs(
      title = level_name,
      x = "",
      y = "mse",
      color = "Window"
    )
}
df_nat <- metrics_level_window %>% filter(node_type == "national")
df_reg <- metrics_level_window %>% filter(node_type == "regional")
df_sta <- metrics_level_window %>% filter(node_type == "station")

plot_mse_level(df_nat, "National")
plot_mse_level(df_reg, "Regional")
plot_mse_level(df_sta, "Station")
