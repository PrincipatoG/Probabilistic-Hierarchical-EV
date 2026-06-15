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


# --- Load the result tables ---
df_all <- readRDS("Output/new_result_CP_all_but_Nested_MNR_all_windows.RDS")
df_naive <- readRDS("Output/new_result_CP_naive_all_windows.RDS")
df_naive$method <- paste0(df_naive$method , "_naive")

df_ACI <- readRDS("Output/result_Adaptive_CP_MNR_all_windows.RDS")

df_all <- df_all %>% rbind(df_naive)

# --- Get the metric of interests ---
df_metrics <- node_metrics(df_all, min_obs = 14)
# df_global  <- global_metrics(df_all)
df_type    <- type_metrics(df_all)

df_metrics_ACI <- node_metrics(df_ACI, min_obs = 14)
# df_global_ACI  <- global_metrics(df_ACI)
df_type_ACI    <- type_metrics(df_ACI)

# -----------------------------
# Figures
# -----------------------------
theme_clean <- theme_minimal(base_size = 18) +
  theme(
    panel.grid = element_blank(),     # pas de grille
    axis.line = element_line(),       # axes visibles
    axis.ticks = element_line(),      # ticks visibles
    plot.title = element_text(face = "bold")
  )
# --- Cov/Len Plots ---
plots <- plot_by_type(df_all, alpha = 0.1)

plots$national
plots$regional
plots$station

# --- Save the figures ---
ggsave("Figures/new_national_cov_len.pdf",plots$national)
ggsave("Figures/new_regional_cov_len.pdf",plots$regional)
ggsave("Figures/new_station_cov_len.pdf",plots$station)

plots_ACI <- plot_by_type(df_ACI, alpha = 0.1)

plots_ACI$national
plots_ACI$regional
plots_ACI$station

# --- Save the figures ---
ggsave("Figures/aci_national_cov_len.pdf",plots_ACI$national)
ggsave("Figures/aci_regional_cov_len.pdf",plots_ACI$regional)
ggsave("Figures/aci_station_cov_len.pdf",plots_ACI$station)

# --- Cov CDF Plots ---
plot_coverage_hist <- function(df, df_type, node_type_value, bin_width = 0.02, alpha = 0.1) {
  
  df_plot <- df %>%
    filter(node_type == node_type_value)
  
  df_vline <- df_type %>%
    filter(node_type == node_type_value)
  
  ggplot(df_plot, aes(x = coverage, weight = n_eff, fill = method)) +
    
    geom_histogram(
      aes(y = after_stat(density)),
      binwidth = bin_width,
      alpha = 0.7,
      color = "black"
    ) +
    
    facet_wrap(~ method, ncol = 1, scales = "fixed") +
    
    # vertical line per method
    geom_vline(
      data = df_vline,
      aes(xintercept = coverage, color = method),
      linetype = "dashed",
      linewidth = 0.8,
      show.legend = FALSE
    ) +
    # target coverage
    geom_vline(xintercept = 1 - alpha,
               linetype = "dashed",
               color = "grey60") +
    coord_cartesian(xlim = c(0, 1)) +
    
    labs(
      x = "Coverage",
      y = "CDF"
    ) +
    
    scale_fill_brewer(palette = "Set2") +
    scale_color_brewer(palette = "Set2") +
    
    theme_minimal() +
    
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold"),
      # --- GRID CLEANING ---
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.4),
      
      panel.background = element_blank(),
      plot.background = element_blank()
    )
}
unique_node_types <- unique(df_metrics$node_type)

plots_cov <- lapply(unique_node_types, function(nt) {
  plot_coverage_hist(df_metrics, df_type, nt)
})

names(plots_cov) <- unique_node_types

plots_cov$national
plots_cov$regional
plots_cov$station

# --- Save the figures ---
ggsave("Figures/national_cov.pdf",plots_cov$national)
ggsave("Figures/regional_cov.pdf",plots_cov$regional)
ggsave("Figures/station_cov.pdf",plots_cov$station)

unique_node_types_ACI <- unique(df_metrics_ACI$node_type)

plots_cov_ACI <- lapply(unique_node_types, function(nt) {
  plot_coverage_hist(df_metrics_ACI, df_type_ACI, nt)
})

names(plots_cov_ACI) <- unique_node_types_ACI

plots_cov_ACI$national
plots_cov_ACI$regional
plots_cov_ACI$station

# --- Save the figures ---
ggsave("Figures/national_cov_ACI.pdf",plots_cov_ACI$national)
ggsave("Figures/regional_cov_ACI.pdf",plots_cov_ACI$regional)
ggsave("Figures/station_cov_ACI.pdf",plots_cov_ACI$station)

# --- Len CDF Plots ---
plot_length_hist <- function(df, df_type, node_type_value, bin_width = 500, alpha = 0.1) {
  
  df_plot <- df %>%
    filter(node_type == node_type_value)
  
  df_vline <- df_type %>%
    filter(node_type == node_type_value)
  
  ggplot(df_plot, aes(x = length, weight = n_eff, fill = method)) +
    
    geom_histogram(
      aes(y = after_stat(density)),
      binwidth = bin_width,
      alpha = 0.7,
      color = "black"
    ) +
    
    facet_wrap(~ method, ncol = 1, scales = "fixed") +
    
    # vertical line per method
    geom_vline(
      data = df_vline,
      aes(xintercept = length, color = method),
      linetype = "dashed",
      linewidth = 0.8,
      show.legend = FALSE
    ) +
    labs(
      x = "Length",
      y = "CDF"
    ) +
    
    scale_fill_brewer(palette = "Set2") +
    scale_color_brewer(palette = "Set2") +
    
    theme_minimal() +
    
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold"),
      # --- GRID CLEANING ---
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.4),
      
      panel.background = element_blank(),
      plot.background = element_blank()
    )
}
unique_node_types <- unique(df_metrics$node_type)

bin_widths <- c(300, 300, 10)
names(bin_widths) <- unique_node_types

plots_len <- lapply(unique_node_types, function(nt) {
  plot_length_hist(df_metrics, df_type, nt, bin_widths[nt])
})

names(plots_len) <- unique_node_types

plots_len$national
plots_len$regional
plots_len$station

# --- Save the figures ---
ggsave("Figures/national_len.pdf",plots_len$national)
ggsave("Figures/regional_len.pdf",plots_len$regional)
ggsave("Figures/station_len.pdf",plots_len$station)

unique_node_types_ACI <- unique(df_metrics_ACI$node_type)

plots_len_ACI <- lapply(unique_node_types, function(nt) {
  plot_length_hist(df_metrics_ACI, df_type_ACI, nt,  bin_widths[nt])
})

names(plots_len_ACI) <- unique_node_types_ACI

plots_len_ACI$national
plots_len_ACI$regional
plots_len_ACI$station

# --- Save the figures ---
ggsave("Figures/national_len_ACI.pdf",plots_len_ACI$national)
ggsave("Figures/regional_len_ACI.pdf",plots_len_ACI$regional)
ggsave("Figures/station_len_ACI.pdf",plots_len_ACI$station)

