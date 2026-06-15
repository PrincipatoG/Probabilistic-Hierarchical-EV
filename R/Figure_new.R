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

results_dir <- "Output_conformal"
figures_dir <- "Figures_new"

alpha <- 0.1

dir.create(figures_dir, showWarnings = FALSE)

# =====================================================
# LOAD ALL RESULTS
# =====================================================

result_files <- list.files(
  results_dir,
  pattern = "\\.RDS$",
  full.names = TRUE
)

result_files <- result_files[
  !grepl("forecast", basename(result_files), ignore.case = TRUE)
]
result_files <- result_files[
  !grepl("result", basename(result_files), ignore.case = TRUE)
]
result_files <- result_files[
  !grepl("number", basename(result_files), ignore.case = TRUE)
]
cat("Detected result files:\n")
print(basename(result_files))

summary_list <- list()
time_list <- list()

for (f in result_files) {
  
  obj <- readRDS(f)
  
  if (!is.list(obj)) next
  
  if (!is.null(obj$summary)) {
    
    tmp <- obj$summary
    tmp$file_source <- basename(f)
    
    summary_list[[length(summary_list) + 1]] <- tmp
  }
  
  if (!is.null(obj$time)) {
    
    tmp <- obj$time
    tmp$file_source <- basename(f)
    
    time_list[[length(time_list) + 1]] <- tmp
  }
}

df_summary <- bind_rows(summary_list)
df_time <- bind_rows(time_list)

# =====================================================
# INFER CONFIGURATION METADATA
# =====================================================

infer_conformal_method <- function(x) {
  
  case_when(
    str_detect(x, "Adaptive_CP_MNR") ~ "Adaptive_CP_MNR",
    str_detect(x, "Nested") ~ "CP_MNR_Nested_star",
    str_detect(x, "naive") ~ "CP_MNR_naive",
    str_detect(x, "CP_MNR") ~ "CP_MNR",
    TRUE ~ "Unknown"
  )
}

infer_projection_method <- function(x) {
  
  case_when(
    str_detect(x, "refined OLS") ~ "refined OLS",
    str_detect(x, "OLS") ~ "OLS",
    str_detect(x, "Direct") ~ "Direct",
    TRUE ~ x
  )
}

df_summary <- df_summary %>%
  mutate(
    conformal_method = infer_conformal_method(file_source),
    projection_method = infer_projection_method(method),
    experiment = paste(
      conformal_method,
      projection_method,
      sep = " | "
    )
  )
# # --- Tests 
# df_nat <- df_summary %>% filter(node_type == "national")
# df_reg <- df_summary %>% filter(node == 464)

# =====================================================
# WEIGHTED QUANTILES
# =====================================================

weighted_quantile <- function(
    x,
    w,
    probs = c(0.1, 0.9)
) {
  
  ok <- !(is.na(x) | is.na(w))
  
  x <- x[ok]
  w <- w[ok]
  
  ord <- order(x)
  
  x <- x[ord]
  w <- w[ord]
  
  cw <- cumsum(w) / sum(w)
  
  sapply(
    probs,
    function(p) {
      x[which(cw >= p)[1]]
    }
  )
}

# =====================================================
# METHOD-LEVEL SUMMARY
# =====================================================

method_summary <- function(df) {
  
  df %>%
    
    mutate(
      
      # Rename conformal algorithms
      conformal_method = case_when(
        
        conformal_method == "CP_MNR_naive" ~ "CP",
        
        conformal_method == "CP_MNR" ~ "CP MNR",
        
        conformal_method == "Adaptive_CP_MNR" ~ "ACI MNR",
        
        conformal_method == "CP_MNR_Nested_star" ~ "CP MNR Nested*",
        
        TRUE ~ conformal_method
      ),
      
      # Nested corresponds to refined OLS
      projection_method = case_when(
        
        projection_method == "Nested" ~ "rOLS",
        
        projection_method == "refined OLS" ~ "rOLS",
        
        TRUE ~ projection_method
      )
    ) %>%
    
    group_by(
      node_type,
      conformal_method,
      projection_method
    ) %>%
    
    summarise(
      
      # Weighted averages
      coverage_mean = weighted.mean(
        coverage,
        w = n_eff,
        na.rm = TRUE
      ),
      
      length_mean = weighted.mean(
        length,
        w = n_eff,
        na.rm = TRUE
      ),
      
      # Weighted interquartile ranges
      coverage_q25 = weighted_quantile(
        coverage,
        n_eff,
        probs = 0.25
      )[1],
      
      coverage_q75 = weighted_quantile(
        coverage,
        n_eff,
        probs = 0.75
      )[1],
      
      length_q25 = weighted_quantile(
        length,
        n_eff,
        probs = 0.25
      )[1],
      
      length_q75 = weighted_quantile(
        length,
        n_eff,
        probs = 0.75
      )[1],
      
      .groups = "drop"
    )
}
# ===================================================== 
# METRICS 
# =====================================================
node_metrics <- function( df, min_obs = 14 ) { 
  group_vars <- c( "node", 
                   "node_type",
                   "conformal_method", 
                   "projection_method" )
  # garder model_tag si présent 
  if ("model_tag" %in% names(df)) { 
    group_vars <- c(group_vars, "model_tag") 
    } 
  
  # garder experiment si présent 
  if ("experiment" %in% names(df)) { 
    group_vars <- c(group_vars, "experiment") 
    }
  
  df %>% dplyr::group_by(across(all_of(group_vars))) %>% 
    dplyr::summarise( coverage = weighted.mean( coverage, w = n_eff, na.rm = TRUE ), 
                      length = weighted.mean( length, w = n_eff, na.rm = TRUE ), 
                      n_eff = sum(n_eff, na.rm = TRUE), 
                      .groups = "drop" ) 
  } 

df_metrics <- node_metrics( df_summary, min_obs = 14 )
  
df_type <- type_metrics(df_summary)  
  
# =====================================================
# PUBLICATION THEME
# =====================================================

theme_aoas <- theme_minimal(base_size = 15) +
  
  theme(
    
    panel.grid.minor = element_blank(),
    
    panel.grid.major = element_blank(),
    
    axis.line = element_line(
      linewidth = 0.35,
      colour = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.35
    ),
    
    strip.text = element_text(
      face = "bold",
      size = 13
    ),
    
    panel.spacing.y = unit(0.15, "cm"),
    
    legend.title = element_blank(),
    
    legend.position = "bottom",
    
    legend.box = "vertical",
    
    legend.spacing.y = unit(0.1, "cm"),
    
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  )
# =====================================================
# VISUAL SETTINGS
# =====================================================

algo_colors <- c(
  "CP" = "#4E79A7",
  "CP MNR" = "#F28E2B",
  "ACI MNR" = "#59A14F",
  "CP MNR Nested*" = "#E15759"
)

algo_shapes <- c(
  "CP" = 16,
  "CP MNR" = 17,
  "ACI MNR" = 15,
  "CP MNR Nested*" = 18
)

projection_colors <- c(
  "Direct" = "#6D597A",   
  "OLS"    = "#F4A261",   
  "rOLS"   = "#2A9D8F" 
)

projection_shapes <- c(
  "OLS" = 16,
  "rOLS" = 17,
  "Direct" = 15
)
conformal_order <- c(
  "CP MNR Nested*",
  "ACI MNR",
  "CP MNR",
  "CP"
)
# =====================================================
# COMPARE CONFORMAL METHODS FOR rOLS
# =====================================================

plot_algorithms_rols <- function(
    df,
    node_type_value,
    y_limits,
    alpha = 0.1
) {
  
  df_plot <- df %>%
    
    filter(
      node_type == node_type_value,
      projection_method == "rOLS"
    )
  
  y_lim <- get_y_limits(
    node_type_value,
    y_limits
  )
  
  ggplot(
    df_plot,
    aes(
      x = coverage_mean,
      y = length_mean,
      colour = conformal_method,
      shape = conformal_method
    )
  ) +
    
    # Horizontal uncertainty bars
    geom_errorbarh(
      aes(
        xmin = coverage_q25,
        xmax = coverage_q75
      ),
      height = 0.015 * diff(y_lim),
      linewidth = 0.5,
      linetype = "dashed"
    ) +
    
    # Vertical uncertainty bars
    geom_errorbar(
      aes(
        ymin = length_q25,
        ymax = length_q75
      ),
      width = 0.003,
      linewidth = 0.5,
      linetype = "dashed"
    ) +
    
    # Central points
    geom_point(
      size = 5
    ) +
    
    geom_vline(
      xintercept = 1 - alpha,
      linetype = "dashed",
      colour = "grey60",
      linewidth = 0.5
    ) +
    
    scale_colour_manual(
      values = algo_colors
    ) +
    
    scale_shape_manual(
      values = algo_shapes
    ) +
    
    coord_cartesian(
      xlim = get_x_limits(node_type_value),
      ylim = y_lim
    ) +
    
    labs(
      x = "Coverage",
      y = "Length",
      title = paste(
        node_type_value,
        "- rOLS"
      )
    ) +
    
    theme_aoas
}
# =====================================================
# COMPARE PROJECTION METHODS
# =====================================================

plot_projections <- function(
    df,
    node_type_value,
    conformal_value,
    y_limits,
    alpha = 0.1
) {
  
  df_plot <- df %>%
    
    filter(
      node_type == node_type_value,
      conformal_method == conformal_value
    )
  
  y_lim <- get_y_limits(
    node_type_value,
    y_limits
  )
  
  ggplot(
    df_plot,
    aes(
      x = coverage_mean,
      y = length_mean,
      colour = projection_method,
      shape = projection_method
    )
  ) +
    
    geom_errorbarh(
      aes(
        xmin = coverage_q25,
        xmax = coverage_q75
      ),
      height = 0.015 * diff(y_lim),
      linewidth = 0.5,
      linetype = "dashed"
    ) +
    
    geom_errorbar(
      aes(
        ymin = length_q25,
        ymax = length_q75
      ),
      width = 0.003,
      linewidth = 0.5,
      linetype = "dashed"
    ) +
    
    geom_point(
      size = 5
    ) +
    
    geom_vline(
      xintercept = 1 - alpha,
      linetype = "dashed",
      colour = "grey60",
      linewidth = 0.5
    ) +
    
    scale_colour_manual(
      values = projection_colors
    ) +
    
    scale_shape_manual(
      values = projection_shapes
    ) +
    
    coord_cartesian(
      xlim = get_x_limits(node_type_value),
      ylim = y_lim
    ) +
    
    labs(
      x = "Coverage",
      y = "Length",
      title = paste(
        node_type_value,
        "-",
        conformal_value
      )
    ) +
    
    theme_aoas
}

# =====================================================
# COVERAGE-LENGTH SUMMARY PLOT
# =====================================================
df_summary <- method_summary(df_metrics)

plot_cov_len_summary <- function(
    df,
    node_type_value,
    alpha = 0.1
) {
  
  df_plot <- df %>%
    filter(node_type == node_type_value)
  
  ggplot(
    df_plot,
    aes(
      x = coverage_mean,
      y = length_mean,
      colour = conformal_method,
      shape = projection_method
    )
  ) +
    
    # Horizontal uncertainty bars
    geom_errorbarh(
      aes(
        xmin = coverage_q25,
        xmax = coverage_q75
      ),
      height = 0,
      linewidth = 0.8,
      alpha = 0.35
    ) +
    
    # Vertical uncertainty bars
    geom_errorbar(
      aes(
        ymin = length_q25,
        ymax = length_q75
      ),
      width = 0,
      linewidth = 0.8,
      alpha = 0.35
    ) +
    
    # Central points
    geom_point(
      size = 4.5,
      stroke = 1.1
    ) +
    
    # Target coverage
    geom_vline(
      xintercept = 1 - alpha,
      linetype = "dashed",
      colour = "grey55",
      linewidth = 0.5
    ) +
    
    labs(
      x = "Coverage",
      y = "Average interval length"
    ) +
    
    # Elegant journal-style palette
    scale_colour_manual(
      values = c(
        "CP" = "#4C78A8",
        "CP MNR" = "#F58518",
        "ACI MNR" = "#54A24B",
        "CP MNR Nested*" = "#B279A2"
      )
    ) +
    
    # Distinguishable publication-quality shapes
    scale_shape_manual(
      values = c(
        "OLS" = 16,
        "rOLS" = 17,
        "Direct" = 15
      )
    ) +
    
    guides(
      colour = guide_legend(
        order = 1,
        override.aes = list(size = 4)
      ),
      
      shape = guide_legend(
        order = 2
      )
    ) +
    
    theme_aoas
}

# =====================================================
# GENERATE FIGURES
# =====================================================

df_summary <- method_summary(df_metrics)

# =====================================================
# Y-AXIS LIMITS BY NODE TYPE
# =====================================================

y_limits <- df_summary %>%
  
  group_by(node_type) %>%
  
  summarise(
    
    y_min = min(length_q25, na.rm = TRUE),
    
    y_max = max(length_q75, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  
  mutate(
    
    # Add small margins
    
    y_range = y_max - y_min,
    
    y_min = y_min - 0.05 * y_range,
    
    y_max = y_max + 0.05 * y_range
  )

# =====================================================
# AXIS SETTINGS
# =====================================================

get_y_limits <- function(
    node_type_value,
    y_limits
) {
  
  row <- y_limits %>%
    filter(node_type == node_type_value)
  
  c(row$y_min, row$y_max)
}
get_x_limits <- function(node_type_value) {
  
  if (node_type_value %in% c("national", "regional")) {
    return(c(0.8, 1))
  } else {
    return(c(0.6, 1))
  }
}
# ----------------------------------
# Figure 1: algorithms under rOLS
# ----------------------------------
node_types <- unique(df_metrics$node_type)

plots_algorithms <- lapply(
  node_types,
  function(nt) {
    
    plot_algorithms_rols(
      df_summary,
      node_type_value = nt,
      y_limits = y_limits,
      alpha = alpha
    )
  }
)

names(plots_algorithms) <- node_types

plots_algorithms$station
plots_algorithms$regional
plots_algorithms$national

# ----------------------------------
# Figure 2: projections by algorithm
# ----------------------------------

conformal_methods <- unique(
  df_summary$conformal_method
)

plots_projection <- list()

for (cm in conformal_methods) {
  
  plots_projection[[cm]] <- lapply(
    node_types,
    function(nt) {
      
      plot_projections(
        df_summary,
        node_type_value = nt,
        conformal_value = cm,
        y_limits = y_limits,
        alpha = alpha
      )
    }
  )
  
  names(plots_projection[[cm]]) <- node_types
}

# =====================================================
# SAVE FIGURES
# =====================================================

# --- Algorithms under rOLS

for (nt in node_types) {
  
  ggsave(
    filename = file.path(
      figures_dir,
      paste0(nt, "_algorithms_rOLS.pdf")
    ),
    
    plot = plots_algorithms[[nt]],
    
    width = 7,
    height = 5
  )
}

# --- Projection comparisons

for (cm in conformal_methods) {
  
  for (nt in node_types) {
    
    ggsave(
      filename = file.path(
        figures_dir,
        paste0(
          nt,
          "_",
          gsub("[^A-Za-z0-9]", "_", cm),
          "_projection.pdf"
        )
      ),
      
      plot = plots_projection[[cm]][[nt]],
      
      width = 7,
      height = 5
    )
  }
}


# =====================================================
# CLEAN HISTOGRAM BLOCK (ADD-ON)
# =====================================================

# safety: alias if missing (prevents crash)
theme_clean <- theme_aoas


# =====================================================
# PREPROCESSING (same logic as before)
# =====================================================

preprocess_methods <- function(df) {
  
  df %>%
    mutate(
      
      conformal_method = case_when(
        conformal_method == "CP_MNR_naive" ~ "CP",
        conformal_method == "CP_MNR" ~ "CP MNR",
        conformal_method == "Adaptive_CP_MNR" ~ "ACI MNR",
        conformal_method == "CP_MNR_Nested_star" ~ "CP MNR Nested*",
        TRUE ~ conformal_method
      ),
      
      projection_method = case_when(
        projection_method %in% c("Nested", "refined OLS") ~ "rOLS",
        TRUE ~ projection_method
      ),
      
      conformal_method = factor(conformal_method, levels = conformal_order)
    )
}

# =====================================================
# COVERAGE HISTOGRAM (same style as earlier working version)
# =====================================================
algo_colors_dark <- c(
  "CP" = "#2F5D8A",
  "CP MNR" = "#C46E18",
  "ACI MNR" = "#2E7D32",
  "CP MNR Nested*" = "#8B0000"
)

plot_coverage_hist <- function(df, node_type_value, bin_width = 0.02, alpha = 0.1) {
  
  df_plot <- df %>%
    filter(node_type == node_type_value) %>%
    preprocess_methods()
  
  df_mean <- df_plot %>%
    group_by(conformal_method, projection_method) %>%
    summarise(
      mean_cov = weighted.mean(coverage, n_eff, na.rm = TRUE),
      .groups = "drop"
    )
  limit_figure <- get_x_limits(node_type_value)
  ggplot(df_plot, aes(x = coverage, weight = n_eff, fill = conformal_method)) +
    
    geom_histogram(
      aes(y = after_stat(density)),
      binwidth = bin_width,
      alpha = 0.75,
      colour = "white",
      linewidth = 0.3
    ) +
    
    facet_grid(
      conformal_method ~ projection_method,
      drop = FALSE,
      scales = "free_y",
      space = "free_y"#,
      #switch = "y"
    ) +
    
    geom_vline(
      xintercept = 1 - alpha,
      linetype = "dashed",
      colour = "grey55",
      linewidth = 0.5
    ) +
    
    geom_vline(
      data = df_mean,
      aes(
        xintercept = mean_cov,
        colour = conformal_method
      ),
      linewidth = 1.1,
      alpha = 0.95,
      show.legend = FALSE
    ) +
    scale_fill_manual(
      values = algo_colors,
      drop = FALSE,
      breaks = c(
        "CP",
        "CP MNR",
        "ACI MNR",
        "CP MNR Nested*"
        ),
      name = ""
    ) +
    scale_colour_manual(
      values = algo_colors_dark,
      guide = "none"
    ) +
    
    scale_x_continuous(
      limits = limit_figure,
      breaks = seq(min(limit_figure), 1, by = 0.1)
    ) +
    
    scale_y_continuous(
      breaks = function(x) {
        pretty(
          c(0, ceiling(max(x) / 10) * 10),
          n = 3
        )
      }
    ) +
    
    labs(
      x = "Coverage",
      y = "Density",
      title = paste(node_type_value)
    ) +
    
    # guides(
    #   fill = guide_legend(
    #     title = "",
    #     ncol = 1
    #   )
    # ) +
    
    theme_aoas +
    
    theme(
      
      # # REMOVE ROW LABELS
      # strip.text.y = element_blank(),
      
      # RESTORE ROW LABELS (methods on the side)
      strip.text.y = element_text(
        face = "bold",
        size = 12,
        angle = 0
      ),
      
      strip.placement = "outside",
      
      # Keep column labels
      strip.text.x = element_text(
        face = "bold",
        size = 12
      ),
      
      legend.title = element_text(
        face = "bold"
      ),
      
      panel.spacing.y = unit(0.15, "cm"),
      panel.spacing.x = unit(0.25, "cm")
    )
}


# =====================================================
# AXES SETTINGS (UNIQUE SCALES = SAME X-RANGE ACROSS NODE TYPES)
# =====================================================

bin_widths <- c(
  national = 300,
  regional = 300,
  station = 10
)

node_types <- sort(unique(df_metrics$node_type))


# =====================================================
# GENERATE HISTOGRAMS
# =====================================================

plots_cov <- setNames(
  lapply(node_types, function(nt) {
    plot_coverage_hist(
      df_metrics,
      node_type_value = nt,
      bin_width = 0.02,
      alpha = alpha
    )
  }),
  node_types
)
plots_cov$station

# =====================================================
# SAVE HISTOGRAMS
# =====================================================

for (nt in node_types) {
  
  ggsave(
    filename = file.path(figures_dir, paste0(nt, "_coverage_hist.pdf")),
    plot = plots_cov[[nt]],
    width = 12,
    height = 4
  )
  
}

plot_coverage_hist_rols <- function(
    df,
    node_type_value,
    bin_width = 0.02,
    alpha = 0.1
) {
  
  df_plot <- df %>%
    filter(
      node_type == node_type_value,
      projection_method %in% c("refined OLS", "Nested")
    ) %>%
    preprocess_methods() %>%
    mutate(coverage = pmax(coverage, 0.65),
           conformal_method = factor(
             conformal_method,
             levels = c("CP", "CP MNR", "ACI MNR", "CP MNR Nested*")
           ))
  
  df_mean <- df_plot %>%
    group_by(conformal_method) %>%
    summarise(
      mean_cov = weighted.mean(coverage, n_eff, na.rm = TRUE),
      .groups = "drop"
    )
  limit_figure <- get_x_limits(node_type_value)
  
  ggplot(
    df_plot,
    aes(
      x = coverage,
      weight = n_eff,
      fill = conformal_method
    )
  ) +
    
    geom_histogram(
      aes(y = after_stat(density)),
      binwidth = bin_width,
      colour = "white",
      linewidth = 0.3,
      alpha = 0.85
    ) +
    
    facet_wrap(
      ~ conformal_method,
      nrow = 1
    ) +
    
    geom_vline(
      xintercept = 1 - alpha,
      linetype = "dashed",
      colour = "grey50",
      linewidth = 0.5
    ) +
    
    geom_segment(
      data = df_mean,
      inherit.aes = FALSE,
      aes(
        x = mean_cov,
        xend = mean_cov,
        y = -Inf,
        yend = 0,
        colour = conformal_method
      ),
      linewidth = 1.1,
      show.legend = FALSE
    ) +
    scale_colour_manual(
      values = algo_colors_dark,
      guide = "none"
    ) +
    
    scale_fill_manual(values = algo_colors) +
    
    scale_x_continuous(
      limits = c(limit_figure, 1),
      breaks = seq(min(limit_figure), 1, 0.1)
    ) +
    
    scale_y_continuous(
      limits = NULL  # IMPORTANT: shared autoscale
    ) +
    
    labs(
      x = "Coverage",
      y = "Density"
    ) +
    
    guides(
      fill = guide_legend(title = "")
    ) +
    
    theme_aoas +
    
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      
      strip.text = element_blank(),
      strip.background = element_blank(),
      
      panel.spacing.x = unit(0.5, "cm")
    )
}
# =====================================================
# GENERATE HISTOGRAMS
# =====================================================

plots_cov_rols <- setNames(
  lapply(node_types, function(nt) {
    plot_coverage_hist_rols(
      df_metrics,
      node_type_value = nt,
      bin_width = 0.02,
      alpha = alpha
    )
  }),
  node_types
)
plots_cov_rols$station

# =====================================================
# SAVE HISTOGRAMS
# =====================================================

for (nt in node_types) {
  
  ggsave(
    filename = file.path(figures_dir, paste0(nt, "_coverage_hist_rols.pdf")),
    plot = plots_cov_rols[[nt]],
    width = 12,
    height = 4
  )
  
}
