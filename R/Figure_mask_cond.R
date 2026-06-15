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
# LOAD RAW DATA
# =====================================================
# model_names <- list("LOCAL_GAM", "tabICL", "tabICL")
model_names <- list("LOCAL_GAM", "GLOBAL_RF", "GLOBAL_RF")

H <- readRDS("Data/structural.RDS")

res_nat <- readRDS("results/results_scotland.RDS")
res_reg <- readRDS("results/results_regions.RDS")
res_sta <- readRDS("results/results_stations.RDS")

test_period <- res_nat %>%
  filter(type == "test") %>%
  pull(Date)

vectors_test <- build_vectors(
  period = test_period,
  res_nat = res_nat %>% filter(type == "test"),
  res_reg = res_reg %>% filter(type == "test"),
  res_sta = res_sta %>% filter(type == "test"),
  model_names = model_names
)

Y <- do.call(
  rbind,
  lapply(vectors_test, `[[`, "y")
)

M <- is.na(Y) * 1

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

df_time <- df_time %>%
  mutate(
    conformal_method = infer_conformal_method(file_source),
    projection_method = infer_projection_method(method),
    experiment = paste(
      conformal_method,
      projection_method,
      sep = " | "
    ),
    node_type = case_when(
      node == "Y_nat" ~ "national",
      node %in% 1:(ncol(H) - 1) ~ "station", 
      TRUE ~ "regional"
    )
  )

# =====================================================
# MASK CONDITIONAL METRICS (QUANTILE-BASED VERSION)
# =====================================================

get_neighbors <- function(i, H){
  
  region_i <- which(H[i, 3:35] == 1)
  
  neighbors <- which(H[, region_i + 2] == 1)
  
  setdiff(neighbors, i)
}


extract_similar_masking <- function(i, M, H, test_period){
  
  # Filter instants where node i is active
  test_idx <- which(rownames(M) %in% as.character(test_period))
  active_idx <- test_idx[M[test_idx, 1 + 32 + i] == 0]
  
  if(length(active_idx) == 0){
    return(NULL)
  }
  
  # Neighbors
  neighbors <- get_neighbors(i, H)
  
  # Number of active neighbors at each time
  n_active_neighbors <- sapply(
    active_idx,
    function(t) sum(M[t, 1 + 32 + neighbors] == 0)
  )
  
  # =====================================================
  # QUANTILE-BASED STRUCTURE SPLIT
  # =====================================================
  
  q1 <- quantile(n_active_neighbors, 0.25, na.rm = TRUE)
  q3 <- quantile(n_active_neighbors, 0.75, na.rm = TRUE)
  
  low_idx  <- active_idx[n_active_neighbors <= q1]
  mid_idx  <- active_idx[n_active_neighbors > q1 & n_active_neighbors < q3]
  high_idx <- active_idx[n_active_neighbors >= q3]
  
  list(
    low_structure  = rownames(M)[low_idx],
    mid_structure  = rownames(M)[mid_idx],
    high_structure = rownames(M)[high_idx],
    n_active_neighbors = n_active_neighbors
  )
}

extract_similar_masking_bis <- function(i, M, H, test_period){
  
  # Filter instants where node i is active
  test_idx <- which(rownames(M) %in% as.character(test_period))
  active_idx <- test_idx[M[test_idx, 1 + 32 + i] == 0]
  
  if(length(active_idx) == 0){
    return(NULL)
  }
  
  # Neighbors consists in all nodes here
  neighbors <- unique(H$Stations)
  
  # Number of active neighbors at each time
  n_active_neighbors <- sapply(
    active_idx,
    function(t) sum(M[t, 1 + 32 + neighbors] == 0)
  )
  
  # =====================================================
  # QUANTILE-BASED STRUCTURE SPLIT
  # =====================================================
  
  q1 <- quantile(n_active_neighbors, 0.25, na.rm = TRUE)
  q3 <- quantile(n_active_neighbors, 0.75, na.rm = TRUE)
  
  low_idx  <- active_idx[n_active_neighbors <= q1]
  mid_idx  <- active_idx[n_active_neighbors > q1 & n_active_neighbors < q3]
  high_idx <- active_idx[n_active_neighbors >= q3]
  
  list(
    low_structure  = rownames(M)[low_idx],
    mid_structure  = rownames(M)[mid_idx],
    high_structure = rownames(M)[high_idx],
    n_active_neighbors = n_active_neighbors
  )
}


# =====================================================
# BUILD MASK DATAFRAME
# =====================================================

mask_df <- bind_rows(lapply(seq_len(nrow(H)), function(i){
  
  mask <- extract_similar_masking(i, M, H, test_period)
  # mask <- extract_similar_masking_bis(i, M, H, test_period)
  
  if (is.null(mask)) return(NULL)
  
  data.frame(
    node = i,
    date = c(
      mask$low_structure,
      mask$mid_structure,
      mask$high_structure
    ),
    structure_type = c(
      rep("low",  length(mask$low_structure)),
      rep("mid",  length(mask$mid_structure)),
      rep("high", length(mask$high_structure))
    )
  )
}))

mask_df <- mask_df %>%
  mutate(
    node = as.character(node),
    date = as.Date(date)
  )


# =====================================================
# COVERAGE ESTIMATION
# =====================================================

final_df <- df_time %>%
  filter(
    node_type == "station",
    projection_method %in% c("refined OLS", "Nested")
  ) %>%
  inner_join(mask_df, by = c("node", "date")) %>%
  group_by(
    node, node_type,
    conformal_method,
    projection_method,
    experiment,
    structure_type
  ) %>%
  summarise(
    coverage = mean(covered),
    .groups = "drop"
  )


# =====================================================
# PLOT 1: DISTRIBUTION
# =====================================================

g1 <- ggplot(final_df, aes(x = structure_type, y = coverage, fill = structure_type)) +
  geom_violin(alpha = 0.4) +
  geom_boxplot(width = 0.15, outlier.alpha = 0.2) +
  facet_wrap(~ conformal_method) +
  theme_minimal() +
  labs(
    title = "Coverage distribution by structure quantiles",
    x = "",
    y = "Coverage"
  )
g1 
ggsave("Figures_cond/violin_regional_neighbors.pdf", g1)
# =====================================================
# NODE SUMMARY
# =====================================================

node_summary <- final_df %>%
  group_by(node, conformal_method, structure_type) %>%
  summarise(coverage = mean(coverage), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = structure_type,
    values_from = coverage
  )


# =====================================================
# METHOD SUMMARY
# =====================================================

method_summary <- node_summary %>%
  group_by(conformal_method) %>%
  summarise(
    mean_low  = mean(low, na.rm = TRUE),
    mean_mid  = mean(mid, na.rm = TRUE),
    mean_high = mean(high, na.rm = TRUE),
    sd_low    = sd(low, na.rm = TRUE),
    sd_mid    = sd(mid, na.rm = TRUE),
    sd_high   = sd(high, na.rm = TRUE),
    .groups = "drop"
  )


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
# PLOT 2: CALIBRATION MAP
# =====================================================

method_summary <- method_summary %>%
  mutate(
    conformal_method = dplyr::recode(
      conformal_method,
      "CP_MNR_naive" = "CP",
      "CP_MNR" = "CP MNR",
      "Adaptive_CP_MNR" = "ACI MNR",
      "CP_MNR_Nested_star" = "CP MNR Nested*"
    ),
    conformal_method = factor(
      conformal_method,
      levels = rev(conformal_order)
    )
  )

g2 <- ggplot(
  method_summary,
  aes(
    x = mean_low,
    y = mean_high,
    color = conformal_method,
    shape = conformal_method
  )
) +
  
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    linewidth = 0.4,
    colour = "grey50"
  ) +
  
  geom_point(
    size = 4,
    alpha = 0.9
  ) +
  
  scale_color_manual(
    values = algo_colors,
    breaks = rev(conformal_order)
  ) +
  
  scale_shape_manual(
    values = algo_shapes,
    breaks = rev(conformal_order)
  ) +
  
  coord_equal() +
  
  labs(
    title = "",
    x = "Mean coverage (low masking)",
    y = "Mean coverage (high masking)"
  ) +
  
  theme_aoas

g2
ggsave("Figures_cond/mean_coverage_regional_neighbors.pdf", g2)
