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

results_root <- "Output_conformal"
figures_dir <- "Figures"

alpha <- 0.1

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

# =====================================================
# LOAD ALL RESULTS ACROSS SEEDS
# =====================================================

summary_list <- list()
time_list <- list()

for (seed_dir in seed_dirs) {
  
  seed_id <- str_extract(basename(seed_dir), "\\d+")
  
  result_files <- list.files(
    seed_dir,
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
  
  cat("\nSeed", seed_id, ":\n")
  print(basename(result_files))
  
  for (f in result_files) {
    
    obj <- readRDS(f)
    
    if (!is.list(obj)) next
    
    if (!is.null(obj$summary)) {
      
      tmp <- obj$summary
      tmp$file_source <- basename(f)
      tmp$seed <- as.integer(seed_id)
      
      summary_list[[length(summary_list) + 1]] <- tmp
    }
    
    if (!is.null(obj$time)) {
      
      tmp <- obj$time
      tmp$file_source <- basename(f)
      tmp$seed <- as.integer(seed_id)
      
      time_list[[length(time_list) + 1]] <- tmp
    }
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
infer_base_model <- function(x) {
  x %>%
    str_remove("\\.RDS$") %>%
    str_split("__") %>%
    map_chr(~ paste(tail(.x, 3), collapse = "__"))
}

df_summary <- df_summary %>%
  mutate(
    conformal_method = infer_conformal_method(file_source),
    projection_method = infer_projection_method(method),
    base_model = infer_base_model(file_source),
    experiment = paste(
      conformal_method,
      projection_method,
      base_model,
      sep = " | "
    )
  )

# =====================================================
# VALIDITY RATE METRICS
# =====================================================

min_obs <- 20   # minimal number of observation

# ---------------------------------
# Step 1: empirical coverage by node
# ---------------------------------

df_node_coverage <- df_summary %>%
  group_by(
    seed,
    experiment,
    node_type,
    node,
    conformal_method,
    projection_method,
    base_model
  ) %>%
  summarise(
    total_eff = sum(n_eff, na.rm = TRUE),
    empirical_coverage = weighted.mean(
      coverage,
      w = n_eff,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# ---------------------------------
# Step 2: filter low-support nodes
# ---------------------------------

df_node_coverage <- df_node_coverage %>%
  filter(total_eff >= min_obs)

# ---------------------------------
# Step 3: validity rates per seed
# ---------------------------------

df_validity_seed <- df_node_coverage %>%
  group_by(
    seed,
    experiment,
    conformal_method,
    projection_method,
    base_model,
    node_type
  ) %>%
  summarise(
    n_nodes = n(),
    
    strict_validity_rate =
      mean(empirical_coverage >= (1 - alpha)),
    
    relaxed_validity_rate =
      mean(empirical_coverage >= (1 - 2 * alpha)),
    
    .groups = "drop"
  )

# ---------------------------------
# Step 4: aggregate across seeds
# ---------------------------------

df_validity_summary <- df_validity_seed %>%
  group_by(
    experiment,
    conformal_method,
    projection_method,
    base_model,
    node_type
  ) %>%
  summarise(
    n_seeds = n(),
    
    strict_mean = mean(strict_validity_rate),
    strict_sd = sd(strict_validity_rate),
    strict_se = strict_sd / sqrt(n_seeds),
    strict_lower = strict_mean - 1.96 * strict_se,
    strict_upper = strict_mean + 1.96 * strict_se,
    
    relaxed_mean = mean(relaxed_validity_rate),
    relaxed_sd = sd(relaxed_validity_rate),
    relaxed_se = relaxed_sd / sqrt(n_seeds),
    relaxed_lower = relaxed_mean - 1.96 * relaxed_se,
    relaxed_upper = relaxed_mean + 1.96 * relaxed_se,
    
    .groups = "drop"
  )

library(stringr)
library(knitr)
library(kableExtra)

df_validity_table <- df_validity_summary %>%
  filter(
    projection_method %in% c("refined OLS", "Nested")
  ) %>%
  mutate(
    method = case_when(
      conformal_method == "CP_MNR_naive" ~ "CP",
      conformal_method == "CP_MNR" ~ "CP-MNR",
      conformal_method == "CP_MNR_Nested_star" ~ "CP-MNR-Nested$^{\\star}$",
      conformal_method == "Adaptive_CP_MNR" ~ "Adaptive CP-MNR",
      TRUE ~ conformal_method
    )
  ) %>%
  mutate(
    method = factor(
      method,
      levels = c(
        "CP",
        "CP-MNR",
        "CP-MNR-Nested$^{\\star}$",
        "Adaptive CP-MNR"
      )
    )
  )

# =====================================================
# FLAG BEST METHODS
# =====================================================

df_validity_table <- df_validity_table %>%
  group_by(node_type) %>%
  mutate(
    best_strict = strict_mean == max(strict_mean, na.rm = TRUE),
    best_relaxed = relaxed_mean == max(relaxed_mean, na.rm = TRUE)
  ) %>%
  ungroup()

# =====================================================
# FORMAT CELLS
# =====================================================

df_validity_table <- df_validity_table %>%
  mutate(
    strict_display = sprintf(
      "%.1f \\%% $\\pm$ %.1f \\%%",
      100 * strict_mean,
      100 * strict_se
    ),
    relaxed_display = sprintf(
      "%.1f \\%% $\\pm$ %.1f \\%%",
      100 * relaxed_mean,
      100 * relaxed_se
    ),
    
    strict_display = ifelse(
      best_strict,
      paste0("\\textbf{", strict_display, "}"),
      strict_display
    ),
    
    relaxed_display = ifelse(
      best_relaxed,
      paste0("\\textbf{", relaxed_display, "}"),
      relaxed_display
    )
  )

# =====================================================
# WIDE FORMAT
# =====================================================

table_strict <- df_validity_table %>%
  dplyr::select(
    node_type,
    conformal_method,
    strict_display
  ) %>%
  pivot_wider(
    names_from = conformal_method,
    values_from = strict_display,
    names_glue = "{conformal_method}__strict"
  )

table_relaxed <- df_validity_table %>%
  dplyr::select(
    node_type,
    conformal_method,
    relaxed_display
  ) %>%
  pivot_wider(
    names_from = conformal_method,
    values_from = relaxed_display,
    names_glue = "{conformal_method}__relaxed"
  )

# Merge
latex_table <- table_strict %>%
  left_join(
    table_relaxed,
    by = "node_type"
  )

rank_df <- df_validity_table %>%
  group_by(node_type) %>%
  mutate(
    rank_strict = rank(-strict_mean, ties.method = "average"),
    rank_relaxed = rank(-relaxed_mean, ties.method = "average")
  ) %>%
  ungroup()
rank_by_node <- rank_df %>%
  group_by(node_type) %>%
  summarise(
    CP_strict = mean(rank_strict[method == "CP"]),
    CP_MNR_strict = mean(rank_strict[method == "CP-MNR"]),
    CP_MNR_Nested_strict = mean(rank_strict[method == "CP-MNR-Nested$^{\\star}$"]),
    Adaptive_strict = mean(rank_strict[method == "Adaptive CP-MNR"]),
    
    CP_relaxed = mean(rank_relaxed[method == "CP"]),
    CP_MNR_relaxed = mean(rank_relaxed[method == "CP-MNR"]),
    CP_MNR_Nested_relaxed = mean(rank_relaxed[method == "CP-MNR-Nested$^{\\star}$"]),
    Adaptive_relaxed = mean(rank_relaxed[method == "Adaptive CP-MNR"]),
    .groups = "drop"
  )
mean_rank_row <- rank_by_node %>%
  summarise(
    node_type = "Mean Rank",
    
    CP_strict = round(mean(CP_strict),2),
    CP_MNR_strict = round(mean(CP_MNR_strict),2),
    CP_MNR_Nested_strict = round(mean(CP_MNR_Nested_strict),2),
    Adaptive_strict = round(mean(Adaptive_strict),2),
    
    CP_relaxed = round(mean(CP_relaxed),2),
    CP_MNR_relaxed = round(mean(CP_MNR_relaxed,2)),
    CP_MNR_Nested_relaxed = round(mean(CP_MNR_Nested_relaxed),2),
    Adaptive_relaxed = round(mean(Adaptive_relaxed),2)
  )
names(mean_rank_row) <- names(latex_table)
latex_table <- latex_table %>% rbind(mean_rank_row)

# Reorder columns
latex_table <- latex_table %>%
  dplyr::select(
    node_type,
    
    CP_MNR_naive__strict,
    CP_MNR_naive__relaxed,
    
    CP_MNR__strict,
    CP_MNR__relaxed,
    
    CP_MNR_Nested_star__strict,
    CP_MNR_Nested_star__relaxed,
    
    Adaptive_CP_MNR__strict,
    Adaptive_CP_MNR__relaxed
  )

# =====================================================
# LATEX OUTPUT
# =====================================================

library(knitr)
library(kableExtra)

latex_output <- kbl(
  latex_table,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  align = "l|cc|cc|cc|cc",
  col.names = c(
    "Node type",
    
    "$1-\\alpha$",
    "$1-2\\alpha$",
    
    "$1-\\alpha$",
    "$1-2\\alpha$",
    
    "$1-\\alpha$",
    "$1-2\\alpha$",
    
    "$1-\\alpha$",
    "$1-2\\alpha$"
  )
) %>%
  add_header_above(
    c(
      " " = 1,
      "CP" = 2,
      "CP-MNR" = 2,
      "CP-MNR Nested" = 2,
      "Adaptive CP-MNR" = 2
    )
  )
cat(latex_output)


df_validity_table <- df_validity_summary %>%
  filter(
    projection_method %in% c("OLS", "refined OLS", "Direct")
  ) %>%
  mutate(
    projection = recode(
      projection_method,
      "OLS" = "OLS",
      "refined OLS" = "rOLS",
      "Direct" = "Direct"
    ),
    method = recode(
      conformal_method,
      "CP_MNR_naive" = "CP",
      "CP_MNR" = "CP-MNR",
      "Adaptive_CP_MNR" = "Adaptive CP-MNR"
    )
  ) %>%
  mutate(
    method = factor(
      method,
      levels = c(
        "CP",
        "CP-MNR",
        "Adaptive CP-MNR"
      )
    )
  )

df_validity_table <- df_validity_table %>%
  group_by(node_type, method) %>%
  mutate(
    best_strict  = strict_mean  == max(strict_mean,  na.rm = TRUE),
    best_relaxed = relaxed_mean == max(relaxed_mean, na.rm = TRUE)
  ) %>%
  ungroup()

df_validity_table <- df_validity_table %>%
  mutate(
    
    strict_display = sprintf(
      "%.1f \\%% $\\pm$ %.1f \\%%",
      100 * strict_mean,
      100 * strict_se
    ),
    
    relaxed_display = sprintf(
      "%.1f \\%% $\\pm$ %.1f \\%%",
      100 * relaxed_mean,
      100 * relaxed_se
    ),
    
    strict_display = ifelse(
      best_strict,
      paste0("\\textbf{", strict_display, "}"),
      strict_display
    ),
    
    relaxed_display = ifelse(
      best_relaxed,
      paste0("\\textbf{", relaxed_display, "}"),
      relaxed_display
    )
  )

make_validity_table <- function(df, value){
  
  tab <- df %>%
    dplyr::select(
      node_type,
      projection,
      method,
      {{ value }}
    ) %>%
    pivot_wider(
      names_from = method,
      values_from = {{ value }}
    ) %>%
    dplyr::select(
      node_type,
      projection,
      CP,
      `CP-MNR`,
      `Adaptive CP-MNR`
    ) %>%
    arrange(
      factor(node_type,
             levels = c("National","Regional","Station")),
      factor(projection,
             levels = c("Direct","OLS","rOLS"))
    )
  
  ## First row of each node type
  tab <- tab %>%
    group_by(node_type) %>%
    mutate(
      node_type = ifelse(row_number()==1,node_type,"")
    ) %>%
    ungroup()
  
  kbl(
    tab,
    format = "latex",
    escape = FALSE,
    booktabs = TRUE,
    align = "llccc",
    col.names = c(
      "Node type",
      "Reconciliation",
      "CP",
      "CP-MNR",
      "Adaptive CP-MNR"
    )
  ) %>%
    kable_styling(latex_options = "hold_position")
}

latex_strict <- make_validity_table(
    df_validity_table,
    strict_display
  )

latex_relaxed <-
  make_validity_table(
    df_validity_table,
    relaxed_display
  )

cat(latex_strict)

cat(latex_relaxed)


# =====================================================
# LENGTH METRIC
# =====================================================

df_node_length <- df_summary %>%
  group_by(
    seed,
    experiment,
    node_type,
    node,
    conformal_method,
    projection_method,
    base_model,
    method
  ) %>%
  summarise(
    total_eff = sum(n_eff, na.rm = TRUE),
    length_mean = weighted.mean(
      length,
      w = n_eff,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

df_node_length <- df_node_length %>%
  filter(total_eff >= min_obs) %>%
  filter(method != "Nested")

df_node_length <- df_node_length %>%
  group_by(
    seed,
    node,
    conformal_method,
    base_model
  ) %>%
  mutate(
    direct_length = length_mean[method == "Direct"]
  ) %>%
  ungroup()

df_node_length <- df_node_length %>%
  mutate(
    win_vs_direct = length_mean < direct_length
  )

df_efficiency_seed <- df_node_length %>%
  filter(method != "Direct") %>%
  group_by(
    seed,
    conformal_method,
    projection_method,
    base_model,
    node_type
  ) %>%
  summarise(
    n_nodes = n(),
    
    efficiency_win_rate = mean(win_vs_direct, na.rm = TRUE),
    
    .groups = "drop"
  )

df_efficiency_summary <- df_efficiency_seed %>%
  group_by(
    conformal_method,
    projection_method,
    base_model,
    node_type
  ) %>%
  summarise(
    n_seeds = n(),
    
    efficiency_mean = mean(efficiency_win_rate),
    efficiency_sd = sd(efficiency_win_rate),
    efficiency_se = efficiency_sd / sqrt(n_seeds),
    
    efficiency_lower = efficiency_mean - 1.96 * efficiency_se,
    efficiency_upper = efficiency_mean + 1.96 * efficiency_se,
    
    .groups = "drop"
  )

# =====================================================
# EFFICIENCY TABLE (FROM df_efficiency_summary ONLY)
# =====================================================

library(dplyr)
library(tidyr)
library(stringr)
library(knitr)
library(kableExtra)

# -----------------------------------------------------
# Step 1: restrict projection methods
# -----------------------------------------------------

df_efficiency_table <- df_efficiency_summary %>%
  filter(
    projection_method %in% c("OLS", "refined OLS")
  ) %>%
  mutate(
    method = case_when(
      conformal_method == "CP_MNR" ~ "CP",
      conformal_method == "CP_MNR_naive" ~ "CP-MNR",
      conformal_method == "CP_MNR_Nested_star" ~ "CP-MNR-Nested$^{\\star}$",
      conformal_method == "Adaptive_CP_MNR" ~ "Adaptive CP-MNR",
      TRUE ~ conformal_method
    )
  ) %>%
  mutate(
    method = factor(
      method,
      levels = c(
        "CP",
        "CP-MNR",
        "CP-MNR-Nested$^{\\star}$",
        "Adaptive CP-MNR"
      )
    ),
    proj = case_when(
      projection_method == "OLS" ~ "OLS",
      projection_method == "refined OLS" ~ "rOLS",
      TRUE ~ projection_method
    )
  )

# -----------------------------------------------------
# Step 2: flag best method per node_type & projection
# -----------------------------------------------------

df_efficiency_table <- df_efficiency_table %>%
  group_by(node_type, conformal_method) %>%
  mutate(
    best_eff = efficiency_mean == max(efficiency_mean, na.rm = TRUE)
  ) %>%
  ungroup()

# -----------------------------------------------------
# Step 3: formatting (same style as validity table)
# -----------------------------------------------------

df_efficiency_table <- df_efficiency_table %>%
  mutate(
    eff_display = sprintf(
      "%.1f \\%% $\\pm$ %.1f \\%%",
      100 * efficiency_mean,
      100 * efficiency_se
    ),
    eff_display = case_when(
      efficiency_mean >= 0.5 & best_eff ~ paste0("\\textbf{", eff_display, "$^{\\star}$}"),
      efficiency_mean >= 0.5 ~ paste0(eff_display, "$^{\\star}$"),
      best_eff ~ paste0("\\textbf{", eff_display, "}"),
      TRUE ~ eff_display
    )
  )

# =====================================================
# RANKING: OLS vs rOLS
# =====================================================

df_rank <- df_efficiency_table %>%
  group_by(node_type, conformal_method) %>%
  mutate(
    rank_proj = rank(-efficiency_mean, ties.method = "average")
  ) %>%
  ungroup()

# =====================================================
# RANK OLS vs rOLS WITHIN (node_type, conformal_method)
# =====================================================

rank_df <- df_efficiency_table %>%
  group_by(node_type, conformal_method) %>%
  mutate(
    rank_eff = rank(
      -efficiency_mean,
      ties.method = "average"
    )
  ) %>%
  ungroup()

# =====================================================
# MEAN RANK BY CONFORMAL METHOD
# =====================================================

mean_rank_rows <- rank_df %>%
  group_by(conformal_method, proj) %>%
  summarise(
    efficiency_mean = mean(rank_eff),
    .groups = "drop"
  ) %>%
  mutate(
    node_type = "Mean Rank",
    n_seeds = NA,
    efficiency_sd = NA,
    efficiency_se = NA,
    efficiency_lower = NA,
    efficiency_upper = NA,
    base_model = NA,
    best_eff = FALSE,
    eff_display = sprintf("%.2f", efficiency_mean)
  )

# =====================================================
# APPEND TO MAIN TABLE
# =====================================================

# df_efficiency_table <- bind_rows(
#   df_efficiency_table,
#   mean_rank_rows
# )

# -----------------------------------------------------
# Step 4: wide format (OLS / rOLS like coverage levels)
# -----------------------------------------------------

table_efficiency <- df_efficiency_table %>%
  dplyr::select(
    node_type,
    conformal_method,
    proj,
    eff_display
  ) %>%
  tidyr::pivot_wider(
    names_from = c(conformal_method, proj),
    values_from = eff_display
  )


# Reorder columns
table_efficiency <- table_efficiency %>%
  dplyr::select(
    node_type,
    
    CP_MNR_naive_OLS,
    CP_MNR_naive_rOLS,
    
    CP_MNR_OLS,
    CP_MNR_rOLS,
    
    Adaptive_CP_MNR_OLS,
    Adaptive_CP_MNR_rOLS
  )



# -----------------------------------------------------
# Step 5: LaTeX output (same structure as validity table)
# -----------------------------------------------------

latex_efficiency <- kbl(
  table_efficiency,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  align = "l|cc|cc|cc",
  col.names = c(
    "Node type",
    "OLS",
    "rOLS",
    "OLS",
    "rOLS",
    "OLS",
    "rOLS"
  )
) %>%
  add_header_above(
    c(
      " " = 1,
      "CP" = 2,
      "CP-MNR" = 2,
      "Adaptive CP-MNR" = 2
    )
  )

cat(latex_efficiency)
