rm(list = ls())

source("R/forecast_function.R")
source("R/reconciliation_function.R")
source("R/metric_function.R")

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(purrr)

# =====================================================
# PARAMETERS
# =====================================================

results_root <- "Output_conformal"
alpha <- 0.1
min_obs <- 20

model_names <- list(
  "Combination",
  "Combination",
  "Combination"
)

# =====================================================
# LOAD STRUCTURAL DATA
# =====================================================

H <- readRDS("Data/structural.RDS")

res_nat <- readRDS("Output_ponctual/Seed_1/results_scotland.RDS")
res_reg <- readRDS("Output_ponctual/Seed_1/results_regions.RDS")
res_sta <- readRDS("Output_ponctual/Seed_1/results_stations.RDS")

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
# CONDITIONAL GROUPING FUNCTIONS
# =====================================================

get_neighbors <- function(i, H) {
  
  region_i <- which(H[i, 3:35] == 1)
  neighbors <- which(H[, region_i + 2] == 1)
  
  setdiff(neighbors, i)
}

extract_structure_groups <- function(i, M, H, test_period) {
  
  test_idx <- which(rownames(M) %in% as.character(test_period))
  active_idx <- test_idx[M[test_idx, 1 + 32 + i] == 0]
  
  if(length(active_idx) == 0) return(NULL)
  
  neighbors <- get_neighbors(i, H)
  
  n_active_neighbors <- sapply(
    active_idx,
    function(t) sum(M[t, 1 + 32 + neighbors] == 0)
  )
  
  q1 <- quantile(n_active_neighbors, 0.25, na.rm = TRUE)
  q2 <- quantile(n_active_neighbors, 0.5, na.rm = TRUE)
  q3 <- quantile(n_active_neighbors, 0.75, na.rm = TRUE)

  # low_idx <- active_idx[n_active_neighbors <= q1]
  # high_idx <- active_idx[n_active_neighbors >= q3]
  low_idx <- active_idx[n_active_neighbors <= q2]
  high_idx <- active_idx[n_active_neighbors >= q2]
  
  data.frame(
    node = as.character(i),
    date = as.Date(c(
      rownames(M)[low_idx],
      rownames(M)[high_idx],
      rownames(M)[active_idx]
    )),
    structure_type = c(
      rep("low", length(low_idx)),
      rep("high", length(high_idx)),
      rep("marginal", length(active_idx))
    )
  )
}


count_neightboor <- function(i, M, H, test_period) {
  
  test_idx <- which(rownames(M) %in% as.character(test_period))
  active_idx <- test_idx[M[test_idx, 1 + 32 + i] == 0]
  
  if(length(active_idx) == 0) return(NULL)
  
  neighbors <- get_neighbors(i, H)
  
  n_active_neighbors <- sapply(
    active_idx,
    function(t) sum(M[t, 1 + 32 + neighbors] == 0)
  )
  
  return(data.frame(Date = test_period[active_idx],
                       number = n_active_neighbors)
         )
}

number_active_410 <- count_neightboor(410, M, H, test_period)

plot(number_active_410$Date, number_active_410$number) # minimal infrastructure is of 45 and 46 nodes (on August, 10th and 11th).
mean(number_active_410$number) # 67

# =====================================================
# BUILD MASK DF (COMMON TO ALL SEEDS)
# =====================================================

mask_df <- bind_rows(
  lapply(seq_len(nrow(H)), function(i) {
    extract_structure_groups(i, M, H, test_period)
  })
)

# =====================================================
# DETECT SEEDS
# =====================================================

seed_dirs <- list.dirs(
  results_root,
  recursive = FALSE,
  full.names = TRUE
)

seed_dirs <- seed_dirs[
  grepl("Seed_", basename(seed_dirs))
]

# =====================================================
# HELPERS
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

# =====================================================
# LOAD ALL TIME RESULTS
# =====================================================

time_list <- list()

for (seed_dir in seed_dirs) {
  
  seed_id <- str_extract(basename(seed_dir), "\\d+")
  
  result_files <- list.files(
    seed_dir,
    pattern = "\\.RDS$",
    full.names = TRUE
  )
  
  result_files <- result_files[
    !grepl("forecast|result|number", basename(result_files))
  ]
  
  for (f in result_files) {
    
    obj <- readRDS(f)
    
    if (is.null(obj$time)) next
    
    tmp <- obj$time
    tmp$file_source <- basename(f)
    tmp$seed <- as.integer(seed_id)
    
    time_list[[length(time_list)+1]] <- tmp
  }
}

df_time <- bind_rows(time_list)

# =====================================================
# METADATA
# =====================================================

df_time <- df_time %>%
  mutate(
    conformal_method = infer_conformal_method(file_source),
    projection_method = infer_projection_method(method)
  )

# =====================================================
# CONDITIONAL COVERAGE PER NODE
# =====================================================
station_nodes <- unique(res_sta$Station.ID)
df_conditional <- df_time %>%
  filter(
    node %in% station_nodes,
    projection_method %in% c("refined OLS", "Nested")
  ) %>%
  inner_join(
    mask_df,
    by = c("node", "date")
  ) %>%
  group_by(
    seed,
    node,
    conformal_method,
    projection_method,
    structure_type
  ) %>%
  summarise(
    n_obs = n(),
    empirical_coverage = mean(covered),
    .groups = "drop"
  ) %>%
  filter(n_obs >= min_obs)

# =====================================================
# VALIDITY RATE PER SEED
# =====================================================

df_validity_seed <- df_conditional %>%
  group_by(
    seed,
    conformal_method,
    projection_method,
    structure_type
  ) %>%
  summarise(
    n_nodes = n(),
    validity_rate =
      mean(empirical_coverage >= (1 - alpha), na.rm = T),
    .groups = "drop"
  )

# =====================================================
# SUMMARY ACROSS SEEDS
# =====================================================

df_validity_summary <- df_validity_seed %>%
  group_by(
    conformal_method,
    projection_method,
    structure_type
  ) %>%
  summarise(
    mean_validity = mean(validity_rate),
    sd_validity = sd(validity_rate),
    se_validity = sd_validity / sqrt(n()),
    lower = mean_validity - 1.96 * se_validity,
    upper = mean_validity + 1.96 * se_validity,
    .groups = "drop"
  )

df_plot <- df_validity_seed %>%
  mutate(
    conformal_method = case_when(
      conformal_method == "CP_MNR" ~ "CP",
      conformal_method == "CP_MNR_naive" ~ "CP-MNR",
      conformal_method == "Adaptive_CP_MNR" ~ "Adaptive CP-MNR",
      conformal_method == "CP_MNR_Nested_star" ~ "CP-MNR-Nested*",
      TRUE ~ conformal_method
    )
  )

df_summary_plot <- df_validity_summary %>%
  mutate(
    conformal_method = case_when(
      conformal_method == "CP_MNR_naive" ~ "CP",
      conformal_method == "CP_MNR" ~ "CP-MNR",
      conformal_method == "Adaptive_CP_MNR" ~ "Adaptive CP-MNR",
      conformal_method == "CP_MNR_Nested_star" ~ "CP-MNR-Nested*",
      TRUE ~ conformal_method
    )
  )
method_levels <- c(
  "CP",
  "CP-MNR",
  "CP-MNR-Nested*",
  "Adaptive CP-MNR"
)

df_summary_plot$conformal_method <- factor(
  df_summary_plot$conformal_method,
  levels = method_levels
)

method_colors <- c(
  "CP" = "#FDB0AA",
  "CP-MNR" = "#7570B3",
  "Adaptive CP-MNR" = "#59A14F",
  "CP-MNR-Nested*" = "#E15759"
)

shape_map <- c(
  "marginal" = 23,  
  "low"      = 25,  # triangle down
  "high"     = 24   # triangle up
)

library(ggplot2)

df_summary_plot$conformal_method <- factor(
  df_summary_plot$conformal_method,
  levels = c("CP", "CP-MNR", "Adaptive CP-MNR", "CP-MNR-Nested*")
)

df_summary_plot$structure_type <- factor(
  df_summary_plot$structure_type,
  levels = c("marginal", "low", "high")
)

theme_aoas <- function(base_size = 13, base_family = "serif") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = rel(1.2),
        hjust = 0,
        margin = margin(t = 0, b = 8) 
      ),
      plot.subtitle = element_text(
        size = rel(0.95),
        color = "grey30"
      ),
      axis.title = element_text(
        face = "bold"
      ),
      axis.text = element_text(
        color = "black"
      ),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(
        color = "grey85",
        linewidth = 0.35
      ),
      axis.line = element_line(
        color = "black",
        linewidth = 0.3
      ),
      legend.position = "top",
      legend.title = element_blank(),
      plot.margin = margin(10, 15, 10, 10)
    )
}

library(ggplot2)
library(patchwork)
library(cowplot)   # for get_legend

p_main <- ggplot(
  df_summary_plot,
  aes(
    x = conformal_method,
    y = mean_validity,
    color = conformal_method,
    shape = structure_type,
    group = structure_type
  )
) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.5,
    linewidth = 0.75,
    position = position_dodge(width = 0.6)
  ) +
  geom_point(
    size = 4,
    stroke = 2,
    position = position_dodge(width = 0.6)
  ) +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = shape_map) +
  theme_minimal(base_size = 32) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none"
  ) +
  labs(y = "",) 

p_main <- p_main +
  scale_y_continuous(
    labels = label_percent(accuracy = 1)
  )
p_main
ggsave(
  filename = "Figures/conditional.pdf",
  plot = p_main,
  width = 14,
  height = 6
)

p_method_for_legend <- ggplot(
  df_summary_plot,
  aes(conformal_method, mean_validity, color = conformal_method)
) +
  geom_point(size = 4) +
  scale_color_manual(values = method_colors) +
  guides(
    color = guide_legend(
      nrow = 1,
      byrow = TRUE,
      title = NULL,
      override.aes = list(shape = 15, size = 5)
    )
  ) +
  theme_aoas() +
  theme(
    legend.position = "bottom",
    
    legend.box = "horizontal",
    legend.margin = margin(0, 0, 0, 0),
    legend.spacing.x = unit(0.4, "cm"),
    legend.key.width = unit(1.2, "cm"),
    legend.text = element_text(
      size = 18,
      margin = margin(l = 1)  
    ),
    plot.margin = margin(0, 0, 0, 0)
  )

p_struct_for_legend <- ggplot(
  df_summary_plot,
  aes(structure_type, mean_validity, shape = structure_type)
) +
  geom_point(
    size = 4,
    stroke = 1.5,
    color = "black"
  ) +
scale_shape_manual(values = shape_map,
                   labels = c(
                     "high" = "High activity",
                     "marginal"  = "Marginal",
                     "low"  = "Low activity"
                   )) +
  guides(
    shape = guide_legend(title = NULL,
                         byrow = TRUE)
  ) +
  theme_aoas() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.margin = margin(0, 0, 0, 0),
    legend.text = element_text(size = 18),
    legend.spacing.y= unit(2.5, "cm"),
    legend.key.height = unit(1, "cm"),
    plot.margin = margin(0, 0, 0, 0)
  )

leg_method <- cowplot::get_legend(p_method_for_legend)
leg_struct <- cowplot::get_legend(p_struct_for_legend)

ggsave(
  filename = "Figures/conditional_legend_bottom.pdf",
  plot = leg_method,
  width = 14,
  height = 2
)

ggsave(
  filename = "Figures/conditional_legend_right.pdf",
  plot = leg_struct,
  width = 2,
  height = 12
)