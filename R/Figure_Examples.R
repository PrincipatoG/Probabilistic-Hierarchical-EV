rm(list = ls())

source("R/forecast_function.R")
source("R/reconciliation_function.R")
source("R/metric_function.R")
source('Model_training/R/pretraitements.R')

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

df_time <- df_time %>%
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


# method_colors <- c(
#   "CP" = "#FDB0AA",
#   "CP-MNR" = "#1B9E77",
#   "CP-MNR-Nested*" = "#D95F02",
#   "Adaptive CP-MNR" = "#7570B3"
# )
method_colors <- c(
  "CP" = "#FDB0AA",
  "CP-MNR" = "#7570B3",
  "Adaptive CP-MNR" = "#59A14F",
  "CP-MNR-Nested*" = "#E15759"
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

# =====================================================
# NATIONAL NODE
# =====================================================
df_national <- df_time %>% filter(node == "Y_nat")
national_data_path <- "Data/dataset_scotland_main.csv"

my_windows <- generate_rolling_windows(SEED = 1)  
test_windows <- bind_rows(my_windows) %>% filter(type =="test")
length(unique(df_national$date))

dataset_scotland <- read.csv(national_data_path) %>% 
  mutate(date = as.Date(Date)) %>% preparation_data_scotland() %>%
  filter(date %in% test_windows$Date)

# =====================================================
# PREDICTION INTERVALS
# =====================================================
df_national <- df_national %>%
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

df_interval <- df_national %>%
  filter(
    projection_method %in% c("refined OLS", "Nested"),
    !is.na(lower),
    !is.na(upper)
  ) %>%
  
  mutate(
    conformal_method = recode(
      conformal_method,
      "CP" = "CP",
      "CP_MNR" = "CP-MNR",
      "CP_MNR_Nested_star" = "CP-MNR-Nested*",
      "Adaptive_CP_MNR" = "Adaptive CP-MNR"
    )
  ) %>%
  filter(seed==7) %>% # Just consider the intervals for the first seed
  # Aggregate by date and method (across all nodes)
  group_by(conformal_method, date) %>%
  summarise(
    lower_mean = mean(lower),
    upper_mean = mean(upper),
    center     = mean((lower + upper) / 2),
    length_mean = mean(upper - lower),
    .groups = "drop"
  ) 
df_interval <- df_interval %>%
  left_join(dataset_scotland, by = join_by(date)) %>%
  filter(!is.na(Consumed_kWh)) %>%
  mutate(covered = (Consumed_kWh >= lower_mean & Consumed_kWh <= upper_mean),
)

y_limits <- range(
  c(df_interval$lower_mean,
    df_interval$upper_mean,
    df_interval$Consumed_kWh),
  na.rm = TRUE
)
library(ggplot2)

theme_clean <- theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.4),
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )

plot_method <- function(method_name) {
  
  df_m <- df_interval %>%
    filter(conformal_method == method_name)
  
  method_col <- method_colors[method_name]
  
  ggplot(df_m, aes(x = date)) +
    
    # intervalle principal
    geom_ribbon(
      aes(ymin = lower_mean, ymax = upper_mean),
      fill = method_col,
      alpha = 0.30
    ) +
    # bornes visibles
    geom_line(
      aes(y = lower_mean),
      color = method_col,
      linewidth = 0.6
    ) +
    geom_line(
      aes(y = upper_mean),
      color = method_col,
      linewidth = 0.6
    ) +
    # Observed value
    geom_line(
      aes(y = Consumed_kWh),
      color = "black",
      linewidth = 0.5,
      alpha = 0.4
    ) +
    geom_point(
      aes(y = Consumed_kWh),
      color = "black",
      size = 0.5,
      alpha = 0.5
    ) +
    coord_cartesian(ylim = y_limits) +
    labs(
      title = method_name,
      x = "Date",
      y = "Consumed kWh"
    ) +
    theme_aoas()
  }

p2 <- plot_method("CP-MNR")
p3 <- plot_method("CP-MNR-Nested*")
p4 <- plot_method("Adaptive CP-MNR")

ggsave("Figures/interval_CP_MNR.pdf", p2, width = 4, height = 3)
ggsave("Figures/interval_CP_MNR_Nested.pdf", p3, width = 4, height = 3)
ggsave("Figures/interval_Adaptive_CP_MNR.pdf", p4, width = 4, height = 3)

p2
p3
p4


# =====================================================
# TIME AVERAGE
# =====================================================
df_national <- df_national %>%
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

df_summary_nat <- df_national %>%
  filter(
    projection_method %in% c("refined OLS", "Nested"),
    !is.na(covered),
    !is.na(length)
  ) %>%
  
  mutate(
    conformal_method = recode(
      conformal_method,
      "CP_MNR_naive" = "CP",
      "CP_MNR" = "CP-MNR",
      "CP_MNR_Nested_star" = "CP-MNR-Nested*",
      "Adaptive_CP_MNR" = "Adaptive CP-MNR"
    )
  ) %>% mutate(
    conformal_method = factor(
      conformal_method,
      levels = c(
        "CP",
        "CP-MNR",
        "Adaptive CP-MNR",
        "CP-MNR-Nested*"
      )
    )
  ) %>% group_by(conformal_method, seed) %>%
  summarise(
    coverage = mean(covered),
    avg_length = mean(length),
    .groups = "drop"
  ) %>%
    group_by(conformal_method) %>%
  summarise(
    mean_coverage = mean(coverage),
    se_coverage   = sd(coverage) / sqrt(n()),
    mean_length   = mean(avg_length),
    se_length     = sd(avg_length) / sqrt(n()),
    n_seeds = n(),
    .groups = "drop"
  ) %>%
  filter(conformal_method != "CP")

library(ggplot2)

theme_clean <- theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.4),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold")
  )
library(ggplot2)
library(cowplot)

g1 <- ggplot(
  df_summary_nat,
  aes(
    x = mean_coverage,
    y = mean_length,
    color = conformal_method,
    shape = conformal_method
  )
) +
  geom_errorbarh(
    aes(
      xmin = mean_coverage - 1.96 * se_coverage,
      xmax = mean_coverage + 1.96 * se_coverage
    ),
    height = 250,
    linewidth = 0.5
  ) +
  geom_errorbar(
    aes(
      ymin = mean_length - 1.96 * se_length,
      ymax = mean_length + 1.96 * se_length
    ),
    width = 0.001,
    linewidth = 0.5
  ) +
  geom_point(size = 3, stroke = 1.2) +
  scale_shape_manual(values = c(16, 17, 15)) +
  scale_color_manual(values = method_colors) +
  scale_y_continuous(
    breaks = function(x) {
      rng <- range(x, na.rm = TRUE)
      pretty(seq(rng[1], rng[2], length.out = 100), n = 3)
    },
    minor_breaks = NULL
  ) +
  geom_vline(
    xintercept = 1 - alpha,
    linetype = "dashed",
    linewidth = 0.7,
    color = "grey"
  ) +
  labs(
    x = "Coverage",
    y = "Length",
    title = "Scotland"
  ) +
  coord_cartesian(xlim = c(0.885, 0.915)) +
  theme_aoas() +
  theme(
    legend.position = "none"
    )
g1

ggsave("Figures/CovLen_Scotland.pdf", g1, width = 4, height = 2.5)

# =====================================================
# GLASGOW CITY NODE
# =====================================================
df_glasgow <- df_time %>% filter(node == "Glasgow City")
regional_data_path <- "Data/dataset_region_main.csv"
dataset_glasgow <- read.csv(regional_data_path) %>% 
  mutate(date = as.Date(Date)) %>% preparation_data_regions_stations(group_col=Region) %>%
  filter(Region == "Glasgow City") %>%
  filter(date %in% test_windows$Date)

df_glasgow <- df_glasgow %>%
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
# PREDICTION INTERVALS
# =====================================================
df_interval <- df_glasgow %>%
  filter(
    projection_method %in% c("refined OLS", "Nested"),
    !is.na(lower),
    !is.na(upper)
  ) %>%
  
  mutate(
    conformal_method = recode(
      conformal_method,
      "CP" = "CP",
      "CP_MNR" = "CP-MNR",
      "CP_MNR_Nested_star" = "CP-MNR-Nested*",
      "Adaptive_CP_MNR" = "Adaptive CP-MNR"
    )
  ) %>%
  filter(seed==7) %>% # Just consider the intervals for the first seed
  # Aggregate by date and method (across all nodes)
  group_by(conformal_method, date) %>%
  summarise(
    lower_mean = mean(lower),
    upper_mean = mean(upper),
    center     = mean((lower + upper) / 2),
    length_mean = mean(upper - lower),
    .groups = "drop"
  ) 
df_interval <- df_interval %>%
  left_join(dataset_glasgow, by = join_by(date)) %>%
  filter(!is.na(Consumed_kWh)) %>%
  mutate(covered = (Consumed_kWh >= lower_mean & Consumed_kWh <= upper_mean),
  )

y_limits <- range(
  c(df_interval$lower_mean,
    df_interval$upper_mean,
    df_interval$Consumed_kWh),
  na.rm = TRUE
)
library(ggplot2)

theme_clean <- theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.4),
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )

plot_method <- function(method_name) {
  
  df_m <- df_interval %>%
    filter(conformal_method == method_name)
  
  method_col <- method_colors[method_name]
  
  ggplot(df_m, aes(x = date)) +
    
    # intervalle principal
    geom_ribbon(
      aes(ymin = lower_mean, ymax = upper_mean),
      fill = method_col,
      alpha = 0.30
    ) +
    # bornes visibles
    geom_line(
      aes(y = lower_mean),
      color = method_col,
      linewidth = 0.6
    ) +
    geom_line(
      aes(y = upper_mean),
      color = method_col,
      linewidth = 0.6
    ) +
    # Observed value
    geom_line(
      aes(y = Consumed_kWh),
      color = "black",
      linewidth = 0.5,
      alpha = 0.4
    ) +
    geom_point(
      aes(y = Consumed_kWh),
      color = "black",
      size = 0.5,
      alpha = 0.5
    ) +
    coord_cartesian(ylim = y_limits) +
    labs(
      title = method_name,
      x = "Date",
      y = "Consumed kWh"
    ) +
    theme_aoas()
}

p2 <- plot_method("CP-MNR")
p3 <- plot_method("CP-MNR-Nested*")
p4 <- plot_method("Adaptive CP-MNR")

ggsave("Figures/interval_Glasgow_CP_MNR.pdf", p2, width = 4, height = 3)
ggsave("Figures/interval_Glasgow_CP_MNR_Nested.pdf", p3, width = 4, height = 3)
ggsave("Figures/interval_Glasgow_Adaptive_CP_MNR.pdf", p4, width = 4, height = 3)

p2
p3
p4


# =====================================================
# TIME AVERAGE
# =====================================================

df_summary_glasgow <- df_glasgow %>%
  filter(
    projection_method %in% c("refined OLS", "Nested"),
    !is.na(covered),
    !is.na(length)
  ) %>%
  
  mutate(
    conformal_method = recode(
      conformal_method,
      "CP_MNR_naive" = "CP",
      "CP_MNR" = "CP-MNR",
      "CP_MNR_Nested_star" = "CP-MNR-Nested*",
      "Adaptive_CP_MNR" = "Adaptive CP-MNR"
    )
  ) %>% mutate(
    conformal_method = factor(
      conformal_method,
      levels = c(
        "CP",
        "CP-MNR",
        "Adaptive CP-MNR",
        "CP-MNR-Nested*"
      )
    )
  ) %>% group_by(conformal_method, seed) %>%
  summarise(
    coverage = mean(covered),
    avg_length = mean(length),
    .groups = "drop"
  ) %>%
  group_by(conformal_method) %>%
  summarise(
    mean_coverage = mean(coverage),
    se_coverage   = sd(coverage) / sqrt(n()),
    mean_length   = mean(avg_length),
    se_length     = sd(avg_length) / sqrt(n()),
    n_seeds = n(),
    .groups = "drop"
  ) %>%
  filter(conformal_method != "CP")

library(ggplot2)

theme_clean <- theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.4),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold")
  )

g2 <- ggplot(
  df_summary_glasgow,
  aes(
    x = mean_coverage,
    y = mean_length,
    color = conformal_method,
    shape = conformal_method
  )
) +
  
  # horizontal SE
  geom_errorbarh(
    aes(
      xmin = mean_coverage - 1.96 * se_coverage,
      xmax = mean_coverage + 1.96 * se_coverage
    ),
    height = 45,
    linewidth = 0.5
  ) +
  
  # vertical SE
  geom_errorbar(
    aes(
      ymin = mean_length - 1.96 * se_length,
      ymax = mean_length + 1.96 * se_length
    ),
    width = 0.004,
    linewidth = 0.5
  ) +
  
  # points
  geom_point(size = 3, stroke = 1.2) +
  
  # shapes distinctes et pleines
  scale_shape_manual(
    values = c(16, 17, 15)
  ) +
  scale_color_manual(values = method_colors) +
  
  # target coverage
  geom_vline(
    xintercept = 1 - alpha,
    linetype = "dashed",
    linewidth = 0.7,
    color = "grey"
  ) +
  
  labs(
    x = "Coverage",
    y = "Length",
    title = "Glasgow"
  ) +
  
  coord_cartesian(xlim = c(0.825, 0.975)) +
  theme_aoas() +
  theme(
    legend.position = "none"
  )
g2

ggsave("Figures/CovLen_Glasgow.pdf", g2, width = 4, height = 2.5)

# =====================================================
# STATION 640 NODE
# =====================================================
df_station <- df_time %>% filter(node == 640)
station_data_path <- "Data/dataset_address_main.csv"
dataset_station <- read.csv(station_data_path) %>% 
  mutate(date = as.Date(Date)) %>% preparation_data_regions_stations(group_col=Station.ID) %>%
  filter(Station.ID == 640) %>%
  filter(date %in% test_windows$Date)

df_station <- df_station %>%
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
# PREDICTION INTERVALS
# =====================================================
df_interval <- df_station %>%
  filter(
    projection_method %in% c("refined OLS", "Nested"),
    !is.na(lower),
    !is.na(upper)
  ) %>%
  
  mutate(
    conformal_method = recode(
      conformal_method,
      "CP_MNR_naive" = "CP",
      "CP_MNR" = "CP-MNR",
      "CP_MNR_Nested_star" = "CP-MNR-Nested*",
      "Adaptive_CP_MNR" = "Adaptive CP-MNR"
    )
  ) %>%
  filter(seed==7) %>% # Just consider the intervals for one arbitrary seed (the choice does not matter)
  # aggregation by date and method (over all nodes) 
  group_by(conformal_method, date) %>%
  summarise(
    lower_mean = mean(lower),
    upper_mean = mean(upper),
    center     = mean((lower + upper) / 2),
    length_mean = mean(upper - lower),
    .groups = "drop"
  ) 
df_interval <- df_interval %>%
  left_join(dataset_station, by = join_by(date)) %>%
  filter(!is.na(Consumed_kWh)) %>%
  mutate(covered = (Consumed_kWh >= lower_mean & Consumed_kWh <= upper_mean),
  )

y_limits <- range(
  c(df_interval$lower_mean,
    df_interval$upper_mean,
    df_interval$Consumed_kWh),
  na.rm = TRUE
)
library(ggplot2)

theme_clean <- theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.4),
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )

plot_method <- function(method_name) {
  
  df_m <- df_interval %>%
    filter(conformal_method == method_name)
  
  method_col <- method_colors[method_name]
  
  ggplot(df_m, aes(x = date)) +
    
    # main interval
    geom_ribbon(
      aes(ymin = lower_mean, ymax = upper_mean),
      fill = method_col,
      alpha = 0.30
    ) +
    # visible bounds
    geom_line(
      aes(y = lower_mean),
      color = method_col,
      linewidth = 0.6
    ) +
    geom_line(
      aes(y = upper_mean),
      color = method_col,
      linewidth = 0.6
    ) +
    # Observed value
    geom_line(
      aes(y = Consumed_kWh),
      color = "black",
      linewidth = 0.5,
      alpha = 0.4
    ) +
    geom_point(
      aes(y = Consumed_kWh),
      color = "black",
      size = 0.5,
      alpha = 0.5
    ) +
    coord_cartesian(ylim = y_limits) +
    labs(
      title = method_name,
      x = "Date",
      y = "Consumed kWh"
    ) +
    theme_aoas()
}
p1 <- plot_method("CP")
p2 <- plot_method("CP-MNR")
p3 <- plot_method("CP-MNR-Nested*")
p4 <- plot_method("Adaptive CP-MNR")

ggsave("Figures/interval_Station_CP.pdf", p1, width = 4, height = 3)
ggsave("Figures/interval_Station_CP_MNR.pdf", p2, width = 4, height = 3)
ggsave("Figures/interval_Station_CP_MNR_Nested.pdf", p3, width = 4, height = 3)
ggsave("Figures/interval_Station_Adaptive_CP_MNR.pdf", p4, width = 4, height = 3)

p1
p2
p3
p4

# =====================================================
# TIME AVERAGE
# =====================================================

df_summary_station <- df_station %>%
  filter(
    projection_method %in% c("refined OLS", "Nested"),
    !is.na(covered),
    !is.na(length)
  ) %>%
  
  mutate(
    conformal_method = recode(
      conformal_method,
      "CP_MNR_naive" = "CP",
      "CP_MNR" = "CP-MNR",
      "CP_MNR_Nested_star" = "CP-MNR-Nested*",
      "Adaptive_CP_MNR" = "Adaptive CP-MNR"
    )
  ) %>% mutate(
    conformal_method = factor(
      conformal_method,
      levels = c(
        "CP",
        "CP-MNR",
        "Adaptive CP-MNR",
        "CP-MNR-Nested*"
      )
    )
  ) %>% group_by(conformal_method, seed) %>%
  summarise(
    coverage = mean(covered),
    avg_length = mean(length),
    .groups = "drop"
  ) %>%
  group_by(conformal_method) %>%
  summarise(
    mean_coverage = mean(coverage),
    se_coverage   = sd(coverage) / sqrt(n()),
    mean_length   = mean(avg_length),
    se_length     = sd(avg_length) / sqrt(n()),
    n_seeds = n(),
    .groups = "drop"
  ) 

library(ggplot2)

theme_clean <- theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.4),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold")
  )

g3 <- ggplot(
  df_summary_station,
  aes(
    x = mean_coverage,
    y = mean_length,
    color = conformal_method,
    shape = conformal_method
  )
) +
  
  # horizontal SE
  geom_errorbarh(
    aes(
      xmin = mean_coverage - 1.96 * se_coverage,
      xmax = mean_coverage + 1.96 * se_coverage
    ),
    height = 6,
    linewidth = 0.5
  ) +
  
  # vertical SE
  geom_errorbar(
    aes(
      ymin = mean_length - 1.96 * se_length,
      ymax = mean_length + 1.96 * se_length
    ),
    width = 0.01,
    linewidth = 0.5
  ) +
  # points
  geom_point(size = 3, stroke = 1.2) + 
  scale_shape_manual(
    values = c(18, 16, 17, 15)
  ) +
  scale_color_manual(values = method_colors) +
  
  # target coverage
  geom_vline(
    xintercept = 1 - alpha,
    linetype = "dashed",
    linewidth = 0.7,
    color = "grey"
  ) +
  scale_y_continuous(
    breaks = scales::breaks_pretty(n = 3)
  ) +
  labs(
    x = "Coverage",
    y = "Length",
    title = "Station 640"
  ) +
  
  coord_cartesian(xlim = c(0.675, 0.925)) +
  # theme_clean +
  theme_aoas() +
  theme(
    legend.position = "none"
  )
g3

ggsave("Figures/CovLen_Station.pdf", g3, width = 4, height = 2.5)

legend_plot <- ggplot(
  df_summary_station,
  aes(
    x = mean_coverage,
    y = mean_length,
    color = conformal_method,
    shape = conformal_method
  )
) +
  geom_point(size = 3, stroke = 1.2) +
  scale_shape_manual(values = c(18, 16, 17, 15)) +
  scale_color_manual(values = method_colors) +
  theme_aoas() +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.spacing.y = unit(0.05, "cm"),
    legend.key.height = unit(0.8, "cm"),
    legend.text = element_text(size = 14)
  )

legend <- get_legend(legend_plot)

ggsave(
  "Figures/Legend_Methods.pdf",
  legend,
  width = 8,
  height = 1,
)

legend_plot <- ggplot(
  df_summary_station,
  aes(
    x = mean_coverage,
    y = mean_length,
    color = conformal_method,
    shape = conformal_method
  )
) +
  geom_point(size = 3, stroke = 1.2) +
  scale_shape_manual(values = c(18, 16, 17, 15)) +
  scale_color_manual(values = method_colors) +
  theme_aoas() +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    legend.spacing.y = unit(0.05, "cm"),
    legend.key.height = unit(0.8, "cm"),
    legend.text = element_text(size = 14)
  )

legend <- get_legend(legend_plot)

ggsave(
  "Figures/Legend_right_Methods.pdf",
  legend,
  width = 2,
  height = 4,
)

# =====================================================
# STATION 406 NODE
# =====================================================
df_station <- df_time %>% filter(node == 406)
station_data_path <- "Data/dataset_address_main.csv"
dataset_station <- read.csv(station_data_path) %>% 
  mutate(date = as.Date(Date)) %>% preparation_data_regions_stations(group_col=Station.ID) %>%
  filter(Station.ID == 406) %>%
  filter(date %in% test_windows$Date)

df_station <- df_station %>%
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
# PREDICTION INTERVALS
# =====================================================
df_interval <- df_station %>%
  filter(
    projection_method %in% c("refined OLS", "Nested"),
    !is.na(lower),
    !is.na(upper)
  ) %>%
  
  mutate(
    conformal_method = recode(
      conformal_method,
      "CP_MNR_naive" = "CP",
      "CP_MNR" = "CP-MNR",
      "CP_MNR_Nested_star" = "CP-MNR-Nested*",
      "Adaptive_CP_MNR" = "Adaptive CP-MNR"
    )
  ) %>%
  filter(seed==7) %>% # Just consider the intervals for the first seed
  # Aggregate by date and method (across all nodes)
  group_by(conformal_method, date) %>%
  summarise(
    lower_mean = mean(lower),
    upper_mean = mean(upper),
    center     = mean((lower + upper) / 2),
    length_mean = mean(upper - lower),
    .groups = "drop"
  ) 
df_interval <- df_interval %>%
  left_join(dataset_station, by = join_by(date)) %>%
  filter(!is.na(Consumed_kWh)) %>%
  mutate(covered = (Consumed_kWh >= lower_mean & Consumed_kWh <= upper_mean),
  )

y_limits <- range(
  c(df_interval$lower_mean,
    df_interval$upper_mean,
    df_interval$Consumed_kWh),
  na.rm = TRUE
)
library(ggplot2)

theme_clean <- theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.4),
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )

plot_method <- function(method_name) {
  
  df_m <- df_interval %>%
    filter(conformal_method == method_name)
  
  method_col <- method_colors[method_name]
  
  ggplot(df_m, aes(x = date)) +
    
    # main interval
    geom_ribbon(
      aes(ymin = lower_mean, ymax = upper_mean),
      fill = method_col,
      alpha = 0.30
    ) +
    # visible bounds
    geom_line(
      aes(y = lower_mean),
      color = method_col,
      linewidth = 0.6
    ) +
    geom_line(
      aes(y = upper_mean),
      color = method_col,
      linewidth = 0.6
    ) +
    # Observed value
    geom_line(
      aes(y = Consumed_kWh),
      color = "black",
      linewidth = 0.5,
      alpha = 0.4
    ) +
    geom_point(
      aes(y = Consumed_kWh),
      color = "black",
      size = 0.5,
      alpha = 0.5
    ) +
    coord_cartesian(ylim = y_limits) +
    labs(
      title = method_name,
      x = "Date",
      y = "Consumed kWh"
    ) +
    theme_aoas()
}
p1 <- plot_method("CP")
p2 <- plot_method("CP-MNR")
p3 <- plot_method("CP-MNR-Nested*")
p4 <- plot_method("Adaptive CP-MNR")

ggsave("Figures/interval_Station_CP.pdf", p1, width = 4, height = 3)
ggsave("Figures/interval_Station_CP_MNR.pdf", p2, width = 4, height = 3)
ggsave("Figures/interval_Station_CP_MNR_Nested.pdf", p3, width = 4, height = 3)
ggsave("Figures/interval_Station_Adaptive_CP_MNR.pdf", p4, width = 4, height = 3)

p1
p2
p3
p4

# =====================================================
# TIME AVERAGE
# =====================================================

df_summary_station <- df_station %>%
  filter(
    projection_method %in% c("refined OLS", "Nested"),
    !is.na(covered),
    !is.na(length)
  ) %>%
  
  mutate(
    conformal_method = recode(
      conformal_method,
      "CP_MNR_naive" = "CP",
      "CP_MNR" = "CP-MNR",
      "CP_MNR_Nested_star" = "CP-MNR-Nested*",
      "Adaptive_CP_MNR" = "Adaptive CP-MNR"
    )
  ) %>% mutate(
    conformal_method = factor(
      conformal_method,
      levels = c(
        "CP",
        "CP-MNR",
        "Adaptive CP-MNR",
        "CP-MNR-Nested*"
      )
    )
  ) %>% group_by(conformal_method, seed) %>%
  summarise(
    coverage = mean(covered),
    avg_length = mean(length),
    .groups = "drop"
  ) %>%
  group_by(conformal_method) %>%
  summarise(
    mean_coverage = mean(coverage),
    se_coverage   = sd(coverage) / sqrt(n()),
    mean_length   = mean(avg_length),
    se_length     = sd(avg_length) / sqrt(n()),
    n_seeds = n(),
    .groups = "drop"
  ) 

library(ggplot2)

theme_clean <- theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.4),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold")
  )

g3 <- ggplot(
  df_summary_station,
  aes(
    x = mean_coverage,
    y = mean_length,
    color = conformal_method,
    shape = conformal_method
  )
) +
  
  # horizontal SE
  geom_errorbarh(
    aes(
      xmin = mean_coverage - 1.96 * se_coverage,
      xmax = mean_coverage + 1.96 * se_coverage
    ),
    height = 6,
    linewidth = 0.5
  ) +
  
  # vertical SE
  geom_errorbar(
    aes(
      ymin = mean_length - 1.96 * se_length,
      ymax = mean_length + 1.96 * se_length
    ),
    width = 0.01,
    linewidth = 0.5
  ) +
  # points
  geom_point(size = 3, stroke = 1.2) + 
  scale_shape_manual(
    values = c(18, 16, 17, 15)
  ) +
  scale_color_manual(values = method_colors) +
  
  # target coverage
  geom_vline(
    xintercept = 1 - alpha,
    linetype = "dashed",
    linewidth = 0.7,
    color = "grey"
  ) +
  scale_y_continuous(
    breaks = scales::breaks_pretty(n = 3)
  ) +
  labs(
    x = "Coverage",
    y = "Length",
    title = "Station 406"
  ) +
  
  coord_cartesian(xlim = c(0.675, 0.925)) +
  # theme_clean +
  theme_aoas() +
  theme(
    legend.position = "none"
  )
g3

ggsave("Figures/CovLen_Station.pdf", g3, width = 4, height = 2.5)

legend_plot <- ggplot(
  df_summary_station,
  aes(
    x = mean_coverage,
    y = mean_length,
    color = conformal_method,
    shape = conformal_method
  )
) +
  geom_point(size = 3, stroke = 1.2) +
  scale_shape_manual(values = c(18, 16, 17, 15)) +
  scale_color_manual(values = method_colors) +
  theme_aoas() +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.spacing.y = unit(0.05, "cm"),
    legend.key.height = unit(0.8, "cm"),
    legend.text = element_text(size = 14)
  )

legend <- get_legend(legend_plot)

ggsave(
  "Figures/Legend_Methods.pdf",
  legend,
  width = 8,
  height = 1,
)

legend_plot <- ggplot(
  df_summary_station,
  aes(
    x = mean_coverage,
    y = mean_length,
    color = conformal_method,
    shape = conformal_method
  )
) +
  geom_point(size = 3, stroke = 1.2) +
  scale_shape_manual(values = c(18, 16, 17, 15)) +
  scale_color_manual(values = method_colors) +
  theme_aoas() +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    legend.spacing.y = unit(0.05, "cm"),
    legend.key.height = unit(0.8, "cm"),
    legend.text = element_text(size = 14)
  )

legend <- get_legend(legend_plot)

ggsave(
  "Figures/Legend_right_Methods.pdf",
  legend,
  width = 2,
  height = 4,
)

# =====================================================
# STATION 386 NODE
# =====================================================
df_station <- df_time %>% filter(node == 386)
station_data_path <- "Data/dataset_address_main.csv"
dataset_station <- read.csv(station_data_path) %>% 
  mutate(date = as.Date(Date)) %>% preparation_data_regions_stations(group_col=Station.ID) %>%
  filter(Station.ID == 386) %>%
  filter(date %in% test_windows$Date)

df_station <- df_station %>%
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
# PREDICTION INTERVALS
# =====================================================
df_interval <- df_station %>%
  filter(
    projection_method %in% c("refined OLS", "Nested"),
    !is.na(lower),
    !is.na(upper)
  ) %>%
  
  mutate(
    conformal_method = recode(
      conformal_method,
      "CP_MNR_naive" = "CP",
      "CP_MNR" = "CP-MNR",
      "CP_MNR_Nested_star" = "CP-MNR-Nested*",
      "Adaptive_CP_MNR" = "Adaptive CP-MNR"
    )
  ) %>%
  filter(seed == 7) %>% # To only give one example
  # Aggregate by date and method (across all nodes)
  group_by(conformal_method, date) %>%
  summarise(
    lower_mean = mean(lower),
    upper_mean = mean(upper),
    center     = mean((lower + upper) / 2),
    length_mean = mean(upper - lower),
    .groups = "drop"
  ) 
df_interval <- df_interval %>%
  left_join(dataset_station, by = join_by(date)) %>%
  filter(!is.na(Consumed_kWh)) %>%
  mutate(covered = (Consumed_kWh >= lower_mean & Consumed_kWh <= upper_mean),
  )

y_limits <- range(
  c(df_interval$lower_mean,
    df_interval$upper_mean,
    df_interval$Consumed_kWh),
  na.rm = TRUE
)
library(ggplot2)

df_observed <- df_interval %>%
  filter(conformal_method == "CP") %>%
  dplyr::select(date, Consumed_kWh) %>%
  distinct()

plot_method <- function(method_name) {
  test_dates <- sort(unique(test_windows$Date))
  display_start <- as.Date("2024-05-27")
  
  df_m <- df_interval %>%
    filter(conformal_method == method_name)
  
  method_col <- method_colors[method_name]
  
  ggplot(df_m, aes(x = date)) +
    scale_x_date(
      limits = range(test_dates)
    ) +
    # main interval
    geom_ribbon(
      aes(ymin = lower_mean, ymax = upper_mean),
      fill = method_col,
      alpha = 0.30
    ) +
    # visible bounds
    geom_line(
      aes(y = lower_mean),
      color = method_col,
      linewidth = 0.6
    ) +
    geom_line(
      aes(y = upper_mean),
      color = method_col,
      linewidth = 0.6
    ) +
    # Observed value
    geom_line(
      data = df_observed,
      aes(y = Consumed_kWh),
      color = "black",
      linewidth = 0.5,
      alpha = 0.4
    ) +
    geom_point(
      data = df_observed,
      aes(y = Consumed_kWh),
      color = "black",
      size = 0.5,
      alpha = 0.5
    ) +
    coord_cartesian(ylim = y_limits) +
    labs(
      title = method_name,
      x = "Date",
      y = "Consumed kWh"
    ) +
    theme_aoas()
}
p1 <- plot_method("CP")
p2 <- plot_method("CP-MNR")
p3 <- plot_method("CP-MNR-Nested*")
p4 <- plot_method("Adaptive CP-MNR")

ggsave("Figures/interval_Station_bis_CP.pdf", p1, width = 4, height = 3)
ggsave("Figures/interval_Station_bis_CP_MNR.pdf", p2, width = 4, height = 3)
ggsave("Figures/interval_Station_bis_CP_MNR_Nested.pdf", p3, width = 4, height = 3)
ggsave("Figures/interval_Station_bis_Adaptive_CP_MNR.pdf", p4, width = 4, height = 3)

p1
p2
p3
p4

# =====================================================
# TIME AVERAGE
# =====================================================

df_summary_station <- df_station %>%
  filter(
    projection_method %in% c("refined OLS", "Nested"),
    !is.na(covered),
    !is.na(length)
  ) %>%
  
  mutate(
    conformal_method = recode(
      conformal_method,
      "CP_MNR_naive" = "CP",
      "CP_MNR" = "CP-MNR",
      "CP_MNR_Nested_star" = "CP-MNR-Nested*",
      "Adaptive_CP_MNR" = "Adaptive CP-MNR"
    )
  ) %>% mutate(
    conformal_method = factor(
      conformal_method,
      levels = c(
        "CP",
        "CP-MNR",
        "Adaptive CP-MNR",
        "CP-MNR-Nested*"
      )
    )
  ) %>% group_by(conformal_method, seed) %>%
  summarise(
    coverage = mean(covered),
    avg_length = mean(length),
    .groups = "drop"
  ) %>%
  group_by(conformal_method) %>%
  summarise(
    mean_coverage = mean(coverage),
    se_coverage   = sd(coverage) / sqrt(n()),
    mean_length   = mean(avg_length),
    se_length     = sd(avg_length) / sqrt(n()),
    n_seeds = n(),
    .groups = "drop"
  ) 

library(ggplot2)

theme_clean <- theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.4),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold")
  )

g3 <- ggplot(
  df_summary_station,
  aes(
    x = mean_coverage,
    y = mean_length,
    color = conformal_method,
    shape = conformal_method
  )
) +
  
  # horizontal SE
  geom_errorbarh(
    aes(
      xmin = mean_coverage - 1.96 * se_coverage,
      xmax = mean_coverage + 1.96 * se_coverage
    ),
    height = 8,
    linewidth = 0.5
  ) +
  
  # vertical SE
  geom_errorbar(
    aes(
      ymin = mean_length - 1.96 * se_length,
      ymax = mean_length + 1.96 * se_length
    ),
    width = 0.02,
    linewidth = 0.5
  ) +
  # points
  geom_point(size = 3, stroke = 1.2) + 
  scale_shape_manual(
    values = c(18, 16, 17, 15)
  ) +
  scale_color_manual(values = method_colors) +
  
  # target coverage
  geom_vline(
    xintercept = 1 - alpha,
    linetype = "dashed",
    linewidth = 0.7,
    color = "grey"
  ) +
  scale_y_continuous(
    breaks = scales::breaks_pretty(n = 3)
  ) +
  labs(
    x = "Coverage",
    y = "Length",
    title = "Station 386"
  ) +
  
  coord_cartesian(xlim = c(0., 0.925)) +
  # theme_clean +
  theme_aoas() +
  theme(
    legend.position = "none"
  )
g3

ggsave("Figures/CovLen_Station_bis.pdf", g3, width = 4, height = 2.5)

legend_plot <- ggplot(
  df_summary_station,
  aes(
    x = mean_coverage,
    y = mean_length,
    color = conformal_method,
    shape = conformal_method
  )
) +
  geom_point(size = 3, stroke = 1.2) +
  scale_shape_manual(values = c(18, 16, 17, 15)) +
  scale_color_manual(values = method_colors) +
  theme_clean +
  theme(
    legend.position = "bottom",
    legend.spacing.y = unit(0.05, "cm"),
    legend.key.height = unit(0.8, "cm"),
    legend.text = element_text(size = 14)
  )

legend <- get_legend(legend_plot)

ggsave(
  "Figures/Legend_Methods.pdf",
  legend,
  width = 8,
  height = 1,
)