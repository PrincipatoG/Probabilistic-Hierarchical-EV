#!/usr/bin/env Rscript

source('Model_training/R/pretraitements.R')
source('Model_training/R/modelisation.R')

library(argparser)

# inputs
p <- arg_parser("Scotland EV model script")

p <- add_argument(p, "--seed",      help="Seed for generating windows", default=40, type="integer")
p <- add_argument(p, "--input",     help="Input dataset path",           default="Data/dataset_scotland_main.csv")
p <- add_argument(p, "--output",    help="Output file path prefix", default="results_new_period/results_scotland.RDS")
p <- add_argument(p, "--parallel",  help="Enable parallel computation",          default=TRUE, type="logical")

argv <- parse_args(p)

seed_windows  <- argv$seed
raw_data_path <- argv$input
output_path   <- argv$output
bool_parallel <- argv$parallel

# Windows
mes_fenetres <- generate_rolling_windows(SEED = seed_windows)  

# Dataset
dataset_scotland <- read.csv(raw_data_path) %>% 
        mutate(Date = as.Date(Date)) %>% preparation_data_scotland()

# Formulas 
gam_formula <- "Consumed_kWh ~  Weekday_Holiday_BIS + 
                                s(Posan, bs='cc') + 
                                Lag1 + 
                                Lag7 + 
                                Lag_Mean_1_7 + 
                                n_cp_paid_lag1:max_paid_cp_30days + 
                                s(tmpf_max)"

rf_formula <- "Consumed_kWh ~ Weekday_Holiday_BIS + 
                              Posan + 
                              Lag1 + 
                              Lag7 + 
                              Lag_Mean_1_7 + 
                              max_paid_cp_30days + 
                              n_cp_paid_lag1 +
                              tmpf_max + 
                              T_15"

mes_variables <-  c("Weekday_Holiday_BIS", "Posan", "Lag1", "Lag7", "Lag_Mean_1_7", "n_cp_lag1",
                    "n_cp_paid_lag1", "max_paid_cp_30days", "tmpf_max", "T_15")

# Run
results_scotland <- compute_scotland_forecasts(
    dataset = dataset_scotland, 
    windows = mes_fenetres, 
    param = list(gam_formula = gam_formula, 
                 rf_formula = rf_formula, 
                 mes_variables = mes_variables),
    parallel_run = bool_parallel)

# Export 
saveRDS(results_scotland, output_path)
