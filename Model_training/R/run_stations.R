#!/usr/bin/env Rscript

source('Model_training/R/pretraitements.R')
source('Model_training/R/modelisation.R')

library(argparser)

# inputs
p <- arg_parser("Station EV models script")

p <- add_argument(p, "--seed", help="Seed for generating windows", default=40, type="integer")
p <- add_argument(p, "--input", help="Input dataset path", default="Data/dataset_address_main.csv")
p <- add_argument(p, "--output", help="Output file path prefix", default="results_new_period/results_stations")
p <- add_argument(p, "--parallel", help="Enable parallel computation", default=TRUE, type="logical")
p <- add_argument(p, "--type", help="Model type: global, local, or all", default="global")
p <- add_argument(p, "--lags", help="Include lag variables", default=TRUE, type="logical")
p <- add_argument(p, "--max_pct_na", help="Remove series with too many NAs from training", default=1, type="numeric")
p <- add_argument(p, "--nbr_series_train", help="Reduce the number of training series",  default=NULL, type="integer")

argv <- parse_args(p)

seed_windows  <- argv$seed
raw_data_path <- argv$input
output_path   <- argv$output
bool_parallel <- argv$parallel
type          <- argv$type
lags          <- argv$lags
max_pct_na    <- argv$max_pct_na
nbr_series_train <- argv$nbr_series_train

print(nbr_series_train)

# Windows
mes_fenetres <- generate_rolling_windows(SEED = seed_windows) 

# Dataset
dataset_stations <- data.table::fread(raw_data_path) %>% 
        mutate(Date = as.Date(Date)) %>% preparation_data_regions_stations(group_col=Station.ID)

# Formulas

# locals
gam_formula="Consumed_kWh ~  Weekday_Holiday_BIS + 
                             s(Posan, bs='cc') +
                             indic_taux_cp_paid_max30 +
                             s(tmpf_max) + 
                             relh_min + 
                             max_cp_30days + 
                             Lag_Mean_1_7 + 
                             Lag1 + 
                             Lag7"

rf_formula="Consumed_kWh ~ Weekday_Holiday_BIS + 
                           Posan + 
                           indic_taux_cp_paid_max30 + 
                           tmpf_max + 
                           relh_min + 
                           max_cp_30days + 
                           Lag_Mean_1_7 + 
                           Lag1 + 
                           Lag7"

mes_variables =  c("Weekday_Holiday_BIS", "Posan", "indic_taux_cp_paid_max30", 
                   "tmpf_max", "relh_min", "max_cp_30days", "Lag_Mean_1_7", "n_cp_lag1",
                   "Lag1", "Lag7")
                    
list_local_param = list(gam_formula_local = gam_formula,
                        rf_formula_local = rf_formula,
                        mes_variables_local = mes_variables)

# globals
if(lags){
    gam_formula="log_Consumed_kWh ~  Weekday_Holiday_BIS + 
                                     s(Posan, bs='cc') + 
                                     indic_taux_cp_paid_max30 + 
                                     log_max_cp_30days + 
                                     taux_lent + 
                                     taux_rapide + 
                                     te(Latitude, tmpf_max) + 
                                     relh_min + 
                                     log_Lag_Mean_1_7 + 
                                     log_ncp_lag1 + 
                                     log_lag_1 + 
                                     log_lag_7"

    rf_formula="log_Consumed_kWh ~ Weekday_Holiday_BIS + 
                                   Posan + 
                                   indic_taux_cp_paid_max30 + 
                                   log_max_cp_30days + 
                                   taux_lent + 
                                   taux_rapide + 
                                   tmpf_max + 
                                   relh_min + 
                                   Latitude + 
                                   Longitude +
                                   log_Lag_Mean_1_7 + 
                                   log_lag_1 + 
                                   log_lag_7 + 
                                   log_ncp_lag1"

    mes_variables =  c("Weekday_Holiday_BIS", "Posan", "indic_taux_cp_paid_max30", 
                       "log_max_cp_30days", "taux_lent", "taux_rapide", "tmpf_max",
                       "relh_min", "Latitude", "Longitude", "log_Lag_Mean_1_7", 
                       "log_ncp_lag1",  "log_lag_1", "log_lag_7") 


} else{
    gam_formula="log_Consumed_kWh ~  Weekday_Holiday_BIS + 
                                     s(Posan, bs='cc') + 
                                     indic_taux_cp_paid_max30 + 
                                     log_max_cp_30days + 
                                     taux_lent + 
                                     taux_rapide + 
                                     te(Latitude, tmpf_max) + 
                                     relh_min"

    rf_formula="log_Consumed_kWh ~ Weekday_Holiday_BIS + 
                                   Posan + 
                                   indic_taux_cp_paid_max30 + 
                                   log_max_cp_30days + 
                                   taux_lent + 
                                   taux_rapide + 
                                   tmpf_max + 
                                   relh_min + 
                                   Latitude + 
                                   Longitude"

    mes_variables =  c("Weekday_Holiday_BIS", "Posan", "indic_taux_cp_paid_max30", 
                       "log_max_cp_30days", "taux_lent", "taux_rapide", "tmpf_max",
                       "relh_min", "Latitude", "Longitude") 
}

list_global_param <- list(gam_formula_global = gam_formula,
                          rf_formula_global = rf_formula,
                          xgb_variables_global = mes_variables, 
                          mes_variables_global = mes_variables,
                          target = "log_Consumed_kWh",
                          max_pct_na = max_pct_na,
                          nbr_series_train = nbr_series_train) 

# Parameters selection based on type
switch(type,
  "local"  = { list_global_param <- NULL },
  "global" = { list_local_param <- NULL },
  "all"    = { }, 
  stop("Unknown type: must be local, global, or all")
)
# Export 
if(is.na(nbr_series_train)){
  nbr_series_train <- dataset_stations %>% dplyr::select(Station.ID) %>% unique() %>% pull(Station.ID) %>% length()
}
my_folder <- paste0(output_path, "_type_", type, "_lags_", lags , "_", nbr_series_train)

# Run 
results_stations <- compute_stations_forecasts(
                            dataset = dataset_stations,   
                            windows = mes_fenetres,
                            output_path = my_folder,
                            local_param = list_local_param,
                            global_param = list_global_param,
                            parallel_run = bool_parallel)