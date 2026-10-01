#!/usr/bin/env Rscript

source('Model_training/R/pretraitements.R')
source('Model_training/R/modelisation.R')

library(argparser)

# inputs
p <- arg_parser("Regional EV models script")

p <- add_argument(p, "--seed",      help="Seed for generating windows", default=40, type="integer")
p <- add_argument(p, "--input",     help="Input dataset path",           default="Data/dataset_region_main.csv")
p <- add_argument(p, "--output",    help="Output file path prefix", default="results_new_period/results_regions")
p <- add_argument(p, "--parallel",  help="Enable parallel computation",          default=TRUE, type="logical")
p <- add_argument(p, "--type",      help="Model type: global, local, or all", default="global")
p <- add_argument(p, "--transform", help="Apply the log1p transformation or normalization",    default='log')
p <- add_argument(p, "--lags",      help="Include lag variables",        default=TRUE, type="logical")

argv <- parse_args(p)

seed_windows  <- argv$seed
raw_data_path <- argv$input
output_path   <- argv$output
bool_parallel <- argv$parallel
type          <- argv$type
transform     <- tolower(argv$transform)
lags          <- argv$lags

# Windows
mes_fenetres <- generate_rolling_windows(SEED = seed_windows) 

# Dataset
dataset_regions <- read.csv(raw_data_path) %>%
        mutate(Date = as.Date(Date)) %>% preparation_data_regions_stations(group_col=Region)

# Formulas

# locals
gam_formula="Consumed_kWh ~  Weekday_Holiday_BIS + 
                                s(Posan, bs='cc') + 
                                Lag1 + 
                                Lag7 + 
                                Lag_Mean_1_7 + 
                                n_cp_paid_lag1:max_paid_cp_30days + 
                                s(tmpf_max)"

rf_formula="Consumed_kWh ~ Weekday_Holiday_BIS + 
                            Posan + 
                            Lag1 + 
                            Lag7 + 
                            Lag_Mean_1_7 + 
                            max_paid_cp_30days + 
                            n_cp_paid_lag1 +
                            tmpf_max + 
                            T_15"

mes_variables =  c("Weekday_Holiday_BIS", "Posan", "Lag1", "Lag7", "Lag_Mean_1_7",
                    "n_cp_paid_lag1", "max_paid_cp_30days", "tmpf_max", "T_15")
                    
list_local_param = list(gam_formula_local = gam_formula,
                        rf_formula_local = rf_formula,
                        mes_variables_local = mes_variables)

# globals
if(transform == 'log'){
    if(lags){
        gam_formula= "log_Consumed_kWh ~  Weekday_Holiday_BIS + 
                                          s(Posan, bs='cc') + 
                                          te(tmpf_max, Latitude) + 
                                          relh_mean + 
                                          taux_cp_paid_max30 +
                                          s(taux_rapide, k=5) +
                                          log_max_cp_30days + 
                                          log_Lag_Mean_1_7 + 
                                          log_ncp_lag1 + 
                                          log_lag_7"

        rf_formula="log_Consumed_kWh ~ Weekday_Holiday_BIS + 
                                       Posan + 
                                       tmpf_max + 
                                       Latitude + 
                                       relh_mean + 
                                       log_max_cp_30days + 
                                       taux_cp_paid_max30 + 
                                       taux_lent + 
                                       taux_rapide + 
                                       taux_accelere + 
                                       log_Lag_Mean_1_7 + 
                                       log_ncp_lag1 + 
                                       log_lag_7 + 
                                       log_lag_1"

        mes_variables =  c("Posan", "Weekday_Holiday_BIS", "tmpf_max", 
                           "Latitude", "relh_mean", "taux_cp_paid_max30", 
                           "log_max_cp_30days", "taux_lent", "taux_rapide", 
                           "taux_accelere", "log_Lag_Mean_1_7", 
                           "log_ncp_lag1", "log_lag_7", "log_lag_1")
    
    } else{
        gam_formula= "log_Consumed_kWh ~  Weekday_Holiday_BIS + 
                                          s(Posan, bs='cc') + 
                                          te(tmpf_max, Latitude) + 
                                          relh_mean + 
                                          taux_cp_paid_max30 +
                                          s(taux_rapide, k=5) +
                                          log_max_cp_30days"

        rf_formula= "log_Consumed_kWh ~ Weekday_Holiday_BIS + 
                                       Posan + 
                                       tmpf_max + 
                                       Latitude + 
                                       relh_mean + 
                                       log_max_cp_30days + 
                                       taux_cp_paid_max30 + 
                                       taux_lent + 
                                       taux_rapide + 
                                       taux_accelere"

        mes_variables =  c("Posan", "Weekday_Holiday_BIS", "tmpf_max", 
                           "Latitude", "relh_mean", "taux_cp_paid_max30", 
                           "log_max_cp_30days", "taux_lent", "taux_rapide", 
                           "taux_accelere")
    }

    list_global_param <- list(gam_formula_global = gam_formula,
                                rf_formula_global = rf_formula,
                                xgb_variables_global = mes_variables, 
                                mes_variables_global = mes_variables,
                                target = "log_Consumed_kWh")

}

if(transform == 'no'){
    if(lags){
        gam_formula = "Consumed_kWh ~ Weekday_Holiday_BIS + 
                                        s(Posan, bs='cc') + 
                                        te(tmpf_max, Latitude) + 
                                        s(max_cp_30days, max_paid_cp_30days) +
                                        taux_lent + 
                                        taux_rapide + 
                                        n_cp_lag1:Weekday_Holiday_BIS + 
                                        Lag1 + 
                                        Lag7 + 
                                        Lag_Mean_1_7"

        rf_formula = "Consumed_kWh ~ Weekday_Holiday_BIS + 
                                        Posan + 
                                        tmpf_max + 
                                        Latitude + 
                                        T_15 + 
                                        max_cp_30days + 
                                        max_paid_cp_30days + 
                                        taux_lent + 
                                        taux_rapide + 
                                        taux_accelere + 
                                        n_cp_lag1 + 
                                        Lag1 + 
                                        Lag7 + 
                                        Lag_Mean_1_7"

        mes_variables =  c("Posan", "Weekday_Holiday_BIS", "tmpf_max", "Latitude", "T_15",
                            "max_cp_30days", "max_paid_cp_30days", "taux_lent", "taux_rapide",
                            "taux_accelere", "n_cp_lag1", "Lag1", "Lag7", "Lag_Mean_1_7",
                            "relh_mean")
    } else {
        gam_formula = "Consumed_kWh ~ Weekday_Holiday_BIS + 
                                        s(Posan, bs='cc') + 
                                        te(tmpf_max, Latitude) + 
                                        s(max_cp_30days, max_paid_cp_30days) +
                                        taux_lent + 
                                        taux_rapide"

        rf_formula = "Consumed_kWh ~ Weekday_Holiday_BIS + 
                                        Posan + 
                                        tmpf_max + 
                                        Latitude + 
                                        T_15 + 
                                        max_cp_30days + 
                                        max_paid_cp_30days + 
                                        taux_lent + 
                                        taux_rapide + 
                                        taux_accelere"

        mes_variables =  c("Posan", "Weekday_Holiday_BIS", "tmpf_max", "Latitude", "T_15",
                            "max_cp_30days", "max_paid_cp_30days", "taux_lent", "taux_rapide",
                            "taux_accelere", "relh_mean")
    }
    list_global_param <- list(gam_formula_global = gam_formula,
                                rf_formula_global = rf_formula,
                                xgb_variables_global = mes_variables, 
                                mes_variables_global = mes_variables,
                                target = 'Consumed_kWh')
} 

if(transform == 'norm'){
    if(lags){
        gam_formula="normalized_max_cp_consumed_kWh ~  Weekday_Holiday_BIS + 
                                                       s(Posan, bs='cc') +
                                                       te(tmpf_max, Latitude) + 
                                                       s(taux_lent, taux_cp_paid_max30) + 
                                                       s(taux_rapide) + 
                                                       s(taux_accelere) + 
                                                       n_cp_paid_lag1 + 
                                                       s(norm_lag_1) + 
                                                       norm_lag_7 + 
                                                       norm_Lag_Mean_1_7"

        rf_formula="normalized_max_cp_consumed_kWh ~ Weekday_Holiday_BIS +
                                                     Posan + 
                                                     tmpf_max + 
                                                     Latitude + 
                                                     taux_cp_paid_max30 + 
                                                     taux_lent + 
                                                     taux_rapide + 
                                                     taux_accelere + 
                                                     n_cp_paid_lag1 + 
                                                     norm_lag_1 + 
                                                     norm_lag_7 + 
                                                     norm_Lag_Mean_1_7"

        mes_variables =  c("Posan", "Weekday_Holiday_BIS", "tmpf_max", "Latitude", "taux_cp_paid_max30",
                           "taux_lent", "taux_rapide", "taux_accelere", "n_cp_paid_lag1", "norm_lag_1", "norm_lag_7", 
                           "norm_Lag_Mean_1_7")

    } else {
        gam_formula="normalized_max_cp_consumed_kWh ~  Weekday_Holiday_BIS + 
                                                       s(Posan, bs='cc') +
                                                       te(tmpf_max, Latitude) + 
                                                       s(taux_lent, taux_cp_paid_max30) + 
                                                       s(taux_rapide) + 
                                                       s(taux_accelere)"

        rf_formula="normalized_max_cp_consumed_kWh ~ Weekday_Holiday_BIS +
                                                     Posan + 
                                                     tmpf_max + 
                                                     Latitude + 
                                                     taux_cp_paid_max30 + 
                                                     taux_lent + 
                                                     taux_rapide + 
                                                     taux_accelere"

        mes_variables =  c("Posan", "Weekday_Holiday_BIS", "tmpf_max", "Latitude", "taux_cp_paid_max30",
                           "taux_lent", "taux_rapide", "taux_accelere")

    }
    list_global_param <- list(gam_formula_global = gam_formula,
                                rf_formula_global = rf_formula,
                                xgb_variables_global = mes_variables, 
                                mes_variables_global = mes_variables,
                                target = 'normalized_max_cp_consumed_kWh')
} 


# Parameters selection based on type
switch(type,
  "local"  = { list_global_param <- NULL },
  "global" = { list_local_param <- NULL },
  "all"    = { }, 
    stop("Unknown type: must be local, global, or all")
)

# Run 
results_regions <- compute_regions_forecasts(
                            dataset = dataset_regions,   
                            windows = mes_fenetres,
                            local_param = list_local_param,
                            global_param = list_global_param,
                            parallel_run = bool_parallel)

if((transform == 'log') & !is.null(results_regions$global)){
    cols_to_fix <- c("GLOBAL_GAM", "GLOBAL_RF", "GLOBAL_XGB", "GLOBAL_MIX")
    
    results_regions$global <- results_regions$global %>% 
                                    mutate(across(all_of(cols_to_fix), expm1))
}

if((transform == 'norm') & !is.null(results_regions$global)){
    cols_to_fix <- c("GLOBAL_GAM", "GLOBAL_RF", "GLOBAL_XGB", "GLOBAL_MIX")
    
    results_regions$global <- results_regions$global %>% 
                                    mutate(across(all_of(cols_to_fix), function(x) return(x*max_cp_30days)))
}

# Export 
file_name <- ifelse(type == 'local', paste0(output_path, '_local.RDS'),
               paste0(output_path, '_', type, '_transform_', transform, "_lags_", lags ,'.RDS'))
saveRDS(results_regions, file_name)