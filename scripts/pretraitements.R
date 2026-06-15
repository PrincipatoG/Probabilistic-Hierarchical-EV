library(dtplyr)
library(dplyr) 
library(tidyr) 
library(readr) 
library(stringr)
library(lubridate) 
library(scales)
library(purrr) 
library(magrittr)
library(broom) 
library(stats)
library(data.table)
library(lars)
library(slider)

#################################################################################################
################################### DECOUPAGE TRAIN / TEST / CALIBRATION ########################
#################################################################################################

generate_rolling_windows_old <- function(SEED = 40) {
  
  set.seed(SEED)
  
  first_train_start <- as.Date("2022-10-01") 
  first_test_start  <- as.Date("2024-09-30") # Lundi
  end_all           <- as.Date("2025-10-01")
  
  train_duration_days <- as.numeric(first_test_start - first_train_start)
  
  # Pool de LUNDIS pour la calibration historique
  pool_dates <- seq(first_train_start, first_test_start - 21, by = "day")
  mondays    <- pool_dates[wday(pool_dates) == 2]
  
  selected_monday_offsets <- as.numeric(sample(mondays, 20) - first_train_start)

  # Fenêtres de test (tous les 14 jours)
  test_starts <- seq(first_test_start, end_all - 13, by = "14 days")
  
  # Génération glissante
  list_of_dfs <- map(test_starts, function(t_start) {
    t_end <- t_start + 13 # Fin du test (Dimanche en semaine 2)
    
    current_train_start <- t_start - train_duration_days
    
    # --- CALIBRATION ---
    # 20 semaines aléatoires (historiques)
    calib_hist_starts <- current_train_start + selected_monday_offsets
    
    # La dernière semaine juste avant le test
    # (du lundi t-7 au dimanche t-1)
    calib_extra_start <- t_start - 7
    
    # Fusion des points de départ de calibration
    all_calib_starts <- c(calib_hist_starts, calib_extra_start)
    
    # Génération de toutes les dates de calibration (7 jours par bloc)
    all_calib_dates <- map(all_calib_starts, ~ seq(.x, .x + 6, by = "day")) %>% 
      reduce(c)
    
    df <- tibble(Date = seq(current_train_start, t_end, by = "day")) %>%
      mutate(type = case_when(
        # Test : les 14 jours cibles
        Date >= t_start & Date <= t_end ~ "test",
        # Calibration : les blocs de 7 jours identifiés
        Date %in% all_calib_dates ~ "calibration",
        # Le reste est du train
        TRUE ~ "train"
      ))
    
    return(df)
  })
  
  return(list_of_dfs)
}

generate_rolling_windows <- function(
    SEED = 40,
    test_start = as.Date("2024-02-05"), # Monday
    end_all = as.Date("2024-08-28"), # Sunday
    train_start = as.Date("2022-10-08"),
    calibration_fraction = 0.2
) {
  
  set.seed(SEED)
  
  # Length of the training period before the first test window
  train_duration_days <- as.numeric(test_start - train_start)
  
  # Candidate Mondays for historical calibration blocks
  pool_dates <- seq(train_start, test_start - 21, by = "day")
  mondays <- pool_dates[wday(pool_dates) == 2]
  
  # Number of available weekly calibration blocks
  n_available_blocks <- length(mondays)
  
  # Number of calibration blocks selected as a fraction
  n_calibration_blocks <- max(
    1,
    round(calibration_fraction * n_available_blocks)
  )
  
  selected_monday_offsets <- as.numeric(
    sample(
      mondays,
      size = n_calibration_blocks,
      replace = FALSE
    ) - train_start
  )
  
  # Test windows every 14 days
  test_starts <- seq(test_start, end_all - 13, by = "14 days")
  
  list_of_dfs <- map(test_starts, function(t_start) {
    
    t_end <- t_start + 13
    
    current_train_start <- t_start - train_duration_days
    
    # Historical calibration blocks
    calib_hist_starts <- current_train_start + selected_monday_offsets
    
    # Most recent week before the test window
    calib_extra_start <- t_start - 7
    
    all_calib_starts <- c(
      calib_hist_starts,
      calib_extra_start
    )
    
    # Expand each calibration block to 7 consecutive days
    all_calib_dates <- map(
      all_calib_starts,
      ~ seq(.x, .x + 6, by = "day")
    ) %>%
      reduce(c)
    
    tibble(
      Date = seq(current_train_start, t_end, by = "day")
    ) %>%
      mutate(
        type = case_when(
          Date >= t_start & Date <= t_end ~ "test",
          Date %in% all_calib_dates ~ "calibration",
          TRUE ~ "train"
        )
      )
  })
  
  return(list_of_dfs)
}

##################################################################################################
#################################### PREPARATION COVARIABLES #####################################
##################################################################################################
preparation_data_scotland <- function(dataset){
    df_full <- dataset %>%
        mutate(
            Weekday_Holiday_BIS = ifelse(Weekday_Holiday %in% c('Tuesday', 'Wednesday',
                'Thursday', 'Friday'), "Weekday", Weekday_Holiday),
            T_15 = (tmpf_max >= 60), 
            max_paid_cp_30days = if_else(is.na(max_paid_cp_30days), n_cp_paid_lag1, max_paid_cp_30days))

    # ON REMPLIT TOUJOURS LES NA 
    df_full <- df_full %>% 
        arrange(Date) %>% 
        mutate(
            across(starts_with("Consumed_kWh_lag_"), 
                   .fns = ~ zoo::na.locf(.x, na.rm = FALSE),
                   .names = "{gsub('Consumed_kWh_lag_', 'Lag', .col)}")) %>%
        mutate(
            Lag_Mean_1_7 = rowMeans(across(num_range("Lag", 1:7)), na.rm = TRUE),
            n_cp_paid_lag1 = zoo::na.locf(n_cp_paid_lag1, na.rm = FALSE), 
            Consumed_kWh_SARIMA = if_else(is.na(Consumed_kWh), Lag7, Consumed_kWh),
            Consumed_kWh_SARIMA = zoo::na.locf(Consumed_kWh_SARIMA, na.rm = FALSE)) 
    return(df_full)
}


preparation_data_regions_stations <- function(dataset, group_col){

    df_full <- dataset %>%
        mutate(
            # VARIABLES DE BASE
            Weekday_Holiday_BIS = ifelse(Weekday_Holiday %in% c('Tuesday', 'Wednesday',
                    'Thursday', 'Friday'), "Weekday", Weekday_Holiday),
            T_15 = (tmpf_max >= 60),
            taux_lent = Lent / nbr_cp_infrastructure,
            taux_accelere = Accelere / nbr_cp_infrastructure,
            taux_rapide = (Rapide + UltraRapide) / nbr_cp_infrastructure,
            
            # TRANSFORMATION LOG 
            log_Consumed_kWh = log1p(Consumed_kWh)
        
        ) %>%
  
        # REMPLISSAGE DES NA
        arrange({{group_col}}, Date) %>% 
        group_by({{group_col}}) %>%
        mutate(
            # CONSOMMATIONS DE BASE
            across(starts_with("Consumed_kWh_lag_"), 
                   .fns = ~ zoo::na.locf(.x, na.rm = FALSE),
                   .names = "{gsub('Consumed_kWh_lag_', 'Lag', .col)}"),
            n_cp_lag1 = zoo::na.locf(n_cp_lag1, na.rm = FALSE),
            n_cp_paid_lag1 = zoo::na.locf(n_cp_paid_lag1, na.rm = FALSE),
            Consumed_kWh_SARIMA = if_else(is.na(Consumed_kWh), Lag7, Consumed_kWh),
            Consumed_kWh_SARIMA = zoo::na.locf(Consumed_kWh_SARIMA, na.rm = FALSE),
            
            # LOG CONSOMMATIONS
            across(log_Consumed_kWh, 
                   .fns = purrr::map(1:7, ~ function(x) lag(x, .x)), 
                   .names = "log_lag_{1:7}")) %>%
            mutate(across(starts_with("log_lag_"), .fns = ~ zoo::na.locf(.x, na.rm = FALSE))) %>% 
        
        ungroup() %>%

        # AJOUT DE COLONNES APRES REMPLISSAGE
        mutate(
            # REMPLISSAGE NA 
            max_paid_cp_30days = if_else(is.na(max_paid_cp_30days), n_cp_paid_lag1, max_paid_cp_30days),
            max_cp_30days = if_else(is.na(max_cp_30days), n_cp_lag1, max_cp_30days),
            # TAUX
            taux_cp_paid = n_cp_paid/n_cp,
            taux_cp_paid_lag1 = n_cp_paid_lag1/n_cp_lag1, 
            taux_cp_paid_max30 = max_paid_cp_30days/max_cp_30days,
            indic_taux_cp_paid_max30 = (taux_cp_paid_max30 != 0),
            taux_utilisation = n_cp_lag1/max_cp_30days,
            # LOGS
            log_max_cp_30days = log1p(max_cp_30days),
            log_ncp_lag1 = log1p(n_cp_lag1),
            # LAGS
            Lag_Mean_1_7 = rowMeans(across(num_range("Lag", 1:7)), na.rm = TRUE),
            log_Lag_Mean_1_7 = rowMeans(across(num_range("log_lag_", 1:7)), na.rm = TRUE),
            diff_log_lag = log_lag_1 - log_lag_7,
            diff_lag = Lag1 - Lag7
        )
    
    col_name <- rlang::as_label(rlang::enquo(group_col)) 
    if(col_name == 'Region'){
        # AJOUT NORMALISATION
         df_full <- df_full %>% 
            mutate(normalized_max_cp_consumed_kWh = ifelse(max_cp_30days > 0, Consumed_kWh / max_cp_30days, Consumed_kWh)) %>% 
            arrange({{group_col}}, Date) %>% 
            group_by({{group_col}}) %>%
            mutate(
                across(normalized_max_cp_consumed_kWh, 
                       .fns = purrr::map(1:7, ~ function(x) lag(x, .x)), 
                       .names = "norm_lag_{1:7}")) %>%
            mutate(across(starts_with("norm_lag_"), .fns = ~ zoo::na.locf(.x, na.rm = FALSE))) %>% 
            ungroup() %>%
            mutate(norm_Lag_Mean_1_7 = rowMeans(pick(starts_with("norm_lag_")), na.rm = TRUE))          
        
        # Ajout latitudes 
         lat_coords <- c( 
            'Shetland Islands' = 60.3,
            'Orkney Islands' = 59.0,
            'Na h-Eileanan an Iar' = 57.8,
            'Moray' = 57.6,
            'Highland' = 57.5,
            'Aberdeenshire' = 57.3,
            'Aberdeen City' = 57.2,
            'Angus' = 56.6,
            'Dundee City' = 56.5,
            'Perth and Kinross' = 56.4,
            'Fife' = 56.2,
            'Stirling' = 56.1,
            'Clackmannanshire' = 56.1,
            'Argyll and Bute' = 56.1,
            'Falkirk' = 56.0,
            'East Lothian' = 56.0,
            'City of Edinburgh' = 55.9,
            'West Dunbartonshire' = 55.9,
            'East Dunbartonshire' = 55.9,
            'Inverclyde' = 55.9,
            'West Lothian' = 55.9,
            'Midlothian' = 55.9,
            'Glasgow City' = 55.9,
            'North Lanarkshire' = 55.9,
            'Renfrewshire' = 55.8,
            'East Renfrewshire' = 55.8,
            'South Lanarkshire' = 55.7,
            'North Ayrshire' = 55.7,
            'Scottish Borders' = 55.7,
            'East Ayrshire' = 55.5,
            'South Ayrshire' = 55.4,
            'Dumfries and Galloway' = 55.0)
        
        df_full <- df_full %>% mutate(Latitude = lat_coords[as.character(Region)])
    }

    return(df_full)
}