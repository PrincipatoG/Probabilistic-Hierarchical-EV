library(dtplyr)
library(dplyr) 
library(tidyr) 
library(readr) 
library(stringr)
library(lubridate)
library(ggplot2) 
library(scales)
library(purrr) 
library(magrittr)
library(broom) 
library(stats)
library(forcats)

library(data.table)
library(patchwork)
library(lars)
library(slider)
library(parallel)
library(forecast)
library(mgcv)
library(ranger)
library(xgboost) 

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


# run_all_local_models <- function(dataset, 
#                                  df_window, 
#                                  gam_formula, 
#                                  rf_formula,
#                                  mes_variables) {
#   
#   # 1. Jointure et préparation des variables
#   df_full <- df_window %>% 
#     left_join(dataset, by = 'Date') %>% 
#     arrange(Date) %>% 
#     mutate(
#       Weekday_Holiday_BIS = ifelse(Weekday_Holiday %in% c('Tuesday', 'Wednesday', 'Thursday', 'Friday'), "Weekday", Weekday_Holiday)) %>% 
#     mutate(
#       across(starts_with("Consumed_kWh_lag_"), 
#              .fns = ~ zoo::na.locf(.x, na.rm = FALSE),
#              .names = "{gsub('Consumed_kWh_lag_', 'Lag', .col)}")) %>%
#     mutate(
#       n_cp_paid_lag1 = zoo::na.locf(n_cp_paid_lag1, na.rm = FALSE), # n_cp_paid_lag1 est parfois à NA, donc il est rempli avec la dernière valeur dispo quand c'est le cas
#       Lag_Mean_1_7 = rowMeans(select(., Lag1:Lag7), na.rm = TRUE),
#       T_15 = (tmpf_max >= 60)
#   )
# 
#   # --- PRÉPARATION SARIMA : Remplissage des trous par J-7 ---
#   # On crée une colonne spécifique pour ARIMA pour ne pas polluer les autres modèles
#   df_full <- df_full %>%
#     mutate(
#       Consumed_kWh_SARIMA = if_else(is.na(Consumed_kWh), Lag7, Consumed_kWh),
#       Consumed_kWh_SARIMA = zoo::na.locf(Consumed_kWh_SARIMA, na.rm = FALSE)
#     )
# 
#   # --- METHODE 1 : ARIMA ---
#   res_arima <- tryCatch({
#     # On utilise la colonne nettoyée au J-7 (pour SARIMA on utilie train et calib pour ne pas casser continuité temporelle)
#     train_vec_clean <- df_full$Consumed_kWh_SARIMA[df_full$type == 'train' | 
#                                                    df_full$type == 'calibration']
#     
#     # Fréquence 7 jours
#     fit_train <- auto.arima(ts(train_vec_clean, frequency = 7))
#     
#     # Application du modèle sur toute la fenêtre (avec les données nettoyées)
#     series_full_clean <- df_full$Consumed_kWh_SARIMA
#     series_full_clean[is.na(series_full_clean)] <- 0 
#     
#     fit_applied <- Arima(series_full_clean, model = fit_train)
#     list(model = fit_train, pred = as.numeric(fitted(fit_applied)), status = "success")
#   }, error = function(e){
#     message("ERREUR ARIMA détectée : ", e$message)
#     list(model=NULL, pred = rep(NA, nrow(df_full)), status = "error")
#   }) 
# 
#   # --- METHODE 2 : GAM ---
#   m_gam <- NULL
#   res_gam <- tryCatch({
#     df_train_gam <- df_full %>% filter(type == "train")
#     m_gam <- gam(as.formula(gam_formula), data = df_train_gam %>% select(Consumed_kWh, mes_variables) %>% drop_na())
#     list(model=m_gam, pred = as.numeric(predict(m_gam, newdata = df_full)), status = "success")
#   }, error = function(e){
#     message("ERREUR GAM détectée : ", e$message)
#     list(model=NULL, pred = rep(NA, nrow(df_full)), status = "error")}) 
# 
#   # --- METHODE 3 : RANDOM FOREST (RANGER) ---
#   rf_fit <- NULL
#   rf_formula <- as.formula(rf_formula)
#   res_rf <- tryCatch({
#     df_train_rf <- df_full %>% filter(type == "train") 
#     
#     rf_fit <- ranger(as.formula(rf_formula), 
#                      data = df_train_rf %>% select(Consumed_kWh, mes_variables) %>% drop_na(), 
#                      num.trees = 500, 
#                      importance = 'impurity')
#     
#     list(model=rf_fit, pred = predict(rf_fit, data = df_full)$predictions, status = "success")
#   }, error = function(e) {
#     message("ERREUR RF détectée : ", e$message)
#     list(model=NULL, pred = rep(NA, nrow(df_full)), status = "error")})
# 
#   # --- RETOUR FINAL ---
#   df_full <- df_full %>%
#     mutate(
#       LOCAL_ARIMA = res_arima$pred,
#       LOCAL_GAM = res_gam$pred,
#       LOCAL_RF = res_rf$pred,
#       Status_ARIMA = res_arima$status,
#       Status_GAM = res_gam$status,
#       Status_RF = res_rf$status
#     ) %>%
#     select(Date, type, Consumed_kWh, Lag1, Lag7, LOCAL_ARIMA, LOCAL_GAM, LOCAL_RF, everything())
# 
#   return(list('all_predict' = df_full,
#               'arima_fit' = res_arima$model,
#               'gam_fit' = res_gam$model,
#               'rf_fit' = res_rf$model))
# }

# run_all_global_models <- function(dataset, 
#                                   df_window, 
#                                   gam_formula, 
#                                   rf_formula,
#                                   xgb_var,
#                                   mes_variables,
#                                   group_col,
#                                   target = "Consumed_kWh") {
#   
#   # 1. Jointure et préparation
#   df_full <- df_window %>% 
#   left_join(dataset, by = 'Date') %>% 
#   arrange(group_col, Date) %>% 
#   group_by(across(all_of(group_col))) %>% 
#   mutate(n_cp_lag7 = lag(n_cp, 7)) %>%
#   mutate(
#     across(starts_with("Consumed_kWh_lag_"), 
#            .fns = ~ zoo::na.locf(.x, na.rm = FALSE),
#            .names = "{gsub('Consumed_kWh_lag_', 'Lag', .col)}"),
#     n_cp_paid_lag1 = zoo::na.locf(n_cp_paid_lag1, na.rm = FALSE),
#     n_cp_lag1 = zoo::na.locf(n_cp_lag1, na.rm = FALSE),
#     n_cp_lag7 = zoo::na.locf(n_cp_lag7, na.rm = FALSE)
#   ) %>% 
#   ungroup() %>% 
#   mutate(
#     Weekday_Holiday_BIS = ifelse(Weekday_Holiday %in% c('Tuesday', 'Wednesday', 'Thursday', 'Friday'), "Weekday", Weekday_Holiday),
#     T_15 = (tmpf_max >= 60),
#     Lag_Mean_1_7 = rowMeans(pick(Lag1:Lag7), na.rm = TRUE),
#     taux_cp_paid = n_cp_paid/n_cp,
#     taux_cp_paid_lag1 = n_cp_paid_lag1/n_cp_lag1, 
#     taux_cp_paid_max30 = max_paid_cp_30days/max_cp_30days,
#     taux_lent = Lent / nbr_cp_infrastructure,
#     taux_accelere = Accelere / nbr_cp_infrastructure,
#     taux_rapide = (Rapide + UltraRapide) / nbr_cp_infrastructure 
#   ) %>%
#   mutate(limit_taux_30days = taux_cp_paid_max30 < 0.25)
#   
# 
#   if (target == 'normalized_lag1_cp_consumed_kWh') {
#     df_full <- df_full %>% 
#       mutate(normalized_lag1_cp_consumed_kWh = ifelse(n_cp_lag1 > 0, Consumed_kWh / n_cp_lag1, Consumed_kWh)) %>% 
#       arrange(across(all_of(group_col)), Date) %>% 
#       group_by(all_of(group_col)) %>%
#       mutate(across(normalized_lag1_cp_consumed_kWh, 
#             .fns = map(1:30, ~ function(x) lag(x, .x)), 
#             .names = "norm_lag_{1:30}")) %>%
#       ungroup() %>%
#       mutate(norm_Lag_Mean_1_7 = rowMeans(pick(norm_lag_1:norm_lag_7), na.rm = TRUE))
#   }
# 
#   if (target == 'normalized_lag7_cp_consumed_kWh') {
#     df_full <- df_full %>% 
#       mutate(normalized_lag7_cp_consumed_kWh = ifelse(n_cp_lag7 > 0, Consumed_kWh / n_cp_lag7, Consumed_kWh)) %>% 
#       arrange(across(all_of(group_col)), Date) %>% 
#       group_by(all_of(group_col)) %>%
#       mutate(across(normalized_lag7_cp_consumed_kWh, 
#             .fns = map(1:30, ~ function(x) lag(x, .x)), 
#             .names = "norm_lag_{1:30}")) %>%
#       ungroup() %>%
#       mutate(norm_Lag_Mean_1_7 = rowMeans(pick(norm_lag_1:norm_lag_7), na.rm = TRUE))
#   }
# 
#   if (target == 'normalized_max_cp_consumed_kWh') {
#     df_full <- df_full %>% 
#       mutate(normalized_max_cp_consumed_kWh = ifelse(max_cp_30days > 0, Consumed_kWh / max_cp_30days, Consumed_kWh)) %>% 
#       arrange(across(all_of(group_col)), Date) %>% 
#       group_by(all_of(group_col)) %>%
#       mutate(across(normalized_max_cp_consumed_kWh, 
#             .fns = map(1:30, ~ function(x) lag(x, .x)), 
#             .names = "norm_lag_{1:30}")) %>%
#       # across(starts_with("norm_lag"), # Si on veut eviter les NA dans les prev, faire ça !!!
#       #      .fns = ~ zoo::na.locf(.x, na.rm = FALSE),
#       #      .names = "{gsub('norm_', 'lag', .col)}")
#       ungroup() %>%
#       mutate(norm_Lag_Mean_1_7 = rowMeans(pick(norm_lag_1:norm_lag_7), na.rm = TRUE))
#   }
# 
# if (target == 'normalized_mean_var_consumed_kWh') {
# 
#     df_full <- df_full %>% 
#       group_by(across(all_of(group_col))) %>% 
#       mutate(
#         # On calcule les stats de référence
#         m_train = mean(Consumed_kWh[type == "train"], na.rm = TRUE),
#         s_train = sd(Consumed_kWh[type == "train"], na.rm = TRUE)) %>%
#       mutate(
#         # On utilise DIRECTEMENT ces stats pour créer la cible
#         normalized_mean_var_consumed_kWh = case_when(
#           is.na(m_train) ~ Consumed_kWh, 
#           is.na(s_train) | s_train == 0 ~ Consumed_kWh - m_train, 
#           TRUE ~ (Consumed_kWh - m_train) / s_train  
#         )
#       ) %>% 
#       ungroup() %>%
#       # On continue avec les lags sur la nouvelle colonne
#       arrange(across(all_of(group_col)), Date) %>% 
#       group_by(across(all_of(group_col))) %>%
#       mutate(across(all_of(target), 
#              .fns = map(1:30, ~ function(x) lag(x, .x)), 
#              .names = "norm_lag_{1:30}")) %>%
#       ungroup() %>%
#       mutate(norm_Lag_Mean_1_7 = rowMeans(pick(num_range("norm_lag_", 1:7)), na.rm = TRUE))
#   }
#   
#   # On ajoute aussi les latitudes 
#   if(group_col == "Region"){
#       lat_coords <- c(
#         'Shetland Islands' = 60.3, 'Orkney Islands' = 59.0, 'Na h-Eileanan an Iar' = 58.2,
#         'Highland' = 57.5, 'Moray' = 57.4, 'Aberdeenshire' = 57.2, 'Aberdeen City' = 57.1,
#         'Angus' = 56.7, 'Perth and Kinross' = 56.5, 'Dundee City' = 56.5, 'Argyll and Bute' = 56.2,
#         'Fife' = 56.2, 'Stirling' = 56.2, 'Clackmannanshire' = 56.1, 'Falkirk' = 56.0,
#         'West Lothian' = 55.9, 'City of Edinburgh' = 55.9, 'East Lothian' = 55.9,
#         'West Dunbartonshire' = 55.9, 'East Dunbartonshire' = 55.9, 'North Lanarkshire' = 55.8,
#         'Glasgow City' = 55.8, 'Renfrewshire' = 55.8, 'Inverclyde' = 55.9, 'Midlothian' = 55.8,
#         'North Ayrshire' = 55.7, 'East Renfrewshire' = 55.7, 'South Lanarkshire' = 55.6,
#         'South Ayrshire' = 55.4, 'East Ayrshire' = 55.4, 'Scottish Borders' = 55.5,
#         'Dumfries and Galloway' = 55.0)
# 
#   df_full <- df_full %>% mutate(Latitude = lat_coords[as.character(Region)]) %>% 
#                          mutate(Zone = case_when(
#                                                   Latitude > 57.0  ~ "NORD",
#                                                   Latitude > 55.7  ~ "CENTRE",
#                                                   TRUE             ~ "SUD"
#                                                 ))
#   }
# 
# 
#   # # --- METHODE 1 : GAM ---
#   #   m_gam <- NULL
#   #   res_gam <- tryCatch({
#   #     df_train_gam <- df_full %>% filter(type == "train") %>% select(all_of(target), all_of(mes_variables)) %>% drop_na()
#   #     m_gam <- bam(as.formula(gam_formula), data = df_train_gam )
#   #   list(model = m_gam, pred = as.numeric(predict(m_gam, newdata = df_full)), status = "success")
#   # }, error = function(e){
#   #   message("ERREUR GAM détectée : ", e$message)
#   #   list(model = NULL, pred = rep(NA, nrow(df_full)), status = "error")}) 
# 
#   # --- METHODE 2 : RANDOM FOREST ---
#   rf_fit <- NULL
#   rf_formula <- as.formula(rf_formula)
#   res_rf <- tryCatch({
#     df_train_rf <- df_full %>% filter(type == "train") %>% select(all_of(target), all_of(mes_variables)) %>% drop_na()
#     
#     rf_fit <- ranger(as.formula(rf_formula), 
#                      data = df_train_rf , 
#                      num.trees = 1000, 
#                      importance = 'impurity')
#     
#     list(model = rf_fit, pred = predict(rf_fit, data = df_full)$predictions, status = "success")
#   }, error = function(e) {
#     message("ERREUR RF détectée : ", e$message)
#     list(model = NULL, pred = rep(NA, nrow(df_full)), status = "error")})
#   
#   # --- METHODE 3 : XGBOOST ----
#   # xgb_fit <- NULL
#   # res_xgb <- tryCatch({
#   #   df_train_xgb <- df_full %>% 
#   #     filter(type == "train") %>% 
#   #     select(all_of(target), all_of(xgb_var)) %>% 
#   #     drop_na()
#   #   formula_ohe <- as.formula(paste("~ -1 +", paste(xgb_var, collapse = " + ")))
#   #   train_x <- model.matrix(formula_ohe, data = df_train_xgb)
#   #   train_y <- df_train_xgb[[target]]
#   #   
#   #   # Matrice pour la prédiction finale (tout le dataset)
#   #   rows_keep <- which(complete.cases(df_full[, xgb_var]))
#   #   
#   #   # 2. On crée la matrice uniquement sur ces lignes
#   #   # Comme on a filtré avant, model.matrix ne supprimera rien de plus
#   #   full_x_clean <- model.matrix(formula_ohe, data = df_full[rows_keep, ])
#   #   
#   #   # 3. Entraînement (déjà fait sur train_x)
#   #   params <- list(
#   #     objective = "reg:squarederror",
#   #     eta = 0.05,              # Apprentissage plus lent et précis
#   #     max_depth = 8,           # Arbres plus profonds pour capturer plus d'interactions
#   #     min_child_weight = 2,    # Évite de créer des feuilles sur des cas trop isolés
#   #     subsample = 0.8,         # Utilise 80% des lignes au hasard par arbre
#   #     colsample_bytree = 0.7,  # Utilise 70% des variables au hasard par arbre
#   #     lambda = 1,              # Régularisation L2
#   #     alpha = 0.5              # Régularisation L1
#   #   )
#   # 
#   #   xgb_fit <- xgboost(
#   #     params = params,
#   #     data = train_x, 
#   #     label = train_y, 
#   #     nrounds = 500,           # Augmenté car eta est plus petit
#   #     verbose = 0, 
#   #     nthread = 4,             # Profite du HPC pour paralléliser le calcul
#   #     early_stopping_rounds = 20 # S'arrête si le modèle n'évolue plus (nécessite un set de validation)
#   #   )
#   #   
#   #   # 4. PRÉDICTION AVEC MAPPING
#   #   # On crée un vecteur de NA de la taille PARFAITE (23808)
#   #   preds_finales <- rep(NA, nrow(df_full))
#   #   
#   #   # On prédit pour les lignes valides
#   #   preds_valides <- predict(xgb_fit, full_x_clean)
#   #   
#   #   # On injecte les prédictions aux bons index
#   #   preds_finales[rows_keep] <- preds_valides
#   # 
#   #  list(model = xgb_fit, pred = preds_finales, status = "success")
#   #   
#   # }, error = function(e) {
#   #   message("ERREUR XGBoost détectée : ", e$message)
#   #   list(model=NULL, pred = rep(NA, nrow(df_full)), status = "error")
#   # })
# 
#   # --- RETOUR FINAL ---
#   # On réassigne les prédictions au dataframe
#   df_full <- df_full %>%
#     mutate(
#       # GLOBAL_GAM = res_gam$pred,
#       GLOBAL_RF = res_rf$pred,
#       # GLOBAL_XGB = res_xgb$pred,
#       # Status_GAM = res_gam$status,
#       Status_RF = res_rf$status,
#       # Status_XGB = res_xgb$status
#     ) %>%
#     # Utilisation de all_of(target) ici aussi pour le select final
#     select(Date, type, all_of(target), GLOBAL_RF, everything()) #, GLOBAL_GAM, GLOBAL_RF, GLOBAL_XGB, everything())
# 
#   return(list('all_predict' = df_full,
#               'rf_fit' = res_rf$model ))
#               # 'gam_fit' = res_gam$model,
#               # 'xgb_fit' = res_xgb$model))
# }

# compute_all_forecasts <- function(dataset, 
#                                   windows, 
#                                   local_param = NULL, 
#                                   global_param = NULL,
#                                   target='Consumed_kWh',
#                                   group_col = NULL) {
#   
#   n_cores <- parallel::detectCores() - 1
#   
#   final_local_df <- NULL
#   if(!is.null(local_param)){
#     if (is.null(group_col)) {
#       # --- I. MODELE SCOTLAND ---
#       message(paste("Lancement du modèle national sur", length(windows), "fenetres..."))
#       
#       results <- parallel::mclapply(windows, function(w) {
#         run_all_local_models(dataset, w, local_param$gam_formula_local, local_param$rf_formula_local, local_param$mes_variables_local)$all_predict
#       }, mc.cores = n_cores)
# 
#       return(dplyr::bind_rows(results, .id = "window_id"))
#       
#     } else {
#       # --- II. MODELES LOCAUX ---
#       message(paste("Modeles locaux par", group_col, "..."))
#       split_data <- dataset %>% split(.[[group_col]])
#       
#       message(paste("Calcul lancé pour", length(split_data), "séries", n_cores, "coeurs..."))
#       
#       # On parallélise sur les séries
#       local_results <- parallel::mclapply(names(split_data), function(name) {
#         entity_data <- split_data[[name]]
#         purrr::map_dfr(windows, ~ run_all_local_models(entity_data, .x, local_param$gam_formula_local, local_param$rf_formula_local, local_param$mes_variables_local)$all_predict,
#         .id = "window_id")
#         
#       }, mc.cores = n_cores)
#       
#       # On réassigne les noms des séries avant de fusionner
#       names(local_results) <- names(split_data)
#       
#       # Fusion finale de toutes les stations en un seul dataframe géant
#       final_local_df <- dplyr::bind_rows(local_results, .id = group_col)
#     }}
# 
#   
#   final_global_df <- NULL
#   if(!is.null(global_param)){
#       # --- III. MODELES GLOBAUX --- 
#       message(paste("Modeles globaux par", group_col, "..."))
# 
#       global_results <- parallel::mclapply(windows, function(window) {
#         run_all_global_models(dataset, window, global_param$gam_formula_global, global_param$rf_formula_global, global_param$xgb_variables_global, global_param$mes_variables_global, group_col, target)$all_predict
#       }, mc.cores = n_cores)
#       
#       # On réassigne les noms des stations avant de fusionner
#       final_global_df <- dplyr::bind_rows(global_results, .id = "window_id")
#   }
#   return(list('local' = final_local_df,
#               'global' = final_global_df))
# }


build_vectors <- function(period,
                          res_nat,
                          res_reg,
                          res_sta,
                          model_names) {
  
  # --- Filtrage une seule fois ---
  res_nat <- res_nat[res_nat$Date %in% period, ]
  res_reg <- res_reg[res_reg$Date %in% period, ]
  res_sta <- res_sta[res_sta$Date %in% period, ]
  
  # Ordonner (important pour alignement temporel)
  res_nat <- res_nat[order(res_nat$Date), ]
  res_reg <- res_reg[order(res_reg$Date), ]
  res_sta <- res_sta[order(res_sta$Date), ]
  
  # --- NATIONAL (vecteurs simples) ---
  Y_nat    <- res_nat$Consumed_kWh
  Yhat_nat <- res_nat[[model_names[[1]]]]
  
  # --- REGIONS (pivot → matrix) ---
  Y_reg <- tidyr::pivot_wider(
    res_reg,
    id_cols = Date,
    names_from = Region,
    values_from = Consumed_kWh
  )
  Yhat_reg <- tidyr::pivot_wider(
    res_reg,
    id_cols = Date,
    names_from = Region,
    values_from = .data[[model_names[[2]]]]
  )
  
  # drop Date + conversion rapide
  Y_reg    <- as.matrix(Y_reg[, -1, drop = FALSE])
  Yhat_reg <- as.matrix(Yhat_reg[, -1, drop = FALSE])
  
  # --- STATIONS ---
  Y_sta <- tidyr::pivot_wider(
    res_sta,
    id_cols = Date,
    names_from = Station.ID,
    values_from = Consumed_kWh
  )
  Yhat_sta <- tidyr::pivot_wider(
    res_sta,
    id_cols = Date,
    names_from = Station.ID,
    values_from = .data[[model_names[[3]]]]
  )
  
  Y_sta    <- as.matrix(Y_sta[, -1, drop = FALSE])
  Yhat_sta <- as.matrix(Yhat_sta[, -1, drop = FALSE])
  
  # --- CONCAT (une seule fois) ---
  Y    <- cbind(Y_nat, Y_reg, Y_sta)
  Yhat <- cbind(Yhat_nat, Yhat_reg, Yhat_sta)
  
  # --- Construction efficace ---
  n <- nrow(Y)
  
  vectors <- vector("list", n)
  for (t in seq_len(n)) {
    vectors[[t]] <- list(
      y = Y[t, ],
      yhat = Yhat[t, ]
    )
  }
  
  names(vectors) <- period
  
  return(vectors)
}
