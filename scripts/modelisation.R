library(parallel)
library(forecast)
library(mgcv)
library(ranger)
library(xgboost) 
library(Matrix)

##################################################################################################
##################################### FONCTIONS ELEMENTAIRES #####################################
##################################################################################################

run_all_local_models <- function(dataset, 
                                 df_window, 
                                 gam_formula, 
                                 rf_formula,
                                 mes_variables,
                                 return_models = FALSE) {
  
  # 1. Jointure initiale
  df_full <- df_window %>% left_join(dataset, by = 'Date')
  
  # --- MÉTHODE 1 : ARIMA ---
  res_arima <- tryCatch({
    train_vec_clean <- df_full$Consumed_kWh_SARIMA[df_full$type %in% c('train', 'calibration')]
    fit_train <- auto.arima(ts(train_vec_clean, frequency = 7))
    
    # On applique sur toute la série pour la continuité
    fit_applied <- Arima(df_full$Consumed_kWh_SARIMA, model = fit_train)
    full_preds <- as.numeric(fitted(fit_applied))
    
    # Filtrage selon le besoin
    pred_out <- if(!return_models) full_preds[df_full$type != "train"] else full_preds
    
    list(model = fit_train, pred = pred_out, status = "success")
  }, error = function(e){
    list(model=NULL, pred = rep(NA, if(!return_models) sum(df_full$type != "train") else nrow(df_full)), status = "error")
  }) 

  # --- MÉTHODE 2 : GAM ---
  res_gam <- tryCatch({
    m_gam <- gam(as.formula(gam_formula), 
                 data = df_full %>% filter(type == "train") %>% drop_na(Consumed_kWh, all_of(mes_variables)))
    
    # Prédiction
    data_to_pred <- if(!return_models) df_full %>% filter(type != "train") else df_full
    list(model=m_gam, pred = as.numeric(predict(m_gam, newdata = data_to_pred)), status = "success")
  }, error = function(e){
    list(model=NULL, pred = rep(NA, if(!return_models) sum(df_full$type != "train") else nrow(df_full)), status = "error")
  }) 

  # --- MÉTHODE 3 : RANDOM FOREST (RANGER) ---
  res_rf <- tryCatch({
    rf_fit <- ranger(as.formula(rf_formula), 
                     data = df_full %>% filter(type == "train") %>% drop_na(Consumed_kWh, all_of(mes_variables)), 
                     num.trees = 500, 
                     num.threads = 1,
                     importance = 'impurity')
    
    # Prédiction à la volée sur le subset
    data_to_pred <- if(!return_models) df_full %>% filter(type != "train") else df_full
    list(model=rf_fit, pred = predict(rf_fit, data = data_to_pred)$predictions, status = "success")
  }, error = function(e) {
    list(model=NULL, pred = rep(NA, if(!return_models) sum(df_full$type != "train") else nrow(df_full)), status = "error")
  })

  # --- ASSEMBLAGE FINAL ---
  # On construit le dataframe final en filtrant df_full si return_models est FALSE
  df_result <- if(!return_models) df_full %>% filter(type != "train") else df_full

  df_result <- df_result %>%
    mutate(
      LOCAL_ARIMA = res_arima$pred,
      LOCAL_GAM   = res_gam$pred,
      LOCAL_RF    = res_rf$pred,
      Status_ARIMA = res_arima$status,
      Status_GAM   = res_gam$status,
      Status_RF    = res_rf$status,
      LOCAL_MIX    = (LOCAL_ARIMA + LOCAL_GAM + LOCAL_RF) / 3
    ) %>%
    dplyr::select(Date, type, Consumed_kWh, LOCAL_ARIMA, LOCAL_GAM, LOCAL_RF, LOCAL_MIX, everything())

  if(!return_models){
    return(list('all_predict' = df_result))
  } else {
    return(list('all_predict' = df_result,
                'arima_fit'   = res_arima$model,
                'gam_fit'     = res_gam$model,
                'rf_fit'      = res_rf$model))
  }
}

run_all_global_models <- function(dataset, 
                                  df_window, 
                                  gam_formula, 
                                  rf_formula,
                                  xgb_var,
                                  mes_variables,
                                  target = "Consumed_kWh",
                                  return_models = FALSE,
                                  num_threads = 1,
                                  max_pct_na = 1,
                                  nbr_series_train = NULL,
                                  SEED = 40) {
  
  message(paste0('window min date test: ', 
          min(df_window %>% filter(type == 'test') %>% pull(Date))))                                  
  df_full <- df_window %>% left_join(dataset, by = 'Date')
  df_full$type_use <- df_full$type

  # --- SAMPLING AVANT LANCEMENT POUR LES STATIONS : OPTIONNEL --- #
  if (!is.null(nbr_series_train) && !is.na(nbr_series_train)){
      set.seed(SEED)
      sample_stations <- dataset %>% 
          group_by(Station.ID) %>% 
          summarize(taux_na = sum(is.na(!!sym(target))) / n()) %>% 
          filter(taux_na <= max_pct_na) %>% 
          pull(Station.ID)
      
      if(length(sample_stations) > nbr_series_train) {
          sample_stations <- sample(sample_stations, nbr_series_train)
      }

      idx_to_change <- which(df_full$type == 'train' & !(df_full$Station.ID %in% sample_stations))
      df_full$type_use[idx_to_change] <- 'test'
  }

  # --- Définition du set de prédiction ---
  df_pred <- if(!return_models) df_full[df_full$type != "train", ] else df_full

  # --- METHODE 1 : GAM ---
  m_gam <- NULL
  res_gam <- tryCatch({
    m_gam <- bam(as.formula(gam_formula), 
                   data = df_full[df_full$type_use == "train", c(target, mes_variables)] %>% drop_na(),
                   discrete = TRUE,
                   nthreads = num_threads)
    list(model = m_gam, pred = as.numeric(predict(m_gam, newdata = df_pred)), status = "success") 
  }, error = function(e){
    message("ERREUR GAM détectée : ", e$message)
    list(model = NULL, pred = rep(NA, nrow(df_pred)), status = "error")}) 
  gc()

  # --- METHODE 2 : RANDOM FOREST ---
  rf_fit <- NULL
  rf_formula <- as.formula(rf_formula)
  res_rf <- tryCatch({  
    rf_fit <- ranger(as.formula(rf_formula), 
                       data = df_full[df_full$type_use == "train", c(target, mes_variables)] %>% drop_na(), 
                       num.trees = 1000,
                       mtry = 5, 
                       num.threads = num_threads,
                       importance = 'impurity')
    list(model = rf_fit, pred = predict(rf_fit, data = df_pred)$predictions, status = "success")
  }, error = function(e) {
    message("ERREUR RF détectée : ", e$message)
    list(model = NULL, pred = rep(NA, nrow(df_pred)), status = "error")})
  gc()

  # --- METHODE 3 : XGBOOST ----
  xgb_fit <- NULL
  res_xgb <- tryCatch({
    df_train_xgb <- df_full[df_full$type_use == "train", c(target, xgb_var)] %>% drop_na()
    formula_ohe <- as.formula(paste("~ -1 +", paste(xgb_var, collapse = " + ")))
    
    train_x <- sparse.model.matrix(formula_ohe, data = df_train_xgb)
    dtrain <- xgb.DMatrix(data = train_x, label = df_train_xgb[[target]])
    
    rm(df_train_xgb, train_x); gc() 

    xgb_fit <- xgb.train(params = list(objective = "reg:squarederror", 
                                       eta = 0.05, 
                                       max_depth = 8, 
                                       min_child_weight = 20,
                                       subsample = 0.8, 
                                       colsample_bytree = 0.6,
                                       lambda = 5,
                                       alpha = 1),
                         data = dtrain,
                         nrounds = 1000)
    rm(dtrain); gc()

    # Prediction
    rows_keep <- which(complete.cases(df_pred[, xgb_var]))
    full_x_sparse <- sparse.model.matrix(formula_ohe, data = df_pred[rows_keep, ])
    preds_valides <- predict(xgb_fit, xgb.DMatrix(full_x_sparse))
    
    preds_finales <- rep(NA, nrow(df_pred))
    preds_finales[rows_keep] <- preds_valides
    
    rm(full_x_sparse); gc()
    list(model = xgb_fit, pred = preds_finales, status = "success")
  }, error = function(e) {
    message("ERREUR XGBoost détectée : ", e$message); list(model=NULL, pred = rep(NA, nrow(df_pred)), status = "error")
  })

  # --- RETOUR FINAL ---
  df_result <- df_pred %>%
    mutate(
      GLOBAL_GAM = res_gam$pred,
      GLOBAL_RF = res_rf$pred,
      GLOBAL_XGB = res_xgb$pred,
      Status_GAM = res_gam$status,
      Status_RF = res_rf$status,
      Status_XGB = res_xgb$status
    ) %>%
    dplyr::select(Date, type, all_of(target), GLOBAL_GAM, GLOBAL_RF, GLOBAL_XGB, everything()) %>% 
    mutate(GLOBAL_MIX = (GLOBAL_GAM + GLOBAL_RF + GLOBAL_XGB)/3)
  
  if(!return_models){
    return(list('all_predict' = df_result))
  } else {
      return(list('all_predict' = df_result,
              'gam_fit' = res_gam$model,
              'rf_fit' = res_rf$model,
              'xgb_fit' = res_xgb$model))
  }
}


###################################################################################################
################################## FONCTIONS MACRO ################################################
###################################################################################################

compute_scotland_forecasts <- function(dataset, windows, param, parallel_run=FALSE){

  if(parallel_run){
    n_cores <- parallel::detectCores() - 1
    message(paste("Lancement du modele national sur", length(windows), "fenetres..."))
      
    results <- parallel::mclapply(windows, function(w) {
      run_all_local_models(dataset, 
                           w,
                           param$gam_formula,
                           param$rf_formula,
                           param$mes_variables)$all_predict}, mc.cores = n_cores)

  }
  else{
    results <- lapply(windows, function(w) {
      run_all_local_models(dataset, 
                           w,
                           param$gam_formula,
                           param$rf_formula,
                           param$mes_variables)$all_predict})
  }
  
  return(dplyr::bind_rows(results, .id = "window_id"))
}


compute_regions_forecasts <- function(dataset, 
                                      windows, 
                                      local_param=NULL,
                                      global_param=NULL,
                                      parallel_run=FALSE){
  

  final_local_df <- NULL
  if(!is.null(local_param)){
    # --- I. MODELES LOCAUX ---
    message(paste("Modeles locaux par Region"))
    split_data <- dataset %>% split(.[['Region']])
    
    if(parallel_run){
      n_cores <- parallel::detectCores() - 1
      message(paste("Lancement pour", length(split_data), "series", n_cores, "coeurs..."))
      
      local_results <- parallel::mclapply(names(split_data), function(name) {
        entity_data <- split_data[[name]]
        purrr::map_dfr(windows, ~ run_all_local_models(entity_data, 
                                                       .x,
                                                       local_param$gam_formula_local, 
                                                       local_param$rf_formula_local,
                                                       local_param$mes_variables_local)$all_predict,
                      .id = "window_id")
      
    }, mc.cores = n_cores)
    }else{
      message(paste("Lancement pour", length(split_data), "series"))
      
      local_results <- lapply(names(split_data), function(name) {
        entity_data <- split_data[[name]]
        purrr::map_dfr(windows, ~ run_all_local_models(entity_data, 
                                                       .x,
                                                       local_param$gam_formula_local, 
                                                       local_param$rf_formula_local,
                                                       local_param$mes_variables_local)$all_predict,
                      .id = "window_id")})
    }
    
    # On réassigne les noms des séries avant de fusionner
    names(local_results) <- names(split_data)
    
    # Fusion finale de toutes les stations en un seul dataframe géant
    final_local_df <- dplyr::bind_rows(local_results, .id = 'Region')
  }

  final_global_df <- NULL
  if(!is.null(global_param)){

      # --- II. MODELES GLOBAUX --- 
      message(paste("Modeles globaux par Region"))
      num_threads = if(parallel_run) parallel::detectCores() - 1 else 1

      n_cores <- parallel::detectCores() - 1
      global_results <- lapply(windows, function(window) {
        run_all_global_models(dataset=dataset, 
                              df_window=window, 
                              gam_formula=global_param$gam_formula_global, 
                              rf_formula=global_param$rf_formula_global,
                              xgb_var=global_param$xgb_variables_global,
                              mes_variables=global_param$mes_variables_global,
                              target=global_param$target,
                              num_threads = num_threads,
                              max_pct_na=global_param$max_pct_na,
                              nbr_series_train = global_param$nbr_series_train,
                              SEED = 40)$all_predict})
      
      # On réassigne les noms des stations avant de fusionner
      final_global_df <- dplyr::bind_rows(global_results, .id = "window_id")
  }
  return(list('local' = final_local_df,
              'global' = final_global_df))
}


compute_stations_forecasts <- function(dataset, 
                                       windows, 
                                       output_path, 
                                       local_param = NULL,
                                       global_param = NULL,
                                       parallel_run = FALSE) {
  
  if(!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

  # --- I. MODÈLES LOCAUX (Un fichier par fenêtre) ---
  if(!is.null(local_param)) {
    message(">>> [LOCAL] Calcul par fenêtres et chunks...")
    
    station_ids <- unique(dataset$Station.ID)
    # Découpage en paquets de 100 stations
    station_chunks <- split(station_ids, ceiling(seq_along(station_ids) / 100))

    for(j in seq_along(windows)) {
      message(paste("   - Fenêtre", j, "/", length(windows)))
      
      for(k in seq_along(station_chunks)) {
        message(paste("     * Chunk station", k, "/", length(station_chunks)))
        ids_chunk <- station_chunks[[k]]
        
        if(parallel_run) {
          # Parallélisation par série au sein du chunk
          res_chunk <- parallel::mclapply(ids_chunk, function(id) {
            entity_data <- dataset[dataset$Station.ID == id, ]
            run_all_local_models(entity_data, windows[[j]], 
                                 local_param$gam_formula_local, 
                                 local_param$rf_formula_local, 
                                 local_param$mes_variables_local)$all_predict
          }, mc.cores = parallel::detectCores() - 1,
             mc.preschedule = FALSE) %>% 
            dplyr::bind_rows()
        } else {
          res_chunk <- purrr::map_dfr(ids_chunk, function(id) {
            entity_data <- dataset[dataset$Station.ID == id, ]
            run_all_local_models(entity_data, windows[[j]], 
                                 local_param$gam_formula_local, 
                                 local_param$rf_formula_local, 
                                 local_param$mes_variables_local)$all_predict
          })
        }

        # Sélection des colonnes et sauvegarde immédiate du chunk
        res_chunk <- res_chunk %>% 
          dplyr::select(Station.ID, Longitude, Latitude, Date, type, Consumed_kWh, starts_with("LOCAL_"))
        
        # Le nom du fichier inclut la fenêtre ET le chunk
        file_name <- paste0("local_win_", j, "_chunk_", k, ".RDS")
        saveRDS(res_chunk, file.path(output_path, file_name))
        
        # Nettoyage
        rm(res_chunk); gc()
      }
    }
  }

  # --- II. MODÈLES GLOBAUX --- 
  if(!is.null(global_param)) {
    message(">>> [GLOBAL] Calcul fenêtre par fenêtre...")

    num_threads = if(parallel_run) parallel::detectCores() - 1 else 1
    for(j in seq_along(windows)) {
      message(paste("   - Fenêtre globale", j, "/", length(windows)))
      res_win <- run_all_global_models(
        dataset = dataset, 
        df_window = windows[[j]], 
        gam_formula = global_param$gam_formula_global, 
        rf_formula = global_param$rf_formula_global,
        xgb_var = global_param$xgb_variables_global,
        mes_variables = global_param$mes_variables_global,
        target = global_param$target,
        num_threads = num_threads,
        max_pct_na=global_param$max_pct_na,
        nbr_series_train = global_param$nbr_series_train,
        SEED = 40
      )$all_predict %>%
        dplyr::select(Station.ID, Longitude, Latitude, Date, type, Consumed_kWh, starts_with("GLOBAL_"))
      
      if(global_param$target == "log_Consumed_kWh"){
        cols_to_fix <- c("GLOBAL_GAM", "GLOBAL_RF", "GLOBAL_XGB", "GLOBAL_MIX")
        res_win <- res_win %>% mutate(across(all_of(cols_to_fix), expm1))
      }

      saveRDS(res_win, file.path(output_path, paste0("global_win_", j, ".RDS")))
      rm(res_win); gc()
    } 
  }
}