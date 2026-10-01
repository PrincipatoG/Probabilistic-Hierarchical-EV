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
  
  # 1. Initial join
  df_full <- df_window %>% left_join(dataset, by = 'Date')
  
  # --- METHOD 1: ARIMA ---
  res_arima <- tryCatch({
    train_vec_clean <- df_full$Consumed_kWh_SARIMA[df_full$type %in% c('train', 'calibration')]
    fit_train <- auto.arima(ts(train_vec_clean, frequency = 7))
    
    # Apply to the full series to preserve continuity
    fit_applied <- Arima(df_full$Consumed_kWh_SARIMA, model = fit_train)
    full_preds <- as.numeric(fitted(fit_applied))
    
    # Filter as needed
    pred_out <- if(!return_models) full_preds[df_full$type != "train"] else full_preds
    
    list(model = fit_train, pred = pred_out, status = "success")
  }, error = function(e){
    list(model=NULL, pred = rep(NA, if(!return_models) sum(df_full$type != "train") else nrow(df_full)), status = "error")
  }) 

  # --- METHOD 2: GAM ---
  res_gam <- tryCatch({
    m_gam <- gam(as.formula(gam_formula), 
                 data = df_full %>% filter(type == "train") %>% drop_na(Consumed_kWh, all_of(mes_variables)))
    
    # Prediction
    data_to_pred <- if(!return_models) df_full %>% filter(type != "train") else df_full
    list(model=m_gam, pred = as.numeric(predict(m_gam, newdata = data_to_pred)), status = "success")
  }, error = function(e){
    list(model=NULL, pred = rep(NA, if(!return_models) sum(df_full$type != "train") else nrow(df_full)), status = "error")
  }) 

  # --- METHOD 3: RANDOM FOREST (RANGER) ---
  res_rf <- tryCatch({
    rf_fit <- ranger(as.formula(rf_formula), 
                     data = df_full %>% filter(type == "train") %>% drop_na(Consumed_kWh, all_of(mes_variables)), 
                     num.trees = 500, 
                     num.threads = 1,
                     importance = 'impurity')
    
    # On-the-fly prediction on the subset
    data_to_pred <- if(!return_models) df_full %>% filter(type != "train") else df_full
    list(model=rf_fit, pred = predict(rf_fit, data = data_to_pred)$predictions, status = "success")
  }, error = function(e) {
    list(model=NULL, pred = rep(NA, if(!return_models) sum(df_full$type != "train") else nrow(df_full)), status = "error")
  })

  # --- FINAL ASSEMBLY ---
  # Build the final data frame by filtering df_full when return_models is FALSE
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

  # --- SAMPLING BEFORE RUNNING FOR THE STATIONS (OPTIONAL) --- #
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

  # --- Define the prediction set ---
  df_pred <- if(!return_models) df_full[df_full$type != "train", ] else df_full

  # --- METHOD 1 : GAM ---
  m_gam <- NULL
  res_gam <- tryCatch({
    m_gam <- bam(as.formula(gam_formula), 
                   data = df_full[df_full$type_use == "train", c(target, mes_variables)] %>% drop_na(),
                   discrete = TRUE,
                   nthreads = num_threads)
    list(model = m_gam, pred = as.numeric(predict(m_gam, newdata = df_pred)), status = "success") 
  }, error = function(e){
    message("GAM ERROR detected: ", e$message)
    list(model = NULL, pred = rep(NA, nrow(df_pred)), status = "error")}) 
  gc()

  # --- METHOD 2 : RANDOM FOREST ---
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
    message("RF ERROR detected: ", e$message)
    list(model = NULL, pred = rep(NA, nrow(df_pred)), status = "error")})
  gc()

  # --- METHOD 3 : XGBOOST ----
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
    message("XGBoost ERROR detected: ", e$message); list(model=NULL, pred = rep(NA, nrow(df_pred)), status = "error")
  })

  # --- FINAL OUTPUT ---
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
################################## MACRO FUNCTION ################################################
###################################################################################################

compute_scotland_forecasts <- function(dataset, windows, param, parallel_run=FALSE){

  if(parallel_run){
    n_cores <- parallel::detectCores() - 1
    message(paste("Launching the national model for", length(windows), "windows..."))
      
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
    # --- I. LOCAL MODELS ---
    message(paste("Local models by region"))
    split_data <- dataset %>% split(.[['Region']])
    
    if(parallel_run){
      n_cores <- parallel::detectCores() - 1
      message(paste("Launching", length(split_data), "series on", n_cores, "cores..."))
      
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
      message(paste("Launching", length(split_data), "series"))
      
      local_results <- lapply(names(split_data), function(name) {
        entity_data <- split_data[[name]]
        purrr::map_dfr(windows, ~ run_all_local_models(entity_data, 
                                                       .x,
                                                       local_param$gam_formula_local, 
                                                       local_param$rf_formula_local,
                                                       local_param$mes_variables_local)$all_predict,
                      .id = "window_id")})
    }
    
    # Reassign series names before merging
    names(local_results) <- names(split_data)
    
    # Final merge of all stations into one large data frame
    final_local_df <- dplyr::bind_rows(local_results, .id = 'Region')
  }

  final_global_df <- NULL
  if(!is.null(global_param)){

      # --- II. GLOBAL MODELS --- 
      message(paste("Global models by region"))
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
      
      # Reassign station names before merging
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

  # --- I. LOCAL MODELS (one file per window) ---
  if(!is.null(local_param)) {
    message(">>> [LOCAL] Computing windows and chunks...")
    
    station_ids <- unique(dataset$Station.ID)
    # Split into batches of 100 stations
    station_chunks <- split(station_ids, ceiling(seq_along(station_ids) / 100))

    for(j in seq_along(windows)) {
      message(paste("   - Window", j, "/", length(windows)))
      
      for(k in seq_along(station_chunks)) {
        message(paste("     * Station chunk", k, "/", length(station_chunks)))
        ids_chunk <- station_chunks[[k]]
        
        if(parallel_run) {
          # Parallelize by series within the chunk
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

        # Select columns and save the chunk immediately
        res_chunk <- res_chunk %>% 
          dplyr::select(Station.ID, Longitude, Latitude, Date, type, Consumed_kWh, starts_with("LOCAL_"))
        
        # The filename includes both the window and the chunk
        file_name <- paste0("local_win_", j, "_chunk_", k, ".RDS")
        saveRDS(res_chunk, file.path(output_path, file_name))
        
        # Cleanup to free memory
        rm(res_chunk); gc()
      }
    }
  }

  # --- II. GLOBAL MODELS ---
  if(!is.null(global_param)) {
    message(">>> [GLOBAL] Computing window by window...")

    num_threads = if(parallel_run) parallel::detectCores() - 1 else 1
    for(j in seq_along(windows)) {
      message(paste("   - Global window", j, "/", length(windows)))
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