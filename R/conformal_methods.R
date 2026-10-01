# -----------------------------
# CONFORMAL PREDICTION METHODS
# -----------------------------

library(dplyr)
library(parallel)

# -----------------------------
# HELPERS
# -----------------------------

compute_projected_predictions <- function(
    Yhat,
    idx,
    method,
    masks,
    H_t_list,
    ncores = 7,
    P = NULL
) {
  
  if (method == "Direct") {
    
    return(Yhat[idx, ])
    
  } else if (method == "OLS") {
    
    return((P %*% t(Yhat[idx, ])) %>% t())
    
  } else if (method == "refined OLS") {
    
    tmp <- parallel::mclapply(idx, function(t) {
      
      nb_active_t <- length(masks[[t]]) - sum(masks[[t]])
      
      P_t <- refined_OLS_projection_fast(
        H_t_list[[t]],
        max(500, nb_active_t)
      )
      
      P_t %*% Yhat[t, ]
      
    }, mc.cores = ncores)
    
    return(
      matrix(
        unlist(tmp),
        nrow = length(idx),
        byrow = TRUE
      )
    )
    
  } else {
    
    stop("Unknown projection method: ", method)
  }
}

compute_residuals <- function(
    Y,
    Yhat,
    idx,
    method,
    masks,
    H_t_list,
    ncores = 7,
    P = NULL
) {
  
  mu <- compute_projected_predictions(
    Yhat = Yhat,
    idx = idx,
    method = method,
    masks = masks,
    H_t_list = H_t_list,
    ncores = ncores,
    P = P
  )
  
  Y[idx, ] - mu
}

build_conformal_output <- function(
    covered,
    length,
    lower,
    upper,
    M_test,
    node_names
) {
  
  list(
    covered = covered,
    length = length,
    coverage_vec = colMeans(covered, na.rm = TRUE),
    coverage_all = mean(covered, na.rm = TRUE),
    length_vec = colMeans(length, na.rm = TRUE),
    length_all = mean(length, na.rm = TRUE),
    lower = lower,
    upper = upper,
    n_eff = nrow(M_test) - colSums(M_test, na.rm = TRUE),
    node_names = node_names 
  )
}

res_to_df <- function(
    res,
    method_name,
    window_id,
    H
) {
  
  node_names <- res$node_names
  
  station_nodes <- colnames(H)[-1]
  # regional_nodes <- grep(
  #   "^Region_",
  #   station_nodes,
  #   value = TRUE
  # )
  # station_only_nodes <- grep(
  #   "^Station_",
  #   station_nodes,
  #   value = TRUE
  # )
  
  data.frame(
    node = node_names,
    coverage = as.numeric(res$coverage_vec),
    length = as.numeric(res$length_vec),
    method = method_name,
    n_eff = as.numeric(res$n_eff)
  ) %>%
    mutate(
      window = window_id,
      node_type = case_when(
        node == "Y_nat" ~ "national",
        # node %in% regional_nodes ~ "regional",
        node %in% 1:(ncol(H) - 1) ~ "station", 
        TRUE ~ "regional"
        # node %in% station_only_nodes ~ "station",
        # TRUE ~ "unknown"
      )
    )
}

res_to_df_window <- function(
    res,
    method_name,
    window_id,
    H
) {
  
  data.frame(
    node = res$node_names,
    coverage = as.numeric(res$coverage_vec),
    length = as.numeric(res$length_vec),
    method = method_name,
    n_eff = as.numeric(res$n_eff)
  ) %>%
    mutate(
      window = window_id,
      node_type = case_when(
        node == "Y_nat" ~ "national",
        node %in% 1:(ncol(H)-1) ~ "station",
        TRUE ~ "regional"
      )
    )
}

res_to_df_time <- function(
    res,
    method_name,
    window_id,
    test_dates
) {
  
  expand.grid(
    date = test_dates,
    node = res$node_names
  ) %>%
    mutate(
      covered = as.vector(res$covered),
      length  = as.vector(res$length),
      lower   = as.vector(res$lower),
      upper   = as.vector(res$upper),
      method  = method_name,
      window  = window_id
    )
}

# -----------------------------
# CP_MNR
# -----------------------------

CP_MNR <- function(
    D_calib,
    D_test,
    Y,
    Yhat,
    masks,
    H_t_list,
    alpha = 0.1,
    method = "Direct",
    ncores = 7,
    P = NULL
) {
  node_names = colnames(Y)
  
  residuals <- compute_residuals(
    Y = Y,
    Yhat = Yhat,
    idx = D_calib,
    method = method,
    masks = masks,
    H_t_list = H_t_list,
    ncores = ncores,
    P = P
  )
  # --- The following step can be used as a safety operation for CP_MNR but is not compatible with CP_MNR_naive
  # # --- masking calibration residuals ---
  # 
  # for (k in seq_along(D_calib)) {
  #   
  #   t <- D_calib[k]
  #   
  #   residuals[k, masks[[t]] == 1] <- NA
  # }
  
  scores <- abs(residuals)
  
  q_vec <- apply(
    scores,
    2,
    quantile,
    probs = 1 - alpha,
    type = 1,
    na.rm = TRUE
  )
  
  mu_test <- compute_projected_predictions(
    Yhat = Yhat,
    idx = D_test,
    method = method,
    masks = masks,
    H_t_list = H_t_list,
    ncores = ncores,
    P = P
  )
  
  lower <- sweep(mu_test, 2, q_vec, "-")
  upper <- sweep(mu_test, 2, q_vec, "+")
  
  M_test <- matrix(
    unlist(masks[D_test]),
    nrow = length(D_test),
    byrow = TRUE
  )
  
  covered <- (Y[D_test, ] >= lower) &
    (Y[D_test, ] <= upper)
  
  covered[M_test == 1] <- NA
  
  length <- upper - lower
  length[M_test == 1] <- NA
  
  build_conformal_output(
    covered = covered,
    length = length,
    lower = lower,
    upper = upper,
    M_test = M_test,
    node_names = node_names 
  )
}

# -----------------------------
# CP_MNR_NAIVE
# -----------------------------

CP_MNR_naive <- function(
    D_calib,
    D_test,
    Y,
    Yhat,
    masks,
    H_t_list,
    alpha = 0.1,
    method = "Direct",
    ncores = 7,
    P = NULL
) {
  node_names = colnames(Y)
  
  ignore <- is.na(Y[, 1])
  
  Y_imp <- Y
  
  Y_imp[!ignore & is.na(as.matrix(Y_imp))] <- 0
  
  CP_MNR(
    D_calib = D_calib,
    D_test = D_test,
    Y = Y_imp,
    Yhat = Yhat,
    masks = masks,
    H_t_list = H_t_list,
    alpha = alpha,
    method = method,
    ncores = ncores,
    P = P
  )
}

# -----------------------------
# CP_MNR_NESTED_STAR
# -----------------------------
CP_MNR_Nested_star <- function(D_calib, 
                               D_test, 
                               Y, 
                               Yhat, 
                               masks,
                               H_t_list, 
                               alpha = 0.1,
                               ncores = 14, 
                               P = NULL,
                               max_diff = 600){
  node_names = colnames(Y)
  
  # -----------------------------
  # OUTPUT MATRICES
  # -----------------------------
  n_test <- length(D_test)
  # --- test mask ---
  M_mat <- matrix(unlist(masks[c(D_calib,D_test)]),
                  nrow = length(D_calib)+length(D_test),
                  byrow = TRUE)
  m <- ncol(M_mat)
  
  lower_nest <- matrix(NA, n_test, m)
  upper_nest <- matrix(NA, n_test, m)
  
  # -----------------------------
  # MAIN LOOP (nested structure)
  # -----------------------------
  for (i in seq_along(D_test)) {
    t <- D_test[i]
    print(t)
    
    # --- similarity filtering ---
    diff_mask <- rowSums(
      sweep(M_mat[D_calib, ], 2, M_mat[t, ], FUN = "!=")
    )
    idx_calib_filtered <- D_calib[diff_mask <= max_diff]
    
    if (length(idx_calib_filtered) == 0) next
    
    M_calib_t <- pmax(M_mat[idx_calib_filtered, , drop = FALSE],
      matrix(M_mat[t, ], nrow = length(idx_calib_filtered), ncol = m, byrow = TRUE))
    # -----------------------------
    # SCORES (UNIFIED PROJECTION)
    # -----------------------------
    n_calib <- length(idx_calib_filtered)
    
    proj_list <- mclapply(seq_len(n_calib), function(ii) {
      
      i <- idx_calib_filtered[ii]
      
      M_i <- M_calib_t[ii, ]
      sum(M_i)
      H_i <- mask_structural(H, M_i)
      nb_active_t <- length(M_i) - sum(M_i)
      
      refined_OLS_projection_fast(H_i, max(500,nb_active_t))
    }, mc.cores = ncores)
    
    Scores <- t(vapply(seq_len(n_calib), function(ii) {
      
      i <- idx_calib_filtered[ii]
      
      P <- proj_list[[ii]]
      
      # 
      y_i <- Y[i, ]
      y_i[is.na(y_i)] <- 0
      
      # score <- P %*% (y_i - Yhat[i, ]) # The projection may apply only to the forecast
      score <-  (y_i - P %*% Yhat[i, ]) 
      score[score == 0] <- NA
      
      as.numeric(score)
      
    }, numeric(m)))
    
    mu_t_proj <- t(vapply(seq_len(n_calib), function(ii) {
      
      P <- proj_list[[ii]]
      
      as.numeric(P %*% Yhat[t, ])
      
    }, numeric(m)))
    
    Scores_up   <- mu_t_proj + abs(Scores)
    Scores_down <- mu_t_proj - abs(Scores)
    
    # It should not change anything 
    Scores_up[M_calib_t == 1] <- NA 
    Scores_down[M_calib_t == 1] <- NA
    
    # -----------------------------
    # NESTED INTERVAL FUNCTION
    # -----------------------------
    nested_component <- function(sd, su, alpha) {
      
      valid <- which(!is.na(sd) & !is.na(su))
      if (length(valid) == 0) return(c(NA, NA))
      
      sd <- sd[valid]; su <- su[valid]
      n <- length(sd)
      
      # Run the sweep line algorithm
      x <- c(sd, su)
      d <- c(rep(1, n), rep(-1, n))
      
      ord <- order(x, -d)
      x <- x[ord]; d <- d[ord]
      
      csum <- cumsum(d)
      inside <- csum >= floor(alpha * (n + 1)) #+ 1) # I believe that the +1 was a computation mistake
      
      inside <- inside | c(inside[-1], FALSE) # Consider the closure 'manually'
       
      if (!any(inside)) return(c(NA, NA))
      
      c(lower = min(x[inside]),
        upper = max(x[inside]))
    }
    
    # -----------------------------
    # FINAL INTERVALS
    # -----------------------------
    result <- lapply(seq_len(m), function(j) {
      nested_component(Scores_down[, j],
                       Scores_up[, j],
                       alpha)
    })
    
    lower_nest[i, ] <- vapply(result, function(x) {
      if (is.null(x) || length(x) == 0 || !is.numeric(x)) return(NA_real_)
      x[["lower"]]
    }, numeric(1))
    
    upper_nest[i, ] <- vapply(result, function(x) {
      if (is.null(x) || length(x) == 0 || !is.numeric(x)) return(NA_real_)
      x[["upper"]]
    }, numeric(1))
  }
  
  # -----------------------------
  # EVALUATION
  # -----------------------------
  y_test <- Y[D_test, ]
  M_test <- M_mat[D_test, ]
  
  covered <- (y_test >= lower_nest) & (y_test <= upper_nest)
  covered[M_test == 1] <- NA
  
  lengths <- upper_nest - lower_nest
  lengths[M_test == 1] <- NA
  
  # list(
  #   lower_nest = lower_nest,
  #   upper_nest = upper_nest,
  #   covered = covered,
  #   coverage_vec = colMeans(covered, na.rm = TRUE),
  #   coverage_all = mean(covered, na.rm = TRUE),
  #   length_vec = colMeans(lengths, na.rm = TRUE),
  #   length_all = mean(lengths, na.rm = TRUE),
  #   n_eff = n_test - colSums(M_test, na.rm = TRUE),
  #   node_names = node_names
  # )
  build_conformal_output(
    covered = covered,
    length = lengths,
    lower = lower_nest,
    upper = upper_nest,
    M_test = M_test,
    node_names = node_names
  )
}
# -----------------------------
# ADAPTIVE_CP_MNR
# -----------------------------

Adaptive_CP_MNR <- function(
    D_calib,
    D_test,
    Y,
    Yhat,
    masks,
    H_t_list,
    alpha = 0.1,
    method = "Direct",
    ncores = 7,
    P = NULL,
    gamma = 0.02,
    a_0,
    alpha_0
) {
  node_names = colnames(Y)
  
  residuals <- compute_residuals(
    Y = Y,
    Yhat = Yhat,
    idx = D_calib,
    method = method,
    masks = masks,
    H_t_list = H_t_list,
    ncores = ncores,
    P = P
  )
  
  for (k in seq_along(D_calib)) {
    
    t <- D_calib[k]
    
    residuals[k, masks[[t]] == 1] <- NA
  }
  
  scores <- abs(residuals)
  
  mu_test <- compute_projected_predictions(
    Yhat = Yhat,
    idx = D_test,
    method = method,
    masks = masks,
    H_t_list = H_t_list,
    ncores = ncores,
    P = P
  )
  
  M_test <- matrix(
    unlist(masks[D_test]),
    nrow = length(D_test),
    byrow = TRUE
  )
  
  n_test <- length(D_test)
  n_calib <- length(D_calib)
  
  m <- ncol(Y)
  
  covered <- matrix(NA, n_test, m)
  length <- matrix(NA, n_test, m)
  
  lower <- matrix(NA, n_test, m)
  upper <- matrix(NA, n_test, m)
  
  a_t <- a_0
  alpha_t <- alpha_0
  
  threshold_vec <- function(s, a_vec) {
    
    mapply(function(col, a) {
      
      if (a == 1) {
        
        0
        
      } else if (a == 0) {
        
        2 * quantile(
          col,
          probs = 1,
          type = 1,
          na.rm = TRUE
        )
        
      } else {
        
        quantile(
          col,
          probs = 1 - a,
          type = 1,
          na.rm = TRUE
        )
      }
      
    }, as.data.frame(s), a_vec)
  }
  
  for (j in seq_along(D_test)) {
    t <- D_test[j]
    
    q_vec <- threshold_vec(scores, a_t)
    
    lower[j, ] <- mu_test[j, ] - q_vec
    upper[j, ] <- mu_test[j, ] + q_vec
    
    covered[j, ] <- (
      Y[t, ] >= lower[j, ]
    ) & (
      Y[t, ] <= upper[j, ]
    )
    
    covered[j, M_test[j, ] == 1] <- NA
    
    length[j, ] <- upper[j, ] - lower[j, ]
    length[j, M_test[j, ] == 1] <- NA
    
    err_t <- abs(Y[t, ] - mu_test[j, ])
    
    alpha_t <- alpha_t - gamma * (
      ifelse(
        is.na(err_t) | is.na(q_vec),
        alpha, 
        ifelse(
          err_t <= q_vec,
          0,
          1
        )
      ) - alpha
    )
    
    rounding <- function(alpha, n_calib) {
      
      1 - ceiling(
        (n_calib + 1) * (1 - alpha)
      ) / (n_calib + 1)
    }
    
    a_t <- ifelse(
      alpha_t >= 0 & alpha_t <= 1,
      rounding(alpha_t, n_calib),
      pmin(pmax(alpha_t, 0), 1)
    )
  }
  
  out <- build_conformal_output(
    covered = covered,
    length = length,
    lower = lower,
    upper = upper,
    M_test = M_test,
    node_names = node_names
  )
  
  out$a_final <- a_t
  out$alpha_final <- alpha_t
  
  out
}

# -----------------------------
# METHOD DISPATCHER
# -----------------------------

run_conformal_method <- function(
    conformal_method,
    ...
) {
  
  switch(
    
    conformal_method,
    
    "CP_MNR" = CP_MNR(...),
    
    "CP_MNR_naive" = CP_MNR_naive(...),
    
    "Adaptive_CP_MNR" = Adaptive_CP_MNR(...),
    
    "CP_MNR_Nested_star" = CP_MNR_Nested_star(...),
    
    stop("Unknown conformal method")
  )
}

# -----------------------------
# RUN WINDOW
# -----------------------------

run_window_CP <- function(
    window,
    i,
    res_nat,
    res_reg,
    res_sta,
    H,
    P_OLS,
    model_names,
    conformal_method,
    projection_methods,
    ncores = 7,
    alpha = 0.1,
    adaptive_state = NULL
) {
  
  print(i)
  
  test_period <- window %>%
    filter(type == "test") %>%
    pull(Date)
  
  calib_period <- window %>%
    filter(type == "calibration") %>%
    pull(Date)
  
  vectors_test <- build_vectors(
    period = test_period,
    res_nat = res_nat %>% filter(window_id == i),
    res_reg = res_reg %>% filter(window_id == i),
    res_sta = res_sta %>% filter(window_id == i),
    model_names = model_names
  )
  
  vectors_calib <- build_vectors(
    period = calib_period,
    res_nat = res_nat %>% filter(window_id == i),
    res_reg = res_reg %>% filter(window_id == i),
    res_sta = res_sta %>% filter(window_id == i),
    model_names = model_names
  )
  
  Y <- do.call(
    rbind,
    lapply(vectors_calib, `[[`, "y")
  )
  
  Yhat <- do.call(
    rbind,
    lapply(vectors_calib, `[[`, "yhat")
  )
  
  Y <- rbind(
    Y,
    do.call(rbind, lapply(vectors_test, `[[`, "y"))
  )
  
  Yhat <- rbind(
    Yhat,
    do.call(rbind, lapply(vectors_test, `[[`, "yhat"))
  )
  
  D_calib <- seq_len(length(calib_period))
  
  D_test <- (
    length(calib_period) + 1
  ):(
    length(calib_period) + length(test_period)
  )
  
  M <- is.na(Y) * 1
  
  masks <- split(M, seq_len(nrow(M)))
  
  H_t_list <- lapply(
    masks,
    function(M_t) {
      mask_structural(H, M_t)
    }
  )
  
  results_df <- list()
  
  updated_state <- adaptive_state
  results_window <- list()
  results_time <- list()
  
  # -----------------------------
  # SPECIAL CASE:
  # CP_MNR_Nested_star
  # -----------------------------
  if (conformal_method == "CP_MNR_Nested_star") {
    
    res <- run_conformal_method(
      conformal_method = conformal_method,
      D_calib = D_calib,
      D_test = D_test,
      Y = Y,
      Yhat = Yhat,
      masks = masks,
      H_t_list = H_t_list,
      alpha = alpha,
      ncores = ncores
    )
    
    summary_df <- res_to_df_window(
      res = res,
      method_name = "Nested",
      window_id = i,
      H = H
    )
    
    time_df <- res_to_df_time(
      res = res,
      method_name = "Nested",
      window_id = i,
      test_dates = test_period
    )
    
    return(list(
      summary_df = summary_df,
      time_df = time_df,
      state = updated_state
    ))
  }
  # if (conformal_method == "CP_MNR_Nested_star") {
  #   
  #   res <- run_conformal_method(
  #     conformal_method = conformal_method,
  #     D_calib = D_calib,
  #     D_test = D_test,
  #     Y = Y,
  #     Yhat = Yhat,
  #     masks = masks,
  #     H_t_list = H_t_list,
  #     alpha = alpha,
  #     ncores = ncores
  #     )
  #   
  #   results_df[["Nested"]] <- res_to_df(
  #     res = res,
  #     method_name = "Nested",
  #     window_id = i,
  #     H = H
  #   )
  #   
  #   return(list(
  #     df = bind_rows(results_df),
  #     state = updated_state
  #   ))
  # }
  
  # -----------------------------
  # STANDARD METHODS
  # -----------------------------
  
  for (method in projection_methods) {
    
    if (conformal_method == "Adaptive_CP_MNR") {
      
      res <- run_conformal_method(
        conformal_method = conformal_method,
        D_calib = D_calib,
        D_test = D_test,
        Y = Y,
        Yhat = Yhat,
        masks = masks,
        H_t_list = H_t_list,
        alpha = alpha,
        method = method,
        ncores = ncores,
        P = P_OLS,
        a_0 = adaptive_state[[method]]$a_t,
        alpha_0 = adaptive_state[[method]]$alpha_t
      )
      
      updated_state[[method]]$a_t <- res$a_final
      
      updated_state[[method]]$alpha_t <- res$alpha_final
      
    } else {
      
      res <- run_conformal_method(
        conformal_method = conformal_method,
        D_calib = D_calib,
        D_test = D_test,
        Y = Y,
        Yhat = Yhat,
        masks = masks,
        H_t_list = H_t_list,
        alpha = alpha,
        method = method,
        ncores = ncores,
        P = P_OLS
      )
    }
    
    # results_df[[method]] <- res_to_df(
    #   res = res,
    #   method_name = method,
    #   window_id = i,
    #   H = H
    # )
    results_window[[method]] <- res_to_df_window(
      res,
      method,
      i,
      H
    )
    
    results_time[[method]] <- res_to_df_time(
      res,
      method,
      i,
      test_period
    )
  }
  
  list(
    summary_df = bind_rows(results_window),
    time_df = bind_rows(results_time),
    state = updated_state
  )
}

# -----------------------------
# RUN EXPERIMENT
# -----------------------------

run_conformal_experiment <- function(
    windows,
    res_nat,
    res_reg,
    res_sta,
    H,
    P_OLS,
    model_names,
    conformal_method,
    projection_methods,
    alpha = 0.1,
    ncores = 7
) {
  
  adaptive_state <- NULL
  
  if (conformal_method == "Adaptive_CP_MNR") {
    
    m <- ncol(P_OLS)
    
    init_state <- function() {
      
      list(
        a_t = rep(alpha, m),
        alpha_t = rep(alpha, m)
      )
    }
    
    adaptive_state <- setNames(
      lapply(projection_methods, function(x) init_state()),
      projection_methods
    )
  }
  
  # all_results <- lapply(seq_along(windows), function(i) {
  #   
  #   out <- run_window_CP(
  #     window = windows[[i]],
  #     i = i,
  #     res_nat = res_nat,
  #     res_reg = res_reg,
  #     res_sta = res_sta,
  #     H = H,
  #     P_OLS = P_OLS,
  #     model_names = model_names,
  #     conformal_method = conformal_method,
  #     projection_methods = projection_methods,
  #     ncores = ncores,
  #     alpha = alpha,
  #     adaptive_state = adaptive_state
  #   )
  #   
  #   if (conformal_method == "Adaptive_CP_MNR") {
  #     
  #     adaptive_state <<- out$state
  #   }
  #   
  #   out$df
  # })
  # 
  # bind_rows(all_results)
  all_summary <- list()
  all_time <- list()
  
  for (i in seq_along(windows)) {
    
    out <- run_window_CP(
      window = windows[[i]],
      i = i,
      res_nat = res_nat,
      res_reg = res_reg,
      res_sta = res_sta,
      H = H,
      P_OLS = P_OLS,
      model_names = model_names,
      conformal_method = conformal_method,
      projection_methods = projection_methods,
      ncores = ncores,
      alpha = alpha,
      adaptive_state = adaptive_state
    )
    
    all_summary[[i]] <- out$summary_df
    all_time[[i]] <- out$time_df
    
    if (conformal_method == "Adaptive_CP_MNR") {
      adaptive_state <- out$state
    }
  }
  
  list(
    summary = bind_rows(all_summary),
    time = bind_rows(all_time)
  )
}