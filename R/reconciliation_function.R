# Initialization functions
build_H <- function(res_nat, res_reg, res_sta) {
  station_levels <- res_sta %>% distinct(Station.ID,Region)
  n <- nrow(station_levels)
  
  ### --- Build structural matrix H --- ###
  # Scotland level: all CPIDs belong to Scotland
  H_scotland <- rep(1, n)
  
  # Region level: one-hot encoding of regions
  region_levels <- sort(unique(station_levels$Region))
  H_region <- map_dfc(region_levels, ~ as.integer(station_levels$Region == .x))
  colnames(H_region) <- paste0("Region_", region_levels)
  
  # Address level: one-hot encoding of Addresss
  Address_levels <- sort(unique(station_levels$Station.ID))
  H_Address <- map_dfc(Address_levels, ~ as.integer(station_levels$Station.ID == .x))
  colnames(H_Address) <- paste0("Station_", Address_levels)
  
  # Assemble final structural matrix (rows = CPID)
  H <- cbind(
    Stations = Address_levels,
    Scotland = H_scotland,
    H_region,
    H_Address
  )
  
  # Output the structural matrix
  return(as.data.frame(H))
}

# Extraction functions
get_mask <- function(date, res_nat, res_reg, res_sta, H){
  ### --- Scotland --- ###
  nat_val <- res_nat %>%
    filter(Date == date) %>%
    pull(Consumed_kWh)
  
  mask_nat <- as.integer(is.na(nat_val))
  
  ### --- Regions --- ###
  region_cols <- names(H) %>% grep("^Region_", ., value = TRUE)
  region_levels <- sub("Region_", "", region_cols)
  
  reg_vals <- res_reg %>%
    filter(Date == date) %>%
    dplyr::select(Region, Consumed_kWh)
  
  mask_reg <- sapply(region_levels, function(reg) {
    val <- reg_vals %>%
      filter(Region == reg) %>%
      pull(Consumed_kWh)
    
    as.integer(length(val) == 0 || is.na(val))
  })
  
  ### --- Stations --- ###
  station_cols <- names(H) %>% grep("^Station_", ., value = TRUE)
  station_levels <- sub("Station_", "", station_cols)
  
  sta_vals <- res_sta %>%
    filter(Date == date) %>%
    dplyr::select(Station.ID, Consumed_kWh)
  
  mask_sta <- sapply(station_levels, function(sta) {
    val <- sta_vals %>%
      filter(Station.ID == sta) %>%
      pull(Consumed_kWh)
    
    as.integer(length(val) == 0 || is.na(val))
  })
  
  ### --- Assemble --- ###
  mask <- c(mask_nat, mask_reg, mask_sta)
  print(sum(mask))
  return(mask)
  
}
mask_structural <- function(H,M_t){
  H_mat <- as.matrix(H[,-1])  # n x m
  active_cols <- (M_t == 0)
  
  # Invalid rows corresponding to inactive nodes
  impacted_rows <- as.vector(H_mat[, M_t == 1, drop = FALSE] %*% rep(1, sum(M_t == 1))) > 0
  
  H_t <- H_mat
  H_t[, M_t == 1] <- 0
  H_t[impacted_rows, ] <- 0
  
  return(H_t)
}
mask_structural_fast <- function(H, M_t){
  
  H_mat <- as.matrix(H[,-1])
  
  active_cols <- which(M_t == 0)
  
  impacted_rows <- rowSums(H_mat[, M_t == 1, drop = FALSE]) > 0
  
  keep_rows <- !impacted_rows
  
  H_red <- H_mat[keep_rows, active_cols, drop = FALSE]
  
  list(
    H_red = H_red,
    keep_rows = keep_rows,
    active_cols = active_cols,
    n = nrow(H_mat),
    m = ncol(H_mat)
  )
}
refined_OLS_projection <- function(H_t){
  t(H_t) %*% ginv(H_t %*% t(H_t)) %*% H_t 
}
refined_OLS_projection_fast <- function(H_t, nb_active = 800){
  
  # Leverage the sparsity of A 
  A <- tcrossprod(Matrix(H_t, sparse = TRUE))
  A <- as(A, "generalMatrix") # A <- as(A, "dgCMatrix") worked but was depreciated
  n <- nrow(A)
  
  # consider only the k higher eigen values
  k <- min(n - 1, nb_active)  # 800 is arbitrary, it should be the number of active nodes
  eg <- tryCatch(
    eigs_sym(A, k = k, which = "LM"),
    error = function(e) NULL
  )
  
  if (is.null(eg)) {
    # fallback robuste
    print("error in the sparse method")
    eg <- eigen(as.matrix(A), symmetric = TRUE)
  }  
  vals <- eg$values
  vecs <- eg$vectors
  
  # Moore-Penrose Pseudo-Inverse
  tol <- max(vals) * .Machine$double.eps
  keep <- vals > tol
  
  if (!any(keep)) {
    return(matrix(0, ncol(H_t), ncol(H_t)))
  }
  
  A_inv <- vecs[, keep, drop = FALSE] %*%
    diag(1 / vals[keep]) %*%
    t(vecs[, keep, drop = FALSE])
  
  # Final projection
  t(H_t) %*% A_inv %*% H_t
}

refined_WLS_projection <- function(H_t,inv_sigma_t){
  t(H_t) %*% ginv(H_t %*% inv_sigma_t %*% t(H_t)) %*% H_t %*% inv_sigma_t
}
# refined_OLS_projectionQR <- function(H_t){
#   qrH <- qr(t(H_t))    
#   Q <- qr.Q(qrH)
#   Q %*% t(Q)
# }
# refined_OLS_projectionSVD <- function(H) {
#   s <- svd(H)
#   U <- s$u
#   P <- U %*% t(U)
#   P
# }
# 
# 
# refined_OLS_projection_fast <- function(H, M_t){
#   
#   obj <- mask_structural_fast(H, M_t)
#   
#   H_red <- obj$H_red
#   n <- obj$n
#   
#   # QR decomposition on reduced system
#   qr_H <- qr(H_red)
#   
#   Q <- qr.Q(qr_H)  # n' x r orthonormal basis
#   
#   # projection in reduced space
#   P_red <- Q %*% t(Q)
#   
#   # embed back into full space
#   P_full <- matrix(0, n, n)
#   P_full[obj$keep_rows, obj$keep_rows] <- P_red
#   
#   return(P_full)
# }
masked_projection <- function(H_t,inv_Sigma_t){
  t(H_t) %*% ginv(H_t %*% inv_Sigma_t %*% t(H_t)) %*% H_t %*% inv_Sigma_t # Check that the pseudo-inverse does not produce nonsense
}
# Reconciliation functions
estimate_inv_sigma_bis <- function(res_nat, res_reg, res_sta){
    # NATIONAL
    err_nat <- res_nat %>%
      dplyr::filter(window_id == i) %>%
      dplyr::filter(Date %in% train_period) %>%
      dplyr::arrange(Date) %>%
      dplyr::transmute(err = Consumed_kWh - .data[[model_names[[1]]]]) %>%
      as.matrix()
    
    # REGIONS
    err_reg <- res_reg %>%
      dplyr::filter(window_id == i) %>%
      dplyr::filter(Date %in% train_period) %>%
      dplyr::mutate(err = Consumed_kWh - .data[[model_names[[2]]]]) %>%
      dplyr::select(Date, Region, err) %>%
      tidyr::pivot_wider(names_from = Region, values_from = err) %>%
      dplyr::arrange(Date) %>%
      dplyr::select(dplyr::all_of(region_levels)) %>%
      as.matrix()
    
    # STATIONS (spécifique fenêtre)
    err_sta <- res_sta %>%
      dplyr::filter(window_id == i) %>%
      dplyr::filter(Date %in% train_period) %>%
      dplyr::mutate(err = Consumed_kWh - .data[[model_names[[3]]]]) %>%
      dplyr::select(Date, Station.ID, err) %>%
      tidyr::pivot_wider(names_from = Station.ID, values_from = err) %>%
      dplyr::arrange(Date) %>%
      dplyr::select(dplyr::all_of(station_levels)) %>%
      as.matrix()
    
    # CONCAT
    E <- cbind(err_nat, err_reg, err_sta)
    
    mean((E[,30]-mean(E[,30], na.rm=T))^2, na.rm=T)
    # VARIANCE (ULTRA RAPIDE)
    var_estimates <- apply(E, 2, function(x) {
      x <- x[!is.na(x)]
      if (length(x) < 2) return(NA_real_)
      mean((x - mean(x))^2)
    })
    zero_var_idx <- which(var_estimates == 0)
    max(var_estimates, na.rm = T)
    
    # Fill NAs with large variance
    is_bad <- is.na(var_estimates)
    var_estimates[is_bad] <- 2 * max(var_estimates, na.rm = TRUE)
  
    inv_Sigma_t <- diag(1 / var_estimates)
    
  return(inv_Sigma_t) 
}
estimate_inv_sigma <- function(res_nat, res_reg, res_sta){
  
  # --------------------
  # NATIONAL
  # --------------------
  nat_df <- res_nat %>%
    dplyr::filter(window_id == i) %>%
    dplyr::filter(Date %in% train_period)
  
  err_nat <- nat_df %>%
    dplyr::filter(Date %in% train_period) %>%
    dplyr::arrange(Date) %>%
    dplyr::transmute(err = Consumed_kWh - .data[[model_names[[1]]]]) %>%
    as.matrix()
  
  # --------------------
  # REGIONS
  # --------------------
  err_reg <- res_reg %>%
    dplyr::filter(window_id == i) %>%
    dplyr::filter(Date %in% train_period) %>%
    dplyr::mutate(err = Consumed_kWh - .data[[model_names[[2]]]]) %>%
    dplyr::select(Date, Region, err) %>%
    tidyr::pivot_wider(names_from = Region, values_from = err) %>%
    dplyr::arrange(Date) %>%
    dplyr::select(dplyr::any_of(region_levels)) %>%   
    as.matrix()
  
  # --------------------
  # STATIONS
  # --------------------
  err_sta <- res_sta %>%
    dplyr::filter(window_id == i) %>%
    dplyr::filter(Date %in% train_period) %>%
    dplyr::mutate(err = Consumed_kWh - .data[[model_names[[3]]]]) %>%
    dplyr::select(Date, Station.ID, err) %>%
    tidyr::pivot_wider(names_from = Station.ID, values_from = err) %>%
    dplyr::arrange(Date) %>%
    dplyr::select(dplyr::any_of(station_levels)) %>%  
    as.matrix()
  
  # --------------------
  # CONCAT
  # --------------------
  E <- cbind(err_nat, err_reg, err_sta)
  
  # --------------------
  # VARIANCE
  # --------------------
  var_estimates <- apply(E, 2, function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 2) return(NA_real_)
    mean((x - mean(x))^2)
  })
  
  # --------------------
  # FALLBACK ROBUSTE
  # --------------------
  is_bad <- is.na(var_estimates)
  var_estimates[is_bad] <- Inf
  
  inv_Sigma_t <- diag(1 / var_estimates)
  
  return(inv_Sigma_t)
}
get_P <- function(H, inv_Sigma, H_list = NULL) {
  
  # Base matrix (without intercept column)
  H0 <- as.matrix(H[, -1])
  
  # -------------------------
  # 1. BASE PROJECTIONS
  # -------------------------
  P_OLS <- masked_projection(H0, diag(ncol(H0)))
  # P_WLS <- masked_projection(H0, inv_Sigma)
  
  # -------------------------
  # 2. REFINED PROJECTIONS
  # -------------------------
  if (!is.null(H_list)) {
    
    P_OLS_ref <- lapply(H_list, function(Ht) {
      masked_projection(Ht, diag(ncol(H0)))
    })
    
    # P_WLS_ref <- lapply(H_list, function(Ht) {
    #   masked_projection(Ht, inv_Sigma)
    # })
    # 
  } else {
    P_OLS_ref <- NULL
    # P_WLS_ref <- NULL
  }
  
  # -------------------------
  # 3. OUTPUT STRUCTURE
  # -------------------------
  list(
    OLS = P_OLS,
    OLS_ref = P_OLS_ref)
  # ,
  #   WLS = P_WLS,
  #   WLS_ref = P_WLS_ref
  # )
}
# Metric functions
gmse <- function(y, yhat) {
  sum((y - yhat)^2, na.rm = TRUE)
}
compute_metric <- function(y, yhat_direct, yhat_global, yhat_dyn) {
  
  list(
    GMSE_direct = sum((y - yhat_direct)^2, na.rm = TRUE),
    GMSE_global = sum((y - yhat_global)^2, na.rm = TRUE),
    GMSE_dyn    = sum((y - yhat_dyn)^2, na.rm = TRUE)
  )
}