# -----------------------------
# PACKAGES
# -----------------------------
library(readr)
library(stringr)
library(purrr)
library(fastmatrix)

# -----------------------------
# FUNCTIONS
# -----------------------------

mask_to_label <- function(mask_vec, var_names = paste0("Y", 1:length(mask_vec))) {
  n_missing <- sum(mask_vec == 0)
  if (all(mask_vec == 1)) return("Y fully observed")
  if (n_missing == n) return(NA) # erase the zero node hierarchy
  missing <- var_names[mask_vec == 0]
  paste0("(", paste(missing, collapse = ", "), ") missing")
}
# To adapt
compute_masked_coverage <- function(covered_mat, mask_mat) {
  # Unique mask
  mask_df <- as.data.frame(mask_mat)
  mask_df$idx <- seq_len(nrow(mask_df))
  unique_masks <- unique(mask_df[, -ncol(mask_df)])
  res <- list()
  # Loop on the mask
  for (i in seq_len(nrow(unique_masks))) {
    pattern <- unique_masks[i, , drop = FALSE]
    # test point on which this mask is applied
    matches <- which(apply(mask_mat, 1, function(x) all(x == unlist(pattern))))
    if (length(matches) == 0) next
    # component-wise coverage
    obs_cols <- which(pattern == 1)
    cov_vals <- covered_mat[matches, obs_cols, drop = FALSE]
    coverage_per_component <- colMeans(cov_vals, na.rm = TRUE)
    
    res[[paste0("mask_", i)]] <- list(
      pattern = pattern,
      indices = matches,
      coverage_vec = coverage_per_component,
      coverage_all = mean(coverage_per_component)
    )
  }
  
  res
}
compute_masked_length <- function(length_mat, mask_mat) {
  mask_df <- as.data.frame(mask_mat)
  mask_df$idx <- seq_len(nrow(mask_df))
  unique_masks <- unique(mask_df[, -ncol(mask_df)])
  res <- list()
  
  for (i in seq_len(nrow(unique_masks))) {
    pattern <- unique_masks[i, , drop = FALSE]
    matches <- which(apply(mask_mat, 1, function(x) all(x == unlist(pattern))))
    if (length(matches) == 0) next
    
    obs_cols <- which(pattern == 1)
    if(length(obs_cols) == 0) next  # skip if all values are missing
    
    len_vals <- length_mat[matches, obs_cols, drop = FALSE]
    
    # Exactly as in compute_masked_coverage
    res[[paste0("mask_", i)]] <- list(
      pattern = pattern,
      indices = matches,
      length_vec = colMeans(len_vals, na.rm = TRUE),  # Mean for each observed column
      length_all = mean(len_vals, na.rm = TRUE)       # Overall mean across observed columns
    )
  }
  
  res
}

CP_MNR <- function(D_calib, D_test, Y, Yhat, M_mat, H_t_list, alpha = 0.1){
  # Calibration scores
  residuals <- Y[D_calib,] - Yhat[D_calib,]
  # Here, we implement "refined OLS" (which is computationnally extensive = 10 minutes but here we only do it ones)
  coherent_residuals <- t(sapply(D_calib, function(t) {
    refined_OLS_projection(H_t_list[[t]]) %*% (Y[t, ] - Yhat[t, ]) # check the dimension of the matrices/vectors 
  }))
  
  # Mask the masked components in the scores
  M_calib <- M_mat[D_calib, ]
  all_children_missing <- rowSums(M_calib[, 2:(m), drop = FALSE]) == 0
  M_calib[all_children_missing, 1] <- 0
  residuals[M_calib == 0] <- NA
  coherent_residuals[M_calib == 0] <- NA
  
  scores <- abs(residuals)
  coherent_scores <- abs(coherent_residuals)  
  q_vec <- apply(scores, 2, function(s) quantile(s, probs = 1 - alpha, type = 1, na.rm = TRUE))
  q_vec_coherent <- apply(coherent_scores, 2, function(s) quantile(s, probs = 1 - alpha, type = 1, na.rm = TRUE))
  
  # Test prediction sets
  mu_test <- Yhat[D_test,]
  
  lower <- sweep(mu_test, 2, q_vec, `-`)
  upper <- sweep(mu_test, 2, q_vec, `+`)
  
  mu_test_coherent <- t(sapply(D_test, function(t) {
    masked_projection(H_array[t,,]) %*% (Yhat[t, ])
  }))
  lower_coherent <- sweep(mu_test_coherent, 2, q_vec_coherent, `-`)
  upper_coherent <- sweep(mu_test_coherent, 2, q_vec_coherent, `+`)
  
  # Empirical coverages
  y_test <- Y[D_test, ]
  M_test <- M_mat[D_test, ]
  
  covered <- (y_test >= lower) & (y_test <= upper)
  covered[M_test == 0] <- NA 
  length <- (upper - lower)
  length[M_test == 0] <- NA 
  coverage_vec <- colMeans(covered, na.rm = TRUE)
  coverage_all <- mean(covered, na.rm = TRUE)
  length_vec   <- colMeans(length, na.rm =TRUE)
  length_all   <- mean(length_vec)
  
  covered_coherent <- (y_test >= lower_coherent) & (y_test <= upper_coherent)
  covered_coherent[M_test == 0] <- NA 
  length_coherent <- (upper_coherent - lower_coherent)
  length_coherent[M_test == 0] <- NA 
  coverage_vec_coherent <- colMeans(covered_coherent, na.rm = TRUE)
  coverage_all_coherent <- mean(covered_coherent, na.rm = TRUE)
  length_vec_coherent   <- colMeans(length_coherent, na.rm = TRUE)
  length_all_coherent   <- mean(length_vec_coherent)
  
  # Mask conditional coverage
  mask_conditional_coverage <- compute_masked_coverage(covered, M_test)
  mask_conditional_coverage_coherent <- compute_masked_coverage(covered_coherent, M_test)
  mask_conditional_length <- compute_masked_length(length, M_test)
  mask_conditional_length_coherent <- compute_masked_length(length_coherent, M_test)
  
  
  list(
    q_vec = q_vec,
    q_vec_coherent = q_vec_coherent,
    lower = lower,
    upper = upper,
    lower_coherent = lower_coherent,
    upper_coherent = upper_coherent,
    covered = covered,
    coverage_vec = coverage_vec,
    coverage_all = coverage_all,
    length_vec = length_vec,
    length_all = length_all,
    covered_coherent = covered_coherent,
    coverage_vec_coherent = coverage_vec_coherent,
    coverage_all_coherent = coverage_all_coherent,
    length_vec_coherent = length_vec_coherent,
    length_all_coherent = length_all_coherent,
    mask_conditional = mask_conditional_coverage,
    mask_conditional_coherent = mask_conditional_coverage_coherent,
    mask_conditional_length = mask_conditional_length,
    mask_conditional_length_coherent = mask_conditional_length_coherent
  )
}

CP_MNR_Nested <- function(D_calib, D_test, Y, Yhat, M_mat, H_array, alpha = 0.1){
  # Initialize the matrices to collect the prediction sets
  n_test <- length(D_test)
  m <- ncol(M_mat)
  lower_nest <- matrix(NA, nrow = n_test, ncol = m)
  upper_nest <- matrix(NA, nrow = n_test, ncol = m)
  # For each test point, we compute new scores
  for (i in seq_along(D_test)) {
    t <- D_test[i]
    print(i)
    # Combine the mask
    M_t <- M_mat[t, ]
    M_calib_t <- 1*(M_mat[D_calib, ] & matrix(M_t, nrow = length(D_calib), ncol = m, byrow = TRUE))
    # Manual fix to deal with the zero node case
    all_children_missing <- rowSums(M_calib_t[, 2:(m), drop = FALSE]) == 0
    M_calib_t[all_children_missing, 1] <- 0
    # Scores MNR Nested
    Scores <- t(sapply(seq_along(D_calib), function(k) {
      masked_projection(mask_structural(H, M_calib_t[k, ])) %*% (Y[D_calib[k], ] - Yhat[D_calib[k], ])
    }))
    
    mu_t_proj <- t(sapply(seq_along(D_calib), function(k) {
      masked_projection(mask_structural(H, M_calib_t[k, ])) %*% (Yhat[t, ])
    }))
    
    Scores_up <- mu_t_proj + abs(Scores)
    Scores_down <- mu_t_proj - abs(Scores)
    
    # We dont take into account the component we masked in the score to compute the quantiles
    Scores_up[M_calib_t == 0] <- NA
    Scores_down[M_calib_t == 0] <- NA
    # We implicitely require that at least some calibration points have each of the test component
    
    # Bound MNR Nested
    upper_nest[i, ] <- apply(Scores_up, 2, function(s) quantile(s, probs = 1 - alpha/2, type = 1, na.rm = TRUE))
    lower_nest[i, ] <- -apply(-Scores_down, 2, function(s) quantile(s, probs = 1 - alpha/2, type = 1, na.rm = TRUE))
  }
  
  # Empirical coverage and length
  y_test <- Y[D_test, ]
  M_test <- M_mat[D_test, ]
  
  # Mask the metrics for the unobserved components
  covered <- (y_test >= lower_nest) & (y_test <= upper_nest)
  covered[M_test == 0] <- NA
  
  lengths <- upper_nest - lower_nest
  lengths[M_test == 0] <- NA
  
  # Marginal coverages and lengths
  coverage_vec <- colMeans(covered, na.rm = TRUE)
  coverage_all <- mean(covered, na.rm = TRUE)
  length_vec   <- colMeans(lengths, na.rm = TRUE)
  length_all   <- mean(length_vec)
  
  # Mask--Conditionnal coverage and length
  mask_conditional_coverage <- compute_masked_coverage(covered, M_test)
  mask_conditional_length   <- compute_masked_length(lengths, M_test)
  
  
  list(
    lower_nest = lower_nest,
    upper_nest = upper_nest,
    covered = covered,
    coverage_vec = coverage_vec,
    coverage_all = coverage_all,
    length_vec = length_vec,
    length_all = length_all,
    mask_conditional = mask_conditional_coverage,
    mask_conditional_length = mask_conditional_length
  )
}

CP_MNR_Nested_star <- function(D_calib, D_test, Y, Yhat, M_mat, H_array, alpha = 0.1, max_diff = 1){
  # Initialize the matrices to collect the prediction sets
  n_test <- length(D_test)
  m <- ncol(M_mat)
  lower_nest <- matrix(NA, nrow = n_test, ncol = m)
  upper_nest <- matrix(NA, nrow = n_test, ncol = m)
  # For each test point, we compute new scores
  for (i in seq_along(D_test)) {
    t <- D_test[i]
    print(i)
    # Filter the calibration point the keep only the similar structures
    diff_mask <- rowSums(M_mat[D_calib, ] != matrix(M_t, nrow = length(D_calib), ncol = m, byrow = TRUE))
    idx_calib_filtered <- D_calib[diff_mask <= max_diff]
    if (length(idx_calib_filtered) == 0) {
      warning(sprintf("Test point %d: no calibration points associated", t))
      next
    }
    M_calib_t <- 1*(M_mat[idx_calib_filtered, ] & matrix(M_t, nrow = length(idx_calib_filtered), ncol = m, byrow = TRUE))
    # Manual fix to deal with the zero node case
    all_children_missing <- rowSums(M_calib_t[, 2:(m), drop = FALSE]) == 0
    M_calib_t[all_children_missing, 1] <- 0
    # Scores MNR Nested
    n_calib <- length(idx_calib_filtered)
    
    Scores <- t(sapply(1:n_calib, function(ii) {
      masked_projection(mask_structural(H, M_calib_t[ii, ])) %*% 
        (Y[idx_calib_filtered[ii], ] - Yhat[idx_calib_filtered[ii], ])
    }))
    
    mu_t_proj <- t(sapply(1:n_calib, function(ii) {
      masked_projection(mask_structural(H, M_calib_t[ii, ])) %*% Yhat[t, ]
    }))
    
    Scores_up <- mu_t_proj + abs(Scores)
    Scores_down <- mu_t_proj - abs(Scores)
    
    # We dont take into account the component we masked in the score to compute the quantiles
    Scores_up[M_calib_t == 0] <- NA
    Scores_down[M_calib_t == 0] <- NA
    
    nested_component <- function(sd, su, alpha) {
      # Keep only valid calibration points for this component
      valid_idx <- which(!is.na(sd) & !is.na(su))
      if(length(valid_idx) == 0) return(c(lower = NA, upper = NA))
      sd <- sd[valid_idx]
      su <- su[valid_idx]
      n <- length(sd)
      # Construct events for sweep line
      events_x     <- c(sd, su)
      # Assign +1 for lower bounds and -1 for upper bounds
      events_delta <- c(rep(1, n), rep(-1, n))
      # Sort events: increasing for lower bounds, decreasing for upper bounds
      events_order <- order(events_x, -events_delta)
      events_x     <- events_x[events_order]
      events_delta <- events_delta[events_order]
      
      # Cumulative sum to count number of temporary sets each point belongs to
      counts <- cumsum(events_delta)
      # Select points that belong to enough calibration intervals
      inside <- counts >= floor(alpha*(n + 1) + 1)
      # Include the last upper bound to ensure intervals are closed
      inside <- inside | c(inside[-1], FALSE)
      
      # Identify contiguous intervals
      diff_inside <- c(TRUE, diff(inside) != 0)
      # Here, we consider the convex hull to get intervals
      interval_starts <- min(events_x[inside])
      interval_ends   <- max(events_x[inside])
      # Return the first interval (or adapt to return all intervals)
      if(length(interval_starts) == 0) return(c(lower = NA, upper = NA))
      c(lower = interval_starts, upper = interval_ends)
    }
    result <- lapply(1:m, function(j) {
      nested_component(Scores_down[, j], Scores_up[, j], alpha)
    })
    # Bound MNR Nested
    for(j in 1:m) {
      lower_nest[i,j] <- result[[j]]["lower"]
      upper_nest[i,j] <- result[[j]]["upper"]
    }
  }
  
  # Empirical coverage and length
  y_test <- Y[D_test, ]
  M_test <- M_mat[D_test, ]
  
  # Mask the metrics for the unobserved components
  covered <- (y_test >= lower_nest) & (y_test <= upper_nest)
  covered[M_test == 0] <- NA
  
  lengths <- upper_nest - lower_nest
  lengths[M_test == 0] <- NA
  
  # Marginal coverages and lengths
  coverage_vec <- colMeans(covered, na.rm = TRUE)
  coverage_all <- mean(covered, na.rm = TRUE)
  length_vec   <- colMeans(lengths, na.rm = TRUE)
  length_all   <- mean(length_vec)
  
  # Mask--Conditionnal coverage and length
  mask_conditional_coverage <- compute_masked_coverage(covered, M_test)
  mask_conditional_length   <- compute_masked_length(lengths, M_test)
  
  
  list(
    lower_nest = lower_nest,
    upper_nest = upper_nest,
    covered = covered,
    coverage_vec = coverage_vec,
    coverage_all = coverage_all,
    length_vec = length_vec,
    length_all = length_all,
    mask_conditional = mask_conditional_coverage,
    mask_conditional_length = mask_conditional_length
  )
}