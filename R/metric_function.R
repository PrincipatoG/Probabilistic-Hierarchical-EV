# -----------------------------
# POINT FORECASTING
# -----------------------------

mse_by_node <- function(df, method = "OLS_ref") {
  
  pred_col <- method
  
  df %>%
    group_by(node) %>%
    summarise(
      n = sum(!is.na(y) & !is.na(.data[[pred_col]])),
      mse = ifelse(
        n == 0,
        NA_real_,
        mean((y - .data[[pred_col]])^2, na.rm = TRUE)
      ),
      .groups = "drop"
    )
}

mse_by_type <- function(df, method = "OLS_ref") {
  
  node_mse <- mse_by_node(df, method)
  
  node_info <- df %>%
    distinct(node, node_type)
  
  node_mse <- node_mse %>%
    left_join(node_info, by = "node")
  
  node_mse %>%
    group_by(node_type) %>%
    summarise(
      total_n = sum(n, na.rm = TRUE),
      mse = sum(mse * n, na.rm = TRUE) / total_n,
      .groups = "drop"
    )
}

mse_global <- function(df, method = "OLS_ref") {
  
  pred_col <- method
  
  valid <- !is.na(df$y) & !is.na(df[[pred_col]])
  
  mean((df$y[valid] - df[[pred_col]][valid])^2)
}

# -----------------------------
# NMAE
# -----------------------------

nmae_by_node <- function(df, method = "OLS_ref") {
  
  pred_col <- method
  
  df %>%
    group_by(node) %>%
    summarise(
      n = sum(!is.na(y) & !is.na(.data[[pred_col]])),
      nmae = ifelse(
        n == 0,
        NA_real_,
        sum(abs(y - .data[[pred_col]]), na.rm = TRUE) /
          sum(abs(y), na.rm = TRUE)
      ),
      .groups = "drop"
    )
}

nmae_by_type <- function(df, method = "OLS_ref") {
  
  node_nmae <- nmae_by_node(df, method)
  
  node_info <- df %>%
    distinct(node, node_type)
  
  node_nmae <- node_nmae %>%
    left_join(node_info, by = "node")
  
  df %>%
    group_by(node_type) %>%
    summarise(
      nmae = sum(abs(y - .data[[method]]), na.rm = TRUE) /
        sum(abs(y), na.rm = TRUE),
      .groups = "drop"
    )
}

nmae_global <- function(df, method = "OLS_ref") {
  
  pred_col <- method
  
  valid <- !is.na(df$y) & !is.na(df[[pred_col]])
  
  sum(abs(df$y[valid] - df[[pred_col]][valid])) /
    sum(abs(df$y[valid]))
}

# -----------------------------
# PROBABILISTIC FORECASTING
# -----------------------------
node_metrics <- function(df, min_obs = 14) {
  
  df %>%
    group_by(node, node_type, method) %>%
    summarise(
      coverage = weighted.mean(coverage, w = n_eff, na.rm = TRUE),
      length   = weighted.mean(length,   w = n_eff, na.rm = TRUE),
      n_eff    = sum(n_eff),
      .groups = "drop"
    ) %>%
    filter(n_eff >= min_obs)
}

global_metrics <- function(df, mask_df) {
  
  df %>%
    dplyr::group_by(method) %>%
    dplyr::summarise(
      coverage = weighted.mean(coverage, w = n_eff, na.rm = TRUE),
      length   = weighted.mean(length,   w = n_eff, na.rm = TRUE),
      .groups = "drop"
    )
}
type_metrics <- function(df) {
  
  df %>%
    dplyr::group_by(node_type, method) %>%
    dplyr::summarise(
      coverage = weighted.mean(coverage, w = n_eff, na.rm = TRUE),
      length   = weighted.mean(length,   w = n_eff, na.rm = TRUE),
      .groups = "drop"
    )
}
plot_by_type <- function(df_all, alpha = 0.1, min_obs = 14) {
  
  df_nodes <- node_metrics(df_all, min_obs)
  df_type  <- type_metrics(df_nodes)
  
  types <- unique(df_nodes$node_type)
  
  plots <- lapply(types, function(tp) {
    
    df_n <- dplyr::filter(df_nodes, node_type == tp)
    df_t <- dplyr::filter(df_type,  node_type == tp)
    
    ggplot() +
      
      # --- Node points (transparent) ---
      geom_point(
        data = df_n,
        aes(
          x = coverage,
          y = length,
          color = method,
          shape = method,
          alpha = n_eff
        ),
        size = 2
      ) +
      
      # --- Mean points ---
      geom_point(
        data = df_t,
        aes(
          x = coverage,
          y = length,
          color = method,
          shape = method,
          
        ),
        size = 5
      ) +
      
      scale_x_continuous(limits = c(0.6, 1)) +
      
      geom_vline(xintercept = 1 - alpha,
                 linetype = "dashed",
                 color = "black") +
      
      scale_alpha(range = c(0.05, 0.2), guide = "none") +
      
      # --- Filled shapes only ---
      scale_shape_manual(values = c(16, 17, 15, 18, 19, 16, 17, 15, 18, 19)) +
      
      labs(
        title = tp,
        x = "Coverage",
        y = "Length"
        ) +
      
      theme_minimal() +
      theme(
        legend.position = "bottom"
      )
  })
  
  names(plots) <- types
  return(plots)
}
