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

build_vectors <- function(period,
                          res_nat,
                          res_reg,
                          res_sta,
                          model_names) {
  
  # --- Filter once ---
  res_nat <- res_nat[res_nat$Date %in% period, ]
  res_reg <- res_reg[res_reg$Date %in% period, ]
  res_sta <- res_sta[res_sta$Date %in% period, ]
  
  # Order (important for temporal alignment)
  res_nat <- res_nat[order(res_nat$Date), ]
  res_reg <- res_reg[order(res_reg$Date), ]
  res_sta <- res_sta[order(res_sta$Date), ]
  
  # --- NATIONAL (simple vectors) ---
  Y_nat    <- res_nat$Consumed_kWh
  Yhat_nat <- res_nat[[model_names[[1]]]]
  
  # --- REGIONS (pivot to matrix) ---
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
  
  # drop Date + fast conversion
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
  
  # --- CONCAT (once) ---
  Y    <- cbind(Y_nat, Y_reg, Y_sta)
  Yhat <- cbind(Yhat_nat, Yhat_reg, Yhat_sta)
  
  # --- Efficient construction ---
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
