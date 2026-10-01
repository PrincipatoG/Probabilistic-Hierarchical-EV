nmae <- function(actual, predicted) {
  mae <- mean(abs(actual - predicted), na.rm = TRUE)
  n_mae <- mae / mean(actual, na.rm = TRUE)
  return(n_mae)
}

nrmse <- function(actual, predicted) {
  rmse <- sqrt(mean((actual - predicted)^2, na.rm = TRUE))
  n_rmse <- rmse / mean(actual, na.rm = TRUE)
  return(n_rmse)
}