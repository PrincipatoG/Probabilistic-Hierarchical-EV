# --- Previous steps are not automatized : they include scrapping and formatting of data ---
rm(list = ls())

log_dir <- "logs"
dir.create(log_dir, showWarnings = FALSE)
# --- Set seed ---
seed_windows <- 1
dir.create(paste0('Data/Seed_', seed_windows), showWarnings = TRUE, recursive = FALSE, mode = "0777")
dir.create(paste0('Output_conformal/Seed_', seed_windows), showWarnings = TRUE, recursive = FALSE, mode = "0777")
dir.create(paste0('Output_ponctual/Seed_', seed_windows), showWarnings = TRUE, recursive = FALSE, mode = "0777")
dir.create(paste0('logs/Seed_', seed_windows), showWarnings = TRUE, recursive = FALSE, mode = "0777")

# # --- Run R prediction files ---
# message("Running: R prediction files")
# 
# message("Running: run_scotland.R")
# 
# system(paste("caffeinate Rscript scripts/run_scotland.R",
#              "--seed", seed_windows,
#              "--output", paste0('Output_ponctual/Seed_', seed_windows,"/results_scotland.RDS") ))
# 
# message("Running: run_regions.R")
# 
# system(paste("caffeinate Rscript scripts/run_regions.R",
#              "--seed", seed_windows,
#              "--output", paste0('Output_ponctual/Seed_', seed_windows,"/results_regions") ))
# 
# message("Running: run_stations.R")
# 
# system(paste("caffeinate Rscript scripts/run_stations.R",
#              "--seed", seed_windows,
#              "--output", paste0('Output_ponctual/Seed_', seed_windows,"/results_stations") ))
# 
# # --- Preparing datasets for Python ---
# message("Preparing: datasets for Python")
# 
# message("Running: prep_data_tabular_national.R")
# 
# system(paste("caffeinate Rscript scripts_tab/prep_data_tabular_national.R",
#              "--seed", seed_windows))
# 
# message("Running: prep_data_tabular_regions.R")
# 
# system(paste("caffeinate Rscript scripts_tab/prep_data_tabular_regions.R",
#              "--seed", seed_windows))
# 
# message("Running: run_stations.R")
# 
# system(paste("caffeinate Rscript scripts_tab/prep_data_tabular_stations.R",
#              "--seed", seed_windows,
#              "--output", 
#              paste0('Output_ponctual/Seed_', seed_windows,"/results_stations") ))
# 
# # --- Run python prediction files ---
# message("Running: python prediction files")
# 
# message("Running: tabICL_national.py")
# 
# system(paste("caffeinate python3.11 scripts_tab/tabICL_national.py",
#              "--seed", seed_windows))
# 
# message("Running: tabICL_regions_GPU.py")
# 
# system(paste("caffeinate python3.11 scripts_tab/tabICL_regions_GPU.py",
#              "--seed", seed_windows))
# 
# message("Running: tabICL_stations_GPU.py")
# 
# system(paste("caffeinate python3.11 scripts_tab/tabICL_stations_GPU.py",
#              "--seed", seed_windows))
# 
# # --- Gather point forecasts and perform forecast reconciliation ---
# 
# message("Performing forecast reconciliation for point forecasting")
# 
# message("Running: forecast_reconciliation.R")
# 
# system(paste("caffeinate Rscript R/forecast_reconciliation.R",
#        "--seed", seed_windows,
#        "--run_gather", TRUE
#        ))


# --- Run conformal scripts ---
methods <- c(
  "CP_MNR"#,
  # "CP_MNR_naive",
  # "Adaptive_CP_MNR",
  # "CP_MNR_Nested_star"
)
log_dir <- paste0('logs/Seed_', seed_windows)
for (m in methods) {

  log_file <- file.path(log_dir, paste0("log_", m, ".txt"))

  cmd <- paste(
    "caffeinate Rscript R/run_conformal.R",
    "--conformal_method", m,
    "--national_model", "Combination",
    "--regional_model", "Combination",
    "--station_model", "Combination",
    "--seed", seed_windows
  )

  message("Running: ", m)
  system(cmd)
}