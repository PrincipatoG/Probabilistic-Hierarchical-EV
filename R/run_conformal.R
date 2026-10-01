#!/usr/bin/env Rscript

rm(list = ls())

source("R/forecast_function.R")
source("R/reconciliation_function.R")
source("R/metric_function.R")
source("R/conformal_methods.R")

library(argparser)
library(dplyr)
library(readr)
library(stringr)
library(purrr)
library(ggplot2)
library(mgcv)
library(MASS)
library(fastmatrix)
library(Matrix)
library(RSpectra)

# -----------------------------
# ARGUMENT PARSER
# -----------------------------

p <- arg_parser("Conformal prediction experiment")

p <- add_argument(p, 
                  "--seed", 
                  help="Seed for generating windows", 
                  default=40, 
                  type="integer")

p <- add_argument(
  p,
  "--conformal_method",
  help = "Conformal prediction method",
  default = "CP_MNR",
  type="character"
)

p <- add_argument(
  p,
  "--projection_methods",
  help = "Comma-separated projection methods",
  default = "Direct,OLS,refined OLS",
  type="character"
)

p <- add_argument(
  p,
  "--national_model",
  help = "National forecasting model",
  default = "Combination",
  type="character"
)

p <- add_argument(
  p,
  "--regional_model",
  help = "Regional forecasting model",
  default = "Combination",
  type="character"
)

p <- add_argument(
  p,
  "--results_national",
  help = "Path to national forecasts",
  default = "results/results_scotland.RDS",
  type="character"
)

p <- add_argument(
  p,
  "--results_regional",
  help = "Path to regional forecasts",
  default = "results/results_regions.RDS",
  type="character"
)

p <- add_argument(
  p,
  "--results_station",
  help = "Path to station forecasts",
  default = "results/results_stations.RDS",
  type="character"
)
p <- add_argument(
  p,
  "--station_model",
  help = "Station forecasting model",
  default = "Combination",
  type="character"
)

p <- add_argument(
  p,
  "--alpha",
  help = "Miscoverage level",
  type = "numeric",
  default = 0.1
)

p <- add_argument(
  p,
  "--ncores",
  help = "Number of cores",
  type = "integer",
  default = 14
)

p <- add_argument(
  p,
  "--output",
  help = "Output path",
  type="character",
  default = NULL
)

p <- add_argument(
  p,
  "--structural_matrix",
  help = "Path to structural hierarchy matrix",
  type="character",
  default = "Data/structural.RDS"
)

argv <- parse_args(p)

# -----------------------------
# PARAMETERS
# -----------------------------

seed_windows  <- argv$seed

model_names <- list(
  argv$national_model,
  argv$regional_model,
  argv$station_model
)

projection_methods <- strsplit(
  argv$projection_methods,
  ","
)[[1]]

conformal_method <- argv$conformal_method

results_national_path <- paste0("Output_ponctual/Seed_", seed_windows, "/results_scotland.RDS")
results_regional_path <- paste0("Output_ponctual/Seed_", seed_windows, "/results_regions.RDS")
results_station_path  <- paste0("Output_ponctual/Seed_", seed_windows, "/results_stations.RDS")

structural_matrix_path <- argv$structural_matrix

# -----------------------------
# BUILD OUTPUT PATH
# -----------------------------

model_tag <- paste(model_names, collapse = "__")

projection_tag <- paste(
  gsub(" ", "_", projection_methods),
  collapse = "__"
)

default_output <- paste0(
  "Output_conformal/Seed_",
  seed_windows,
  "/",
  conformal_method,
  "__",
  projection_tag,
  "__",
  model_tag,
  ".RDS"
)

output_path <- ifelse(
  is.na(argv$output),
  default_output,
  argv$output
)

# -----------------------------
# LOAD DATA
# -----------------------------

res_nat <- readRDS(results_national_path)
res_reg <- readRDS(results_regional_path)
res_sta <- readRDS(results_station_path)

H <- readRDS(structural_matrix_path)

region_levels  <- sub("Region_", "", grep("^Region_", names(H), value = TRUE))
station_levels <- sub("Station_", "", grep("^Station_", names(H), value = TRUE))

P_OLS <- refined_OLS_projection(
  as.matrix(H[, -1])
)

windows <- generate_rolling_windows(seed_windows)

# -----------------------------
# RUN EXPERIMENT
# -----------------------------

df_all <- run_conformal_experiment(
  windows = windows,
  res_nat = res_nat,
  res_reg = res_reg,
  res_sta = res_sta,
  H = H,
  P_OLS = P_OLS,
  model_names = model_names,
  conformal_method = conformal_method,
  projection_methods = projection_methods,
  alpha = argv$alpha,
  ncores = argv$ncores
)

# -----------------------------
# EXPORT
# -----------------------------

saveRDS(df_all, output_path)
