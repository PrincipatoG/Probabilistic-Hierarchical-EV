"""
@author: Yannig
"""

# =============================================================================
# Generic packages
# =============================================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import random
import os
import argparse
import torch
import torch.nn as nn
import torch.optim

# =============================================================================
# Parse command-line arguments
# =============================================================================

parser = argparse.ArgumentParser()
parser.add_argument("--seed", type=int, default=0, help="Random seed")
args = parser.parse_args()

seed = args.seed
print(f"Using seed: {seed}")

# Set random seeds for reproducibility
random.seed(seed)
np.random.seed(seed)
torch.manual_seed(seed)

# =============================================================================
# GPU configuration (Apple Silicon / MPS)
# =============================================================================

if torch.backends.mps.is_available():
    device = torch.device("mps")
    print("MPS GPU detected. Using Apple Silicon GPU.")
else:
    device = torch.device("cpu")
    print("⚠️ MPS GPU not available. Using CPU.")

# =============================================================================
# TabICL imports
# =============================================================================

from tabicl import TabICLRegressor

# =============================================================================
# Metrics
# =============================================================================

def nrmse(actual, predicted):
    actual = np.asarray(actual)
    predicted = np.asarray(predicted)

    rmse_val = np.sqrt(np.nanmean((actual - predicted) ** 2))
    return rmse_val / np.nanmean(actual)


def nmae(actual, predicted):
    actual = np.asarray(actual)
    predicted = np.asarray(predicted)

    mae_val = np.nanmean(np.abs(actual - predicted))
    return mae_val / np.nanmean(actual)

# =============================================================================
# Set working directory
# =============================================================================

os.chdir(f"Data/Seed_{seed}")

# =============================================================================
# Load historical dataset
# =============================================================================

Data0 = pd.read_csv("dataset_regions_for_tab.csv")
Data0["Date"] = pd.to_datetime(Data0["Date"])

print(f"Number of rows in historical dataset: {len(Data0)}")

# =============================================================================
# Define columns
# =============================================================================

cat_cols = ["Weekday_Holiday_BIS", "Region"]

num_cols = [
    "Posan",
    "tmpf_max",
    "Latitude",
    "relh_mean",
    "max_cp_30days",
    "taux_cp_paid_max30",
    "taux_lent",
    "taux_rapide",
    "taux_accelere",
    "Lag_Mean_1_7",
    "n_cp_lag1",
    "Lag7",
    "Lag1"
]

target_col = "Consumed_kWh"

useful_col = cat_cols + num_cols + ["type", target_col]

# =============================================================================
# Ensemble parameters
# =============================================================================

num_models = 5
sample_size = 10000

# List to store all window results
all_results = []

# =============================================================================
# Loop over all windows
# =============================================================================

for window_id in range(1, 15):

    window_path = f"my_windows/window_{window_id}.csv"

    # Skip missing files
    if not os.path.exists(window_path):
        print(f"File {window_path} not found. Skipping.")
        continue

    print("\n=======================================================")
    print(f"Processing Window {window_id}")

    # -------------------------------------------------------------------------
    # Load current window
    # -------------------------------------------------------------------------

    Window = pd.read_csv(window_path)
    Window["Date"] = pd.to_datetime(Window["Date"])

    # Merge historical data with current window
    df_fusionne = pd.merge(Data0, Window, on="Date", how="right")

    # Keep all rows (including NA) for final export
    df_all = df_fusionne[["Date"] + useful_col].copy()

    # Keep only complete rows for training/prediction
    required_cols = cat_cols + num_cols + [target_col, "type"]
    df_propre = df_all.dropna(subset=required_cols).copy()

    # Split train / test
    is_train = df_propre["type"] == "train"
    df_train_full = df_propre.loc[is_train].copy()

    X_test = df_propre.loc[~is_train, cat_cols + num_cols].copy()
    y_test = df_propre.loc[~is_train, target_col].copy()

    # Skip if no test data
    if X_test.empty:
        print(f"No test rows in Window {window_id}. Skipping.")
        continue

    X_test.columns = [str(c) for c in X_test.columns]

    # -------------------------------------------------------------------------
    # Ensemble training
    # -------------------------------------------------------------------------

    all_predictions = []

    for model_idx in range(num_models):

        print(
            f"  -> Training model {model_idx + 1}/{num_models} "
            f"(Window {window_id})"
        )

        # Random subsampling of the training set
        train_sample = df_train_full.sample(
            n=min(sample_size, len(df_train_full)),
            random_state=seed + model_idx
        )

        X_train_sub = train_sample[cat_cols + num_cols].copy()
        y_train_sub = train_sample[target_col].copy()

        X_train_sub.columns = [str(c) for c in X_train_sub.columns]

        # Initialize TabICL model
        reg = TabICLRegressor(
            verbose=False,
            random_state=seed + model_idx,
            n_estimators=1,
            device=device
        )

        # Fit model
        reg.fit(X_train_sub, y_train_sub.values.ravel())

        # Predict on complete test rows
        preds = reg.predict(X_test)

        # Store predictions
        all_predictions.append(preds)

        # Clear GPU cache (important for MPS memory management)
        if torch.backends.mps.is_available():
            torch.mps.empty_cache()

    # -------------------------------------------------------------------------
    # Aggregate ensemble predictions
    # -------------------------------------------------------------------------

    ensemble_preds = np.mean(all_predictions, axis=0)

    # -------------------------------------------------------------------------
    # Build final test dataframe while preserving NA rows
    # -------------------------------------------------------------------------

    # Keep all test rows (including incomplete ones)
    df_test = df_all.loc[df_all["type"] != "train"].copy()

    # Initialize prediction column with NA
    df_test["pred_tabICL"] = np.nan

    # Identify complete rows only
    mask_complete = df_test[cat_cols + num_cols + [target_col]].notna().all(axis=1)
    
    # Assign predictions only to complete rows
    df_test.loc[mask_complete, "pred_tabICL"] = ensemble_preds

    # Add window identifier
    df_test["window_id"] = window_id

    # Keep final column order
    cols_finales = ["window_id", "Date"] + useful_col + ["pred_tabICL"]

    all_results.append(df_test[cols_finales])

# =============================================================================
# Concatenate and export results
# =============================================================================

if all_results:

    df_final = pd.concat(all_results, ignore_index=True)

    output_path = f"../../Output_ponctual/Seed_{seed}/results_regions_tabICL.csv"

    # Create results directory if it does not exist
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    # Export CSV
    df_final.to_csv(output_path, index=False)

    print("\n### Processing completed ###")
    print(f"Results exported to: {output_path}")
    print(f"Total number of test rows: {len(df_final)}")

    # -------------------------------------------------------------------------
    # Compute final metrics (NA-safe)
    # -------------------------------------------------------------------------

    score_nrmse = nrmse(
        df_final[target_col],
        df_final["pred_tabICL"]
    )

    score_nmae = nmae(
        df_final[target_col],
        df_final["pred_tabICL"]
    )

    print("\n### Final scores ###")
    print(f"NRMSE: {score_nrmse:.4f}")
    print(f"NMAE:  {score_nmae:.4f}")

else:
    print("\nNo test data could be processed.")
