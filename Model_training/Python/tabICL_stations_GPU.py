# Generic Packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import random
import os
import sys
from tqdm import tqdm
import joblib
import itertools
import statsmodels.api as sm
from statsmodels.graphics.tsaplots import plot_acf, plot_pacf

# Import torch and configure MPS (Metal Performance Shaders)
import torch
import torch.nn as nn
import torch.optim
import argparse
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

# Detect Apple Silicon chip (MPS)
if torch.backends.mps.is_available():
    device = torch.device("mps")
    print("🚀 MPS GPU detected. PyTorch will use the Apple Silicon GPU.")
else:
    device = torch.device("cpu")
    print("⚠️ MPS GPU not detected. Using CPU.")

# Time series
from statsmodels.tsa.arima.model import ARIMA
# Metric
from sklearn.metrics import root_mean_squared_error as rmse

def nrmse(actual, predicted):
    # Convert to numpy arrays for vectorized operations
    actual = np.asarray(actual)
    predicted = np.asarray(predicted)
    
    # rmse <- sqrt(mean((actual - predicted)^2, na.rm = TRUE))
    rmse_val = np.sqrt(np.nanmean((actual - predicted) ** 2))
    
    # n_rmse <- rmse / mean(actual, na.rm = TRUE)
    n_rmse_val = rmse_val / np.nanmean(actual)
    
    return n_rmse_val

def nmae(actual, predicted):
    # Convert to numpy arrays for vectorized calculations
    actual = np.asarray(actual)
    predicted = np.asarray(predicted)
    
    # mae <- mean(abs(actual - predicted), na.rm = TRUE)
    mae_val = np.nanmean(np.abs(actual - predicted))
    
    # n_mae <- mae / mean(actual, na.rm = TRUE)
    n_mae_val = mae_val / np.nanmean(actual)
    
    return n_mae_val

# FM
from tabicl import TabICLClassifier, TabICLRegressor
# from chronos import BaseChronosPipeline, Chronos2Pipeline

os.chdir(f'Data/Seed_{seed}')

# 1. Load data
print("Loading historical data...")
Data0 = pd.read_csv('dataset_stations_for_tab.csv', on_bad_lines='warn', sep=';')
Data0['Date'] = pd.to_datetime(Data0['Date'])
print(f"Number of rows in Data0: {len(Data0)}")

# Define columns
cat_cols = ['Weekday_Holiday_BIS', 'Region', 'Station.ID']
num_cols = ['Posan', 'tmpf_max', 'Latitude', 'relh_mean', 'max_cp_30days', 'taux_cp_paid_max30',
            'taux_lent', 'taux_rapide', 'taux_accelere', 'Lag_Mean_1_7', 'n_cp_lag1', 'Lag7', 'Lag1']

target_col = 'Consumed_kWh'
useful_col = cat_cols + num_cols + ["type", target_col]

# Ensembling parameters
num_models = 5      # Number of experts
sample_size = 10000 # Sampling size for bagging
batch_size = 2**10 # Inference chunk size (1024); adjust for memory constraints

# List to store results for all windows
all_results = []

# 2. Loop over windows 1 to 26
for window_id in range(1, 15):
    window_path = f"my_windows/window_{window_id}.csv"

    # Check whether the file exists
    if not os.path.exists(window_path):
        print(f"File {window_path} not found. Skipping.")
        continue

    print(f"\n=======================================================")
    print(f"--- Processing Window {window_id} ---")

    # Load and merge
    Window = pd.read_csv(window_path)
    Window["Date"] = pd.to_datetime(Window["Date"])

    df_fusionne = pd.merge(Data0, Window, on="Date", how="right")

    # Clean missing values
    # Keep 'Date' for the final export
    df_propre = df_fusionne[["Date"] + useful_col].dropna().copy()

    # Split Train / Test
    is_train = df_propre["type"] == "train"
    df_train_full = df_propre.loc[is_train].copy()
    
    is_train = df_fusionne["type"] == "train"
    X_test = df_fusionne.loc[~(df_fusionne['type'] =='train'), cat_cols + num_cols].copy()
    y_test = df_fusionne.loc[~(df_fusionne['type'] =='train'), target_col].copy()


    # Skip the window if it has no test data
    if X_test.empty:
        print(f"No test rows in Window {window_id}. Skipping.")
        continue

    X_test.columns = [str(c) for c in X_test.columns]

    # --- START ENSEMBLING FOR THIS WINDOW ---
    all_predictions = []
    
    for i in range(num_models):
        print(f"  -> Training model {i+1}/{num_models} (Window {window_id})...")
        
        # 1. Random sampling
        train_sample = df_train_full.sample(n=sample_size, random_state=42 + i)
        X_train_sub = train_sample[cat_cols + num_cols].copy()
        y_train_sub = train_sample[target_col].copy()
        
        X_train_sub.columns = [str(c) for c in X_train_sub.columns]
        
        # 2. Initialize model
        reg = TabICLRegressor(verbose=False, random_state=30 + i, n_estimators=1, device=device)
        
        # 3. Training
        reg.fit(X_train_sub, y_train_sub.values.ravel())
        
        # 4. Inference (batched to avoid RAM issues)
        pred_tabICL = []
        for j in tqdm(range(0, len(X_test), batch_size), desc="     Prediction", leave=False):
            X_batch = X_test.iloc[j:j+batch_size]
            batch_preds = reg.predict(X_batch)
            pred_tabICL.extend(batch_preds)
            
        # Add predictions to the ensemble
        all_predictions.append(pred_tabICL)
        
        # Clear GPU cache between models
        if torch.backends.mps.is_available():
            torch.mps.empty_cache()

    # 5. Ensembling: average predictions for this window
    ensemble_preds = np.mean(all_predictions, axis=0)
    # --- END ENSEMBLING ---

    # Prepare the final DataFrame for this window
    is_train = df_fusionne["type"] == "train"  
    df_test = df_fusionne.loc[~is_train].copy()
    df_test["pred_tabICL"] = ensemble_preds
    df_test["window_id"] = window_id

    # Arrange columns for output
    cols_finales = ["window_id", "Date"] + useful_col + ["pred_tabICL"]
    all_results.append(df_test[cols_finales])

    nom_fichier_export_provisoire = f"../../Output_ponctual/Seed_{seed}/results_stations_type_tabICL_1398/window_{window_id}.csv"
    df_test[cols_finales].to_csv(nom_fichier_export_provisoire, index=False, sep=';')



# 3. Merge all windows and export to CSV
if all_results:
    df_final = pd.concat(all_results, ignore_index=True)

    # Export to CSV
    nom_fichier_export = f"../../Output_ponctual/Seed_{seed}/results_stations_tabICL.csv"
    
    # Create the 'results' directory if it does not exist
    os.makedirs(os.path.dirname(nom_fichier_export), exist_ok=True)
    
    # Export with a semicolon separator
    df_final.to_csv(nom_fichier_export, index=False, sep=';')

    print(f"\n### Processing completed! ###")
    print(f"File exported to: {nom_fichier_export}")
    print(f"Total number of test rows recorded: {len(df_final)}")
else:
    print("\nNo test data could be processed.")
