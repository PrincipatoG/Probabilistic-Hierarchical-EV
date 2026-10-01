"""
@author: Yannig
"""

#Generic Packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.pyplot as plt
import random
import os
from tqdm import tqdm
import joblib
import itertools
import statsmodels as statsmodels
from statsmodels.graphics.tsaplots import plot_acf, plot_pacf
import argparse
import sys

#Import torch
import torch.nn as nn
import torch.optim
import torch
#Time series
from statsmodels.tsa.arima.model import ARIMA
#Metric
import numpy as np
from sklearn.metrics import mean_squared_error

rmse = lambda y_true, y_pred: np.sqrt(mean_squared_error(y_true, y_pred))
def nrmse(actual, predicted):
    # Convert to numpy arrays for vectorized operations
    actual = np.asarray(actual)
    predicted = np.asarray(predicted)
    
    # rmse <- sqrt(mean((actual - predicted)^2, na.rm = TRUE))
    rmse = np.sqrt(np.nanmean((actual - predicted) ** 2))
    
    # n_rmse <- rmse / mean(actual, na.rm = TRUE)
    n_rmse = rmse / np.nanmean(actual)
    
    return n_rmse


def nmae(actual, predicted):
    # Convert to numpy arrays for vectorized calculations
    actual = np.asarray(actual)
    predicted = np.asarray(predicted)
    
    # mae <- mean(abs(actual - predicted), na.rm = TRUE)
    mae = np.nanmean(np.abs(actual - predicted))
    
    # n_mae <- mae / mean(actual, na.rm = TRUE)
    n_mae = mae / np.nanmean(actual)
    
    return n_mae



#FM
from tabicl import TabICLClassifier, TabICLRegressor
#from chronos import BaseChronosPipeline, Chronos2Pipeline

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

#import os
os.chdir(f'Data/Seed_{seed}')

###########################################################################################
#######Test on one window
###########################################################################################


# Historical data used to train the models
Data0 = pd.read_csv('dataset_scotland_for_tab.csv')
Data0['Date'] = pd.to_datetime(Data0['Date'])
n0 = len(Data0)
print(n0)

Window = pd.read_csv(f"my_windows/window_1.csv")
Window['Date'] = pd.to_datetime(Window['Date'])

df_fusionne = pd.merge(Data0, Window, on='Date', how='right')

###covariate types
cat_cols = ['Weekday_Holiday_BIS']
num_cols = ['Posan', 'Lag1', 'Lag7', 'Lag_Mean_1_7', 'max_paid_cp_30days', 'n_cp_paid_lag1',
            'tmpf_max']
target_col = 'Consumed_kWh'

useful_col= cat_cols + num_cols + ['type', target_col]
df_propre = df_fusionne[useful_col].dropna()


X_train = df_propre.loc[df_propre['type'] == 'train', cat_cols + num_cols].copy()
y_train = df_propre.loc[df_propre['type'] == 'train', target_col].copy()
X_test = df_propre.loc[~(df_propre['type'] =='train'), cat_cols + num_cols].copy()
y_test = df_propre.loc[~(df_propre['type'] =='train'), target_col].copy()

X_train.columns = [str(c) for c in X_train.columns]
X_test.columns = [str(c) for c in X_test.columns]

reg = TabICLRegressor( verbose=True, random_state=30, n_estimators=10)
reg.fit(X_train, y_train.values.ravel())
pred_tabICL = reg.predict(X_test)


score_nrmse = nrmse(y_test, pred_tabICL)
print(f"#########################################OFFline TabICL RMSE: {score_nrmse:.4f}")

score_nmae = nmae(y_test, pred_tabICL)
print(f"#########################################OFFline TabICL RMSE: {score_nmae:.4f}")



###########################################################################################
#######Loop over windows
###########################################################################################

# 1. Load historical data (fixed, outside the loop)
Data0 = pd.read_csv("dataset_scotland_for_tab.csv")
Data0["Date"] = pd.to_datetime(Data0["Date"])

# Define columns
cat_cols = ["Weekday_Holiday_BIS"]
num_cols = [
    "Posan",
    "Lag1",
    "Lag7",
    "Lag_Mean_1_7",
    "max_paid_cp_30days",
    "n_cp_paid_lag1",
    "tmpf_max",
]
target_col = "Consumed_kWh"
useful_col = cat_cols + num_cols + ["type", target_col]

# List to store results for each window
all_results = []

# 2. Loop over windows 1 to 14
for i in range(1, 15):
    window_path = f"my_windows/window_{i}.csv"

    # Check whether the file exists
    if not os.path.exists(window_path):
        print(f"File {window_path} not found. Skipping.")
        continue

    print(f"--- Processing Window {i} ---")

    # Load and merge
    Window = pd.read_csv(window_path)
    Window["Date"] = pd.to_datetime(Window["Date"])

    # Include 'Date' in the merge so it can be kept in the final CSV if needed
    df_fusionne = pd.merge(Data0, Window, on="Date", how="right")

    # Keep ALL rows
    df_all = df_fusionne[["Date"] + useful_col].copy()
    
    # Create a clean version for fit/predict only
    df_propre = df_all.dropna().copy()
    
    # Train on complete rows
    is_train = df_propre["type"] == "train"
    X_train = df_propre.loc[is_train, cat_cols + num_cols].copy()
    y_train = df_propre.loc[is_train, target_col].copy()
    
    # Test only on complete rows
    X_test = df_propre.loc[~is_train, cat_cols + num_cols].copy()

    # Skip the window if it has no test data
    if X_test.empty:
        print(f"No test rows in Window {i}. Skipping.")
        continue

    X_train.columns = [str(c) for c in X_train.columns]
    X_test.columns = [str(c) for c in X_test.columns]

    # Training and prediction
    reg = TabICLRegressor(verbose=False, random_state=30, n_estimators=10)
    reg.fit(X_train, y_train.values.ravel())

    # Filter ONLY the test set for this window's final result
    # Keep ALL test rows, including rows with NA
    df_test = df_all.loc[df_all["type"] != "train"].copy()
    
    # Initialize predictions with NA
    df_test["pred_tabICL"] = np.nan
    
    # Complete rows in the test set
    mask_complete = df_test[cat_cols + num_cols].notna().all(axis=1)
    
    # Predictions only where possible
    X_pred = df_test.loc[mask_complete, cat_cols + num_cols].copy()
    X_pred.columns = [str(c) for c in X_pred.columns]
    
    df_test.loc[mask_complete, "pred_tabICL"] = reg.predict(X_pred)
    df_test["window_id"] = i

    # Arrange columns for output (window ID first)
    cols_finales = ["window_id", "Date"] + useful_col + ["pred_tabICL"]
    all_results.append(df_test[cols_finales])

# 3. Merge all windows and export to CSV
if all_results:
    df_final = pd.concat(all_results, ignore_index=True)

    # Export
    nom_fichier_export = f"../../Output_ponctual/Seed_{seed}/results_scotland_tabICL.csv"
    df_final.to_csv(nom_fichier_export, index=False)

    print(
        f"\n### Processing completed! File exported to: {nom_fichier_export} ###"
    )
    print(f"Total number of test rows recorded: {len(df_final)}")
else:
    print("\nNo test data could be processed.")





df_final.columns


score_nrmse = nrmse(df_final['Consumed_kWh'], df_final['pred_tabICL'])
print(f"#########################################OFFline TabICL RMSE: {score_nrmse:.4f}")


score_nmae = nmae(df_final['Consumed_kWh'], df_final['pred_tabICL'])
print(f"#########################################OFFline TabICL RMSE: {score_nmae:.4f}")

