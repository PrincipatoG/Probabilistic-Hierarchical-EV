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

# Import torch et configuration pour MPS (Metal Performance Shaders)
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

# Détection de la puce Apple Silicon (MPS)
if torch.backends.mps.is_available():
    device = torch.device("mps")
    print("🚀 GPU MPS détecté. PyTorch utilisera la puce Apple Silicon.")
else:
    device = torch.device("cpu")
    print("⚠️ GPU MPS non détecté. Utilisation du CPU.")

# Time series
from statsmodels.tsa.arima.model import ARIMA
# Metric
from sklearn.metrics import root_mean_squared_error as rmse

def nrmse(actual, predicted):
    # Convertit en tableaux numpy pour pouvoir faire des opérations vectorielles
    actual = np.asarray(actual)
    predicted = np.asarray(predicted)
    
    # rmse <- sqrt(mean((actual - predicted)^2, na.rm = TRUE))
    rmse_val = np.sqrt(np.nanmean((actual - predicted) ** 2))
    
    # n_rmse <- rmse / mean(actual, na.rm = TRUE)
    n_rmse_val = rmse_val / np.nanmean(actual)
    
    return n_rmse_val

def nmae(actual, predicted):
    # Convertit en tableaux numpy pour les calculs vectoriels
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

# 1. Chargement des données 
print("Chargement des données historiques...")
Data0 = pd.read_csv('dataset_stations_for_tab.csv', on_bad_lines='warn', sep=';')
Data0['Date'] = pd.to_datetime(Data0['Date'])
print(f"Nombre de lignes dans Data0: {len(Data0)}")

# Définition des colonnes
cat_cols = ['Weekday_Holiday_BIS', 'Region', 'Station.ID']
num_cols = ['Posan', 'tmpf_max', 'Latitude', 'relh_mean', 'max_cp_30days', 'taux_cp_paid_max30',
            'taux_lent', 'taux_rapide', 'taux_accelere', 'Lag_Mean_1_7', 'n_cp_lag1', 'Lag7', 'Lag1']

target_col = 'Consumed_kWh'
useful_col = cat_cols + num_cols + ["type", target_col]

# Paramètres de l'ensembling
num_models = 5      # Le nombre d'experts
sample_size = 10000 # Taille du sampling pour le bagging
batch_size = 2**10 # Taille des chunks d'inférence (16384)

# Liste pour stocker les résultats de toutes les fenêtres
all_results = []

# 2. Boucle sur les fenêtres de 1 à 26
for window_id in range(1, 15):
    window_path = f"my_windows/window_{window_id}.csv"

    # Vérification si le fichier existe
    if not os.path.exists(window_path):
        print(f"Fichier {window_path} introuvable, passage à la suite.")
        continue

    print(f"\n=======================================================")
    print(f"--- Traitement de la Window {window_id} ---")

    # Chargement et fusion
    Window = pd.read_csv(window_path)
    Window["Date"] = pd.to_datetime(Window["Date"])

    df_fusionne = pd.merge(Data0, Window, on="Date", how="right")

    # Nettoyage des valeurs manquantes
    # On garde 'Date' pour l'export final
    df_propre = df_fusionne[["Date"] + useful_col].dropna().copy()

    # Split Train / Test
    is_train = df_propre["type"] == "train"
    df_train_full = df_propre.loc[is_train].copy()
    
    is_train = df_fusionne["type"] == "train"
    X_test = df_fusionne.loc[~(df_fusionne['type'] =='train'), cat_cols + num_cols].copy()
    y_test = df_fusionne.loc[~(df_fusionne['type'] =='train'), target_col].copy()


    # Si pas de données de test pour cette fenêtre, on passe à la suivante
    if X_test.empty:
        print(f"Pas de lignes 'test' pour la Window {window_id}, passage à la suite.")
        continue

    X_test.columns = [str(c) for c in X_test.columns]

    # --- DÉBUT DE L'ENSEMBLING POUR CETTE FENÊTRE ---
    all_predictions = []
    
    for i in range(num_models):
        print(f"  -> Entraînement du modèle {i+1}/{num_models} (Window {window_id})...")
        
        # 1. Échantillonnage aléatoire
        train_sample = df_train_full.sample(n=sample_size, random_state=42 + i)
        X_train_sub = train_sample[cat_cols + num_cols].copy()
        y_train_sub = train_sample[target_col].copy()
        
        X_train_sub.columns = [str(c) for c in X_train_sub.columns]
        
        # 2. Initialisation du modèle
        reg = TabICLRegressor(verbose=False, random_state=30 + i, n_estimators=1, device=device)
        
        # 3. Entraînement
        reg.fit(X_train_sub, y_train_sub.values.ravel())
        
        # 4. Inférence (par batch pour éviter les pb de RAM)
        pred_tabICL = []
        for j in tqdm(range(0, len(X_test), batch_size), desc="     Prédiction", leave=False):
            X_batch = X_test.iloc[j:j+batch_size]
            batch_preds = reg.predict(X_batch)
            pred_tabICL.extend(batch_preds)
            
        # Ajout des prédictions à notre ensemble
        all_predictions.append(pred_tabICL)
        
        # Vidage du cache GPU entre chaque modèle
        if torch.backends.mps.is_available():
            torch.mps.empty_cache()

    # 5. Ensembling : Moyenne des prédictions pour cette fenêtre
    ensemble_preds = np.mean(all_predictions, axis=0)
    # --- FIN DE L'ENSEMBLING ---

    # On prépare le DataFrame final pour cette fenêtre
    is_train = df_fusionne["type"] == "train"  
    df_test = df_fusionne.loc[~is_train].copy()
    df_test["pred_tabICL"] = ensemble_preds
    df_test["window_id"] = window_id

    # Organisation des colonnes pour le rendu
    cols_finales = ["window_id", "Date"] + useful_col + ["pred_tabICL"]
    all_results.append(df_test[cols_finales])

    nom_fichier_export_provisoire = f"../../Output_ponctual/Seed_{seed}/results_stations_tabICL_window_{window_id}.csv"
    df_test[cols_finales].to_csv(nom_fichier_export_provisoire, index=False, sep=';')



# 3. Fusion de toutes les fenêtres et export en CSV
if all_results:
    df_final = pd.concat(all_results, ignore_index=True)

    # Exportation en fichier CSV
    nom_fichier_export = f"../../Output_ponctual/Seed_{seed}/results_stations_tabICL.csv"
    
    # Création du dossier 'results' s'il n'existe pas déjà
    os.makedirs(os.path.dirname(nom_fichier_export), exist_ok=True)
    
    # Export avec séparateur point-virgule
    df_final.to_csv(nom_fichier_export, index=False, sep=';')

    print(f"\n### Traitement terminé ! ###")
    print(f"Le fichier a été exporté sous : {nom_fichier_export}")
    print(f"Nombre total de lignes de test enregistrées : {len(df_final)}")
else:
    print("\nAucune donnée de test n'a pu être traitée.")
