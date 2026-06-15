
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
import tam as ta

# 
import sys
#import os
#sys.path.append('/Users/Yannig/Documents/These_Guillaume_P/Appli_VE/Appli-VE-main/preprocessed_data/')
os.chdir('/Users/Yannig/Documents/These_Guillaume_P/Appli_VE/Appli-VE-main/preprocessed_data/')

#Import torch
import torch.nn as nn
import torch.optim
import torch
#Time series
from statsmodels.tsa.arima.model import ARIMA
#Metric
from sklearn.metrics import root_mean_squared_error as rmse

def nrmse(actual, predicted):
    # Convertit en tableaux numpy pour pouvoir faire des opérations vectorielles
    actual = np.asarray(actual)
    predicted = np.asarray(predicted)
    
    # rmse <- sqrt(mean((actual - predicted)^2, na.rm = TRUE))
    rmse = np.sqrt(np.nanmean((actual - predicted) ** 2))
    
    # n_rmse <- rmse / mean(actual, na.rm = TRUE)
    n_rmse = rmse / np.nanmean(actual)
    
    return n_rmse


def nmae(actual, predicted):
    # Convertit en tableaux numpy pour les calculs vectoriels
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


###########################################################################################
#######Test sur une window
###########################################################################################


# Historical data used to train the models
Data0 = pd.read_csv('dataset_scotland_for_tab.csv')
Data0['Date'] = pd.to_datetime(Data0['Date'])
n0 = len(Data0)
print(n0)

Window = pd.read_csv('mes_fenetres/window_1.csv')
Window['Date'] = pd.to_datetime(Window['Date'])

df_fusionne = pd.merge(Data0, Window, on='Date', how='right')

###covariate types
cat_cols = ['Weekday_Holiday_BIS']
num_cols = ['Posan', 'Lag1', 'Lag7', 'Lag_Mean_1_7', 'max_paid_cp_30days', 'n_cp_paid_lag1',
            'tmpf_max']
target_col = 'Consumed_kWh'

useful_col= cat_cols + num_cols + ['Date', 'type', target_col]
df_propre = df_fusionne[useful_col].dropna()




#gam_formula <- "Consumed_kWh ~  Weekday_Holiday_BIS + 
#                                s(Posan, bs='cc') + 
#                                Lag1 + 
#                                Lag7 + 
#                                Lag_Mean_1_7 + 
#                                n_cp_paid_lag1:max_paid_cp_30days + 
#                                s(tmpf_max)"


###0
tam_formula = "Consumed_kWh ~ s(Posan, k=10, extrapolate='continue') + l(Lag1) + l(Lag7) + " \
"l(Lag_Mean_1_7) + te(l(n_cp_paid_lag1),c(max_paid_cp_30days)) + s(tmpf_max, k=10, extrapolate='continue')"

tam_model = ta.StaticTAM(formula = tam_formula, date_col="Date").fit(df_propre)
df_propre_test = df_propre.loc[~(df_propre['type'] =='train'), ]
pred_tam = tam_model.predict(df_propre_test)['EstimatedConsumed_kWh']

score_nrmse = nrmse(df_propre_test['Consumed_kWh'], pred_tam)
print(f"#########################################OFFline TabICL RMSE: {score_nrmse:.4f}")

score_nmae = nmae(df_propre_test['Consumed_kWh'], pred_tam)
print(f"#########################################OFFline TabICL NMAE: {score_nmae:.4f}")

###1
tam_formula = "Consumed_kWh ~ w(Posan, n_scales=5, extrapolate='continue') + l(Lag1) + l(Lag7) + " \
"l(Lag_Mean_1_7) + l(n_cp_paid_lag1) + c(max_paid_cp_30days) + w(tmpf_max, k=10, extrapolate='continue')"
tam_model = ta.StaticTAM(formula = tam_formula, date_col="Date").fit(df_propre)
df_propre_test = df_propre.loc[~(df_propre['type'] =='train'), ]
pred_tam = tam_model.predict(df_propre_test)['EstimatedConsumed_kWh']

score_nrmse = nrmse(df_propre_test['Consumed_kWh'], pred_tam)
print(f"#########################################OFFline TabICL RMSE: {score_nrmse:.4f}")
