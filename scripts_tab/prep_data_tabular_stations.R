source('scripts/pretraitements.R')
source('scripts/modelisation.R')
library(argparser)

# inputs
p <- arg_parser("Script modèles stations VE")

p <- add_argument(p, "--seed", help="Seed pour la génération des fenêtres", default=40, type="integer")
p <- add_argument(p, "--input", help="Chemin du dataset d'entrée", default="Data/dataset_address_main.csv")
p <- add_argument(p, "--output", help="Chemin (prefixe) du fichier de sortie", default="results/results_stations")
p <- add_argument(p, "--parallel", help="Activer le calcul parallèle", default=TRUE, type="logical")
p <- add_argument(p, "--type", help="Type de modèle : global, local ou all", default="global")
p <- add_argument(p, "--lags", help="Inclure les variables de lags", default=TRUE, type="logical")
p <- add_argument(p, "--max_pct_na", help="Supprimer les séries qui ont trop de NA du train", default=1, type="numeric")
p <- add_argument(p, "--nbr_series_train", help="Réduire le nombre de séries dans le train",  default=NULL, type="integer")

argv <- parse_args(p)

seed_windows  <- argv$seed
raw_data_path <- argv$input
output_path   <- argv$output
bool_parallel <- argv$parallel
type          <- argv$type
lags          <- argv$lags
max_pct_na    <- argv$max_pct_na
nbr_series_train <- argv$nbr_series_train

print(nbr_series_train)

# Fenêtres 
my_windows <- generate_rolling_windows(SEED = seed_windows) 

# Dataset
dataset_stations <- data.table::fread(raw_data_path) %>% 
  mutate(Date = as.Date(Date)) %>% preparation_data_regions_stations(group_col=Station.ID)

write.table(dataset_stations, file=paste0('Data/Seed_', seed_windows, '/dataset_stations_for_tab.csv'), quote = F, row.names = F, sep=';')
# names(dataset_stations)
# 
# 
# dataset_stations[12399+1,]%>%summary
# dataset_stations[12399,]
# 
# dim(dataset_stations)
# nrow(dataset_stations)-1543990
# #1543990
