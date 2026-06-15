source('scripts/pretraitements.R')
source('scripts/modelisation.R')
library(argparser)
# inputs
p <- arg_parser("Script modèle scotland VE")

p <- add_argument(p, "--seed",      help="Seed pour la génération des fenêtres", default=40, type="integer")
p <- add_argument(p, "--input",     help="Chemin du dataset d'entrée",           default="Data/dataset_scotland_main.csv")
p <- add_argument(p, "--output",    help="Chemin (prefixe) du fichier de sortie", default="results/results_scotland.RDS")
p <- add_argument(p, "--parallel",  help="Activer le calcul parallèle",          default=TRUE, type="logical")

argv <- parse_args(p)

seed_windows  <- argv$seed
raw_data_path <- argv$input
output_path   <- argv$output
bool_parallel <- argv$parallel

# Fenêtres 
my_windows <- generate_rolling_windows(SEED = seed_windows)  

# Dataset
dataset_scotland <- read.csv(raw_data_path) %>% 
  mutate(Date = as.Date(Date)) %>% preparation_data_scotland()


write.csv(dataset_scotland, file=paste0('Data/Seed_', seed_windows,'/dataset_scotland_for_tab.csv'), quote = F)
names(dataset_scotland)

dir.create(paste0('Data/Seed_', seed_windows, '/my_windows'), showWarnings = TRUE, recursive = FALSE, mode = "0777")

for(j in c(1:length(my_windows)))
{
  write.csv(my_windows[[j]], file=paste0('Data/Seed_', seed_windows, '/my_windows', '/window_',j, '.csv'), quote = F, row.names = F)
}






