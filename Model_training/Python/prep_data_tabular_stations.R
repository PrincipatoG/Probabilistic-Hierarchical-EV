source('Model_training/R/pretraitements.R')
source('Model_training/R/modelisation.R')

library(argparser)

# inputs
p <- arg_parser("Station EV models script")

p <- add_argument(p, "--seed", help="Seed for generating windows", default=40, type="integer")
p <- add_argument(p, "--input", help="Input dataset path", default="Data/dataset_address_main.csv")
p <- add_argument(p, "--output", help="Output file path prefix", default="results/results_stations")
p <- add_argument(p, "--parallel", help="Enable parallel computation", default=TRUE, type="logical")
p <- add_argument(p, "--type", help="Model type: global, local, or all", default="global")
p <- add_argument(p, "--lags", help="Include lag variables", default=TRUE, type="logical")
p <- add_argument(p, "--max_pct_na", help="Remove series with too many NAs from training", default=1, type="numeric")
p <- add_argument(p, "--nbr_series_train", help="Reduce the number of training series",  default=NULL, type="integer")

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

# Windows
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
