source('Model_training/R/pretraitements.R')
source('Model_training/R/modelisation.R')

library(argparser)

# inputs
p <- arg_parser("Regional EV models script")

p <- add_argument(p, "--seed",      help="Seed for generating windows", default=40, type="integer")
p <- add_argument(p, "--input",     help="Input dataset path",           default="Data/dataset_region_main.csv")
p <- add_argument(p, "--output",    help="Output file path prefix", default="results/results_regions")
p <- add_argument(p, "--parallel",  help="Enable parallel computation",          default=TRUE, type="logical")
p <- add_argument(p, "--type",      help="Model type: global, local, or all", default="global")
p <- add_argument(p, "--transform", help="Apply the log1p transformation or normalization",    default='log')
p <- add_argument(p, "--lags",      help="Include lag variables",        default=TRUE, type="logical")

argv <- parse_args(p)

seed_windows  <- argv$seed
raw_data_path <- argv$input
output_path   <- argv$output
bool_parallel <- argv$parallel
type          <- argv$type
transform     <- tolower(argv$transform)
lags          <- argv$lags

# Windows
my_windows <- generate_rolling_windows(SEED = seed_windows) 

# Dataset
dataset_regions <- read.csv(raw_data_path) %>%
  mutate(Date = as.Date(Date)) %>% preparation_data_regions_stations(group_col=Region)

write.csv(dataset_regions, file=paste0('Data/Seed_', seed_windows,'/dataset_regions_for_tab.csv'), quote = F)
# names(dataset_regions)
# 
# 
# # for(j in c(1:length(my_windows)))
# # {
# #   write.csv(my_windows[[j]], file=paste0('Data/my_windows/window_',j, '.csv'), quote = F, row.names = F)
# # }


