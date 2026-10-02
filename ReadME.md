# Probabilistic Forecasting of EV Charging in a Time-varying Hierarchical Infrastructure

This repository is the official implementation of the article *Probabilistic Forecasting of Electric Vehicle Charging in a Time-varying Hierarchical Infrastructure*.

## Requirements

To install the R packages:

```bash
Rscript install_packages.R
```

>📋  The experiments are primarly run under [R version 4.4.3](https://cran.r-project.org/bin/windows/base/old/4.4.3/).

The experiments also rely on [Python version 3.11](https://www.python.org/downloads/release/python-31116/).

To install the python packages:

```bash
pip install -r python_requirements.txt
```

## Data Collection

The data scrapping rely on a (modified) code from the following [github repository](https://github.com/djordjebatic/GridCharge).

For convenience and reproducibility purpose, the resulting pre-processed datasets can be found here: [**Download the dataset**](https://github.com/PrincipatoG/AOAS-EV/releases/latest/download/Data.zip).

>📋 The complete pipeline used to generate the data strongly depends on the configuration of the user: please contact me if you need the associated code.

## Running Experiments

To perform the complete experiment, run this command:

```
Rscript run_all.R
```

>📋 **Computational requirement** The code can be run on a personal computer, but requires substantial RAM. For reference, the experiments were run on a Mac with 36 GB of RAM. At some point, the memory usage was close to the available capacity. A separate run of a each of its item with well-chosen and device dependant parameters is thus recommended.

### Experiments in Details

The experiments steps can be decomposed as follow:

1. Model training

📊 Code in R: RF, GAM and XGboost 
```
Rscript Model_training/R/run_scotland.R
Rscript Model_training/R/run_regions.R
Rscript Model_training/R/run_stations.R
```

🐍 Code in Python: TabICL
```
python3.11 Model_training/Python/tabICL_national.py
python3.11 Model_training/Python/tabICL_regions_GPU.py
python3.11 Model_training/Python/tabICL_stations_GPU.py
```

2. Forecast reconciliation (for point forecasting)
```
Rscript R/forecast_reconciliation.R
```

3. Conformal prediction with an integrated reconciliation step
```
Rscript R/run_conformal.R
```
The available options for `run_conformal.R` are:

| Argument             | Available options              | Description                         |
| -------------------- | ------------------------------ | ----------------------------------- |
| `--conformal_method` | `CP_MNR`, `CP_MNR_naive`, `Adaptive_CP_MNR`, `CP_MNR_Nested_star`    | Procedure to perform conformal prediction with a time-varying hierarchical structure                |
| `--national_model`   | `LOCAL_GAM`, `LOCAL_RF`, `LOCAL_XBG`, `tabICL`, `Combination`        | Forecasting model used at the national level |
| `--regional_model`   | `GLOBAL_GAM`, `GLOBAL_RF`, `GLOBAL_XGB`, `tabICL`, `Combination`     | Forecasting model used at the regional level (for all nodes) |
| `--station_model`    | `GLOBAL_GAM`, `GLOBAL_RF`, `GLOBAL_XGB`, `tabICL`, `Combination`     | Forecasting model used at the station level (for all nodes) |
| `--seed`             | Integer              | Seed used for the experiment                  |

For example, to run the **CP_MNR_Nested_star** conformal procedure with forecasts being the combination of all the forecasting methods:

```bash
  Rscript R/run_conformal.R ,
    --conformal_method "CP_MNR_Nested_star",
    --national_model "Combination",
    --regional_model "Combination",
    --station_model "Combination",
    --seed 1
```

### From Results to Figures

Four files can be used to get the plot used to illustrate the results of the experiments. 
Some post-treatment in LaTeX are made from the row ggplot2 figures to obtain *journal-ready* figures.
```bash
Rscript R/Figure_Example.R
Rscript R/Figure_multi_seed.R
Rscript R/Figure_multi_seed_cond.R
Rscript R/Figure_multi_seed_ponctual.R
```

## Results

We copy here the main experimental results from the article:

Validity rates, i.e., proportion of valid pairs (node, window) at two coverage levels: the nominal level 1 − α = 90% and a relaxed level 1 − 2α = 80%.
Results are grouped by aggregation level and conformal procedure, and reported as mean ± 1.96 × standard error.

| Node type | Level | CP | CP-MNR | CP-MNR Nested⋆ | Adaptive CP-MNR |
|:---|:---:|---:|---:|---:|---:|
| National | 1 − α | 40.0% ± 16.3% | 40.0% ± 16.3% | **90.0% ± 10.0%** | 0.0% ± 0.0% |
|  | 1 − 2α | **100.0% ± 0.0%** | **100.0% ± 0.0%** | **100.0% ± 0.0%** | **100.0% ± 0.0%** |
| Regional | 1 − α | 51.9% ± 1.2% | 51.9% ± 1.2% | **71.2% ± 1.2%** | 47.8% ± 0.9% |
|  | 1 − 2α | 98.8% ± 0.5% | 98.8% ± 0.5% | 99.1% ± 0.5% | **100.0% ± 0.0%** |
| Station | 1 − α | 38.1% ± 0.4% | 56.2% ± 0.4% | **73.6% ± 0.6%** | 57.6% ± 0.6% |
|  | 1 − 2α | 67.7% ± 0.3% | 93.6% ± 0.3% | 90.8% ± 0.4% | **98.0% ± 0.1%** |

## Additionnal Code

The code located in `./Additional/` can be used to generate the descriptive figures used in the article.
