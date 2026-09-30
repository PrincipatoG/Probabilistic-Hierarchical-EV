# Probabilistic Forecasting of Electric Vehicle Charging in a Time-varying Hierarchical Infrastructure

This repository is the official implement of the article *Probabilistic Forecasting of Electric Vehicle Charging in a Time-varying Hierarchical Infrastructure*.

## Requirements

To install the R packages:

```bash
Rscript install_packages.R
```

>📋  The experiments are primarly run under [R version 4.4.3](https://cran.r-project.org/bin/windows/base/old/4.4.3/).

The experiments also rely on [python version 3.8.18](https://www.python.org/downloads/release/python-3818/).

To install the python packages:

```bash
pip install -r requirements.txt
```

## Data collection

The data scrapping rely on a (modified) code from the following [github repository](https://github.com/djordjebatic/GridCharge).

For convenience and reproducibility purpose, the resulting pre-processed dataset is available at **mettre les .csv en ligne**.

>📋 The complete pipeline used to generate the data strongly depends on the configuration of the user: please contact use if you need the associated code.

## Running Experiments

To perform the complete experiment, run this command:

```
sbatch launcher.sh
```

>📋 This code can be run on a personal computer. The "config" parameter can be adapted to correspond to one of the six configurations considered in the article. The arguments of the "config" parameter in the above example are depth = 3, n = 12, T = $10^6$, N= $10^3$ (not used here) and 1, the configuration index. The complete experiment can be run using a "for loop" on the simulation indices in ./R/study.R (but we prefer the parallelized setup described in the following section).

## Results

We copy here the main experimental results from the article:

### SCP for joint coverage based on ellipsoidal sets

|     Matrix A      | Config. |     Alg. (1)   |     Alg. (2)       |
|-------------------|---------|----------------|--------------------|
| $\widehat{\Sigma}^{-1}$ |    1    | 17.4 ± 0.6     | **16.5 ± 0.55**    |
| $\widehat{\Sigma}^{-1}$ |    2    | 3.36 ± 0.12    | **3.19 ± 0.11**    |

### Component-wise SCP

| Config. | Direct         | OLS            | WLS              | Combi            | MinT             |
|---------|----------------|----------------|------------------|------------------|------------------|
| 1       | 876 ± 254      | 787 ± 226      | 322 ± 131        | 364 ± 101        | **216 ± 47**     |
| 2       | 871 ± 253      | 753 ± 216      | 308 ± 116        | 361 ± 92         | **246 ± 51**     |
| 3       | 3032 ± 467     | 2954 ± 455     | 1869 ± 377       | 1758 ± 395       | **1502 ± 578**   |
| 4       | 3036 ± 479     | 2901 ± 458     | 1581 ± 340       | 1604 ± 349       | **1404 ± 571**   |
| 5       | 10424 ± 885    | 10358 ± 880    | **8861 ± 853**   | 9664 ± 850       | 10613 ± 918      |
| 6       | 10621 ± 889    | 10460 ± 875    | **7673 ± 785**   | 9068 ± 806       | 10503 ± 905      |
