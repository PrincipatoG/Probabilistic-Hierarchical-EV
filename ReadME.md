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
pip install -r python_requirements.txt
```

## Data collection

The data scrapping rely on a (modified) code from the following [github repository](https://github.com/djordjebatic/GridCharge).

For convenience and reproducibility purpose, the resulting pre-processed dataset is available at **mettre les .csv en ligne**.

>📋 The complete pipeline used to generate the data strongly depends on the configuration of the user: please contact use if you need the associated code.

## Running Experiments

To perform the complete experiment, run this command:

```
Rscript run_all.R
```

>📋 This code can be run on a personal computer but actually requires a strong memory capacity (my 36 GB Mac almost took fire at some point). A separate run of a each of its item with well-chosen and device dependant parameters is thus recommended.

## Results

We copy here the main experimental results from the article:

Validity rates at two coverage levels: the nominal level $1-\alpha=90\%$ and a relaxed level $1-2\alpha=80\%$. 
Results are grouped by aggregation level and conformal procedure, and reported as $\operatorname{mean} \pm 1.96 \times \operatorname{standard\ error}$.

| Node type | Level | CP | CP-MNR | CP-MNR Nested$^{\star}$ | Adaptive CP-MNR |
|:---|:---:|---:|---:|---:|---:|
| National | $1-\alpha$ | 40.0% $\pm$ 16.3% | 40.0% $\pm$ 16.3% | **90.0% $\pm$ 10.0%** | 0.0% $\pm$ 0.0% |
|  | $1-2\alpha$ | **100.0% $\pm$ 0.0%** | **100.0% $\pm$ 0.0%** | **100.0% $\pm$ 0.0%** | **100.0% $\pm$ 0.0%** |
| Regional | $1-\alpha$ | 51.9% $\pm$ 1.2% | 51.9% $\pm$ 1.2% | **71.2% $\pm$ 1.2%** | 47.8% $\pm$ 0.9% |
|  | $1-2\alpha$ | 98.8% $\pm$ 0.5% | 98.8% $\pm$ 0.5% | 99.1% $\pm$ 0.5% | **100.0% $\pm$ 0.0%** |
| Station | $1-\alpha$ | 38.1% $\pm$ 0.4% | 56.2% $\pm$ 0.4% | **73.6% $\pm$ 0.6%** | 57.6% $\pm$ 0.6% |
|  | $1-2\alpha$ | 67.7% $\pm$ 0.3% | 93.6% $\pm$ 0.3% | 90.8% $\pm$ 0.4% | **98.0% $\pm$ 0.1%** |
