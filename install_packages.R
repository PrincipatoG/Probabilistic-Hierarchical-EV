# install_packages.R

# Install remotes if necessary
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

packages <- c(
  "argparser",
  "broom",
  "cowplot",
  "data.table",
  "dplyr",
  "dtplyr",
  "fastmatrix",
  "forcats",
  "forecast",
  "ggplot2",
  "ggspatial",
  "kableExtra",
  "knitr",
  "lars",
  "lubridate",
  "magrittr",
  "patchwork",
  "purrr",
  "ranger",
  "readr",
  "RSpectra",
  "scales",
  "sf",
  "slider",
  "stringr",
  "tidyr",
  "xgboost"
)

versions <- c(
  "0.7.3",
  "1.0.8",
  "1.2.0",
  "1.17.2",
  "1.1.4",
  "1.3.1",
  "0.6.6",
  "1.0.0",
  "9.0.0",
  "3.5.2",
  "1.1.10",
  "1.4.0",
  "1.50",
  "1.3",
  "1.9.4",
  "2.0.3",
  "1.3.2",
  "1.0.4",
  "0.18.0",
  "2.1.5",
  "0.16.2",
  "1.4.0",
  "1.0.23",
  "0.3.3",
  "1.5.1",
  "1.3.1",
  "3.2.1.1"
)

for (i in seq_along(packages)) {
  
  package <- packages[i]
  version <- versions[i]
  
  installed <- requireNamespace(
    package,
    quietly = TRUE
  )
  
  if (
    !installed ||
    as.character(packageVersion(package)) != version
  ) {
    
    message(
      "Installing ",
      package,
      " version ",
      version
    )
    
    remotes::install_version(
      package = package,
      version = version,
      repos = "https://cloud.r-project.org"
    )
    
  } else {
    
    message(
      package,
      " ",
      version,
      " already installed."
    )
  }
}