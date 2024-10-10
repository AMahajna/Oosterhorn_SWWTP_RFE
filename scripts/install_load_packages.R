###############################################################################
## Function to install and load required packages
install_and_load <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

## List of required packages
required_packages <- c(
  "FeatureTerminatoR", "ggcorrplot", "DataExplorer", "caret", "mlbench", 
  "randomForest", "tidyverse", "readxl", "psych", "heatmaply", "GGally", 
  "corrplot", "rio", "car", "FSelector", "data.table", "dplyr", "lubridate", 
  "reshape2", "ggplot2", "cowplot","gridExtra"
)

## Install and load the packages
install_and_load(required_packages)