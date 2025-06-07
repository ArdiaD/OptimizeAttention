# File that generate all results in the paper
# !!! Will take several days !!!
rm(list = ls())
sessionInfo()

# Simple setup; generate some results with 1 day of computations
rm(list = ls())
source("01_replicate_epu.R")
source("02_plot_epu.R")
rm(list = ls())
source("03_forecast_inflation.R")
rm(list = ls())
source("04_nowcast_inflation.R")
rm(list = ls())
source("05_measure_performance.R")
rm(list = ls())
source("06_analyze_topics.R")
rm(list = ls())
source("07_plot_sentiment.R")

# Full set of results; !!! takes days !!!
rm(list = ls())
source("01_replicate_epu.R")
source("02_plot_epu.R")
rm(list = ls())
h = 1; lambda_1 = 1; use_lag = 1
source("03_forecast_inflation.R")
rm(list = ls())
h = 3; lambda_1 = 1; use_lag = 1
source("03_forecast_inflation.R")
rm(list = ls())
h = 6; lambda_1 = 1; use_lag = 1
source("03_forecast_inflation.R")
rm(list = ls())
h = 1; lambda_1 = 0; use_lag = 1
source("03_forecast_inflation.R")
rm(list = ls())
h = 3; lambda_1 = 0; use_lag = 1
source("03_forecast_inflation.R")
rm(list = ls())
h = 6; lambda_1 = 0; use_lag = 1
source("03_forecast_inflation.R")

rm(list = ls())
h = 1; lambda_1 = 1; use_lag = 0
source("03_forecast_inflation.R")
rm(list = ls())
h = 3; lambda_1 = 1; use_lag = 0
source("03_forecast_inflation.R")
rm(list = ls())
h = 6; lambda_1 = 1; use_lag = 0
source("03_forecast_inflation.R")

rm(list = ls())
h = 1; lambda_1 = 0; use_lag = 0
source("03_forecast_inflation.R")
rm(list = ls())
h = 3; lambda_1 = 0; use_lag = 0
source("03_forecast_inflation.R")
rm(list = ls())
h = 6; lambda_1 = 0; use_lag = 0
source("03_forecast_inflation.R")

rm(list = ls())
h = 1; lambda_1 = 0; use_lag = 0
source("04_nowcast_inflation.R")
rm(list = ls())
h = 1; lambda_1 = 1; use_lag = 0
source("04_nowcast_inflation.R")
rm(list = ls())
h = 1; lambda_1 = 0; use_lag = 0
source("04_nowcast_inflation.R")
rm(list = ls())
h = 1; lambda_1 = 0; use_lag = 1
source("04_nowcast_inflation.R")

rm(list = ls())
source("05_measure_performance.R")
rm(list = ls())
source("06_analyze_topics.R")
rm(list = ls())
source("07_plot_sentiment.R")