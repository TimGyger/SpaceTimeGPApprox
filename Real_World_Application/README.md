# Real World Application

These scripts are used to reproduce the real world experiments in the paper.

- **`Temperature.ipynb`**  
  Reproduces the experiments in *Section: Short-term prediction of daily maximum temperature*, including rolling five-day-ahead forecasts for 2025, monthly parameter re-estimation, and evaluation using RMSE and CRPS under temporal and spatio-temporal prediction at withheld stations.

- **`Precipitation.ipynb`**  
  Reproduces the experiments in *Section: Short-term prediction of precipitation*, using the exponential space–time covariance model for rainfall and the same operational forecasting framework.

  ## Forecasting Setup

- Expanding-window training starting from 2024.
- Monthly covariance parameter re-estimation.
- 25% of stations withheld for spatial generalization.
- Benchmarks include persistence and fixed-effects regression.

