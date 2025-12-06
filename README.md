# Space-Time Gaussian Process Approximations

This repository provides the R code for the simulation studies and real-world applications presented in the paper *"Scalable non-separable spatio-temporal Gaussian Process Models for large-scale Climate Data."*  

The iterative methods for full-scale approximations are implemented in the **GPBoost** package, available here: [https://github.com/fabsig/GPBoost](https://github.com/fabsig/GPBoost).

## Repository Structure

- **`Data`**  
  Provides scripts for generating data and the actual data for both the simulation studies and real-world experiments.

- **`Simulation_Studies`**  
  Contains scripts to run the simulation experiments. Refer to the `README.md` file in the `simulation_studies` folder for a detailed description of the various simulations and how to execute them.

- **`Real_World_Application`**  
  Includes scripts for the real-world application on daytime land surface temperature data.
