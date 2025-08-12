################################################################################
### The Second Competition on Spatial Statistics for Large Datasets
################################################################################

# Abdulah, Sameh, et al. "The second competition on spatial statistics for large datasets." arXiv preprint arXiv:2211.03119 (2022).

#######
## Packages
#######

library(gpboost)
library(scoringRules)

data_solution <- read.csv("https://raw.githubusercontent.com/TimGyger/SpaceTimeGPApprox/refs/heads/main/Data/2a-solutions.csv")

vec_vecchia_RMSE <- rep(0,9) 
vec_vecchia_corr_RMSE <- rep(0,9)
vec_fitc_RMSE <- rep(0,9)
vec_fitc_grid_RMSE <- rep(0,9)
vec_vif_RMSE <- rep(0,9) 

vec_vecchia_CRPS <- rep(0,9) 
vec_vecchia_corr_CRPS <- rep(0,9)
vec_fitc_CRPS <- rep(0,9)
vec_fitc_grid_CRPS <- rep(0,9)
vec_vif_CRPS <- rep(0,9) 

vec_vecchia_time <- rep(0,9) 
vec_vecchia_corr_time <- rep(0,9)
vec_fitc_time <- rep(0,9)
vec_fitc_grid_time <- rep(0,9)
vec_vif_time <- rep(0,9) 

for (i in 1:9) {
  data_train <- read.csv(paste0("https://raw.githubusercontent.com/TimGyger/SpaceTimeGPApprox/refs/heads/main/Data/2a_",i,"_train.csv"))
  data_test <- read.csv(paste0("https://raw.githubusercontent.com/TimGyger/SpaceTimeGPApprox/refs/heads/main/Data/2a_",i,"_test.csv"))
  coords_train <- cbind(data_train$t,data_train$x,data_train$y)
  y_train <- data_train$z
  
  coords_test <- cbind(data_test$t,data_test$x,data_test$y)
  y_test <- data_solution$z3
  
  # Vecchia (euclidean)
  ptm <- proc.time()
  Linear_gp_model <- fitGPModel(gp_coords =  coords_train,y = y_train, likelihood = "gaussian", 
                                cov_function = "space_time_gneiting", vecchia_ordering = "random", cov_fct_shape = 1.5, matrix_inversion_method = "cholesky", seed = 1, 
                                num_neighbors = 30, gp_approx = "vecchia",params = list(trace = T,estimate_cov_par_index=c(1,1,1,1,1,0,1,1)))#,optimizer_coef = "gradient_descent",optimizer_cov = "gradient_descent"))
  
  pred_Linear_Model <- predict(Linear_gp_model,
                               gp_coords_pred = coords_test, 
                               y = y_train, predict_var = T)
  pred_mu <- pred_Linear_Model$mu
  pred_var <- pred_Linear_Model$var
  vec_vecchia_time[i] <- proc.time() - ptm
  vec_vecchia_RMSE[i] <- sqrt(mean((pred_t-(y_test))^2))
  vec_vecchia_CRPS[i] <- mean(crps_norm(y = y_test,mean = pred_mu,sd = sqrt(pred_var)))
  
  # Vecchia (correlation)
  ptm <- proc.time()
  Linear_gp_model <- fitGPModel(gp_coords =  coords_train,y = y_train, likelihood = "gaussian", 
                                cov_function = "space_time_gneiting", vecchia_ordering = "random", cov_fct_shape = 1.5, matrix_inversion_method = "cholesky", seed = 1, 
                                num_neighbors = 30, gp_approx = "vecchia_correlation_based",params = list(trace = T,estimate_cov_par_index=c(1,1,1,1,1,0,1,1)))#,optimizer_coef = "gradient_descent",optimizer_cov = "gradient_descent"))
  
  pred_Linear_Model <- predict(Linear_gp_model,
                               gp_coords_pred = coords_test, 
                               y = y_train, predict_var = T)
  pred_mu <- pred_Linear_Model$mu
  pred_var <- pred_Linear_Model$var
  vec_vecchia_corr_time[i] <- proc.time() - ptm
  vec_vecchia_corr_RMSE[i] <- sqrt(mean((pred_t-(y_test))^2))
  vec_vecchia_corr_CRPS[i] <- mean(crps_norm(y = y_test,mean = pred_mu,sd = sqrt(pred_var)))
  
  # FITC (kmeans++)
  ptm <- proc.time()
  Linear_gp_model <- fitGPModel(gp_coords =  coords_train,y = y_train, likelihood = "gaussian", 
                                cov_function = "space_time_gneiting",cov_fct_shape = 1.5, matrix_inversion_method = "cholesky", seed = 1, 
                                ind_points_selection = "kmeans++", num_ind_points = 500, gp_approx = "fitc",params = list(trace = T,estimate_cov_par_index=c(1,1,1,1,1,0,1,1)))#,optimizer_coef = "gradient_descent",optimizer_cov = "gradient_descent"))
  
  pred_Linear_Model <- predict(Linear_gp_model,
                               gp_coords_pred = coords_test, 
                               y = y_train, predict_var = T)
  pred_mu <- pred_Linear_Model$mu
  pred_var <- pred_Linear_Model$var
  vec_fitc_time[i] <- proc.time() - ptm
  vec_fitc_RMSE[i] <- sqrt(mean((pred_t-(y_test))^2))
  vec_fitc_CRPS[i] <- mean(crps_norm(y = y_test,mean = pred_mu,sd = sqrt(pred_var)))
  
  # FITC (kmeans++-grid)
  ptm <- proc.time()
  Linear_gp_model <- fitGPModel(gp_coords =  coords_train,y = y_train, likelihood = "gaussian", 
                                cov_function = "space_time_gneiting",cov_fct_shape = 1.5, matrix_inversion_method = "cholesky", seed = 1, 
                                ind_points_selection = "kmedoids", num_ind_points = 500, gp_approx = "fitc",params = list(trace = T,estimate_cov_par_index=c(1,1,1,1,1,0,1,1)))#,optimizer_coef = "gradient_descent",optimizer_cov = "gradient_descent"))
  
  pred_Linear_Model <- predict(Linear_gp_model,
                               gp_coords_pred = coords_test, 
                               y = y_train, predict_var = T)
  pred_mu <- pred_Linear_Model$mu
  pred_var <- pred_Linear_Model$var
  vec_fitc_grid_time[i] <- proc.time() - ptm
  vec_fitc_grid_RMSE[i] <- sqrt(mean((pred_t-(y_test))^2))
  vec_fitc_grid_CRPS[i] <- mean(crps_norm(y = y_test,mean = pred_mu,sd = sqrt(pred_var)))
  
  # VIF
  ptm <- proc.time()
  Linear_gp_model <- fitGPModel(gp_coords =  coords_train,y = y_train, likelihood = "gaussian", 
                                cov_function = "space_time_gneiting",cov_fct_shape = 1.5, matrix_inversion_method = "cholesky", seed = 1, 
                                num_neighbors = 30,ind_points_selection = "kmedoids", num_ind_points = 500, gp_approx = "vif_correlation_based",vecchia_ordering = "random",
                                params = list(trace = T,estimate_cov_par_index=c(1,1,1,1,1,0,1,1)))#,optimizer_coef = "gradient_descent",optimizer_cov = "gradient_descent"))
  
  pred_Linear_Model <- predict(Linear_gp_model,
                               gp_coords_pred = coords_test, 
                               y = y_train, predict_var = T)
  pred_mu <- pred_Linear_Model$mu
  pred_var <- pred_Linear_Model$var
  vec_vif_time[i] <- proc.time() - ptm
  vec_vif_RMSE[i] <- sqrt(mean((pred_t-(y_test))^2))
  vec_vif_CRPS[i] <- mean(crps_norm(y = y_test,mean = pred_mu,sd = sqrt(pred_var)))
  
  
}





