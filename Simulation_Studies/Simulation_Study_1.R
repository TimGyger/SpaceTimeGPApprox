#######################################################
### Packages
#######################################################

### Install/Load Packages
# Package names
packages <- c("Rcpp","RcppEigen", "fields", "gpboost")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

# Random Fields
# Function to check and install a package if not already installed
install_if_missing <- function(package, version = NULL) {
  if (!requireNamespace(package, quietly = TRUE)) {
    if (!is.null(version)) {
      remotes::install_version(package, version)
    } else {
      install.packages(package)
    }
  }
}

#######################################################
### Data
#######################################################
num_rep <- 10
n <- 15000
n_t_steps <- 30
t_steps <- seq(1, n_t_steps, length.out = n_t_steps)

Matern_Gneiting <- function(coord_i,coord_j,c,a,alpha,delta,beta){
  t_lag <- abs(coord_i[1]-coord_j[1])
  s_lag <- sqrt((coord_i[2]-coord_j[2])^2+(coord_i[3]-coord_j[3])^2)
  dist_t <- (a*t_lag^(2*alpha)+1)^(beta/2)
  dist_s <- c*s_lag
  dist <- dist_s/dist_t
  (1+dist)*exp(-dist)/(a*t_lag^(2*alpha)+1)^(beta+delta)
}


# Example parameters
sim_data <- function(c,a,alpha,delta,beta,seed_init){
  # Generate covariance matrix (careful: can be huge!)
  set.seed(10)
  space_coords <- matrix(runif(n/n_t_steps*2), ncol = 2)
  
  coords <- cbind(rep(t_steps,each = n/n_t_steps),rep(space_coords[,1],n_t_steps),rep(space_coords[,2],n_t_steps))
  
  n_pts <- nrow(coords)
  idx <- 1:n
  coords_small <- coords[idx, ]
  
  cov_mat <- matrix(0, length(idx), length(idx))
  for (i in 1:length(idx)) {
    for (j in i:length(idx)) {
      cov_ij <- Matern_Gneiting(coords_small[i, ], coords_small[j, ], c, a, alpha, delta, beta)
      cov_mat[i, j] <- cov_ij
      cov_mat[j, i] <- cov_ij  # symmetric
    }
  }
  
  
  cppFunction(depends = "RcppEigen", code = '
  Eigen::MatrixXd cholesky(const Eigen::MatrixXd &A) {
    return A.llt().matrixL();
  }
')
  
  sigma_marginal <- 1
  cov_mat <- cov_mat *sigma_marginal
  nugget <- 1e-2  # small positive value, adjust as needed
  # Define covariance matrix Sigma
  C <- cholesky(cov_mat) # faster than 'chol' R function
  
  n_pts <- nrow(cov_mat)
  set.seed(seed_init)
  z <- rnorm(n_pts)               # standard normals
  field <- C %*% z   
  y <- field + rnorm(length(field),0,sqrt(nugget))
  
  for (i in 1:n_t_steps) {
    start <- (i-1)*n/n_t_steps + 1 
    end <- start + n/n_t_steps -1
    quilt.plot(coords[start:end,2],coords[start:end,3],y[start:end])
  }
  
  return(cbind(coords,y))
}

#######################################################
### Experiments
#######################################################
para <- c(1e-2,1,0.5,20,0.4,1.5,0.4,0.2)
matrix_fsva <- matrix(0,num_rep,11)
matrix_vec <- matrix(0,num_rep,11)
matrix_vec_corr <- matrix(0,num_rep,11)
matrix_fitc_gird <- matrix(0,num_rep,11)
matrix_fitc <- matrix(0,num_rep,11)
matrixp_fsva <- matrix(0,num_rep,length(para))
matrixp_vec <- matrix(0,num_rep,length(para))
matrixp_vec_corr <- matrix(0,num_rep,length(para))
matrixp_fitc_gird <- matrix(0,num_rep,length(para))
matrixp_fitc <- matrix(0,num_rep,length(para))

for (ii in 1:num_rep) {
      data <- sim_data(20,0.5,0.4,0.4,0.2,ii)
      y <- data[,4]
      coords <- data[,1:3]
      coords_train <- coords[1:10000,]
      coords_test <- coords[10001:15000,]
      
      y_train <- y[1:10000]
      y_test <- y[10001:15000]
      
      gp_model <- fitGPModel(gp_coords = coords_train, cov_function = "space_time_gneiting", 
                             vecchia_ordering = "time_random_space",
                             y = y_train,likelihood = "gaussian",gp_approx = "vecchia",num_neighbors = 30,
                             params = list(trace = F,estimate_cov_par_index=c(1,1,1,1,1,1,1,1),
                                           optimizer_cov = "lbfgs"))
      
      gp_model_pred <- predict(gp_model,gp_coords_pred = coords_test, 
                               y = y_train, predict_var = F)
      
      RMSE <- sqrt(mean((gp_model_pred$mu-y_test)^2))
      vec_RMSE <- rep(0,length(y_test)/(n/n_t_steps))
      for (i in 1:length(vec_RMSE)) {
        start <- (i-1)*n/n_t_steps + 1 
        end <- start + n/n_t_steps -1
        vec_RMSE[i] <- sqrt(mean((gp_model_pred$mu[start:end]-y_test[start:end])^2))
      }
      matrix_vec[ii,] <- c(RMSE,vec_RMSE)
      matrixp_vec[ii,] <- gp_model$get_cov_pars()
      gc()
      gp_model <- fitGPModel(gp_coords = coords_train, cov_function = "space_time_gneiting", 
                             vecchia_ordering = "time_random_space",
                             y = y_train,likelihood = "gaussian",gp_approx = "vecchia_correlation_based",num_neighbors = 30,
                             params = list(estimate_cov_par_index=c(1,1,1,1,1,1,1,1),
                                           optimizer_cov = "lbfgs"))
      
      gp_model_pred <- predict(gp_model,gp_coords_pred = coords_test, 
                               y = y_train, predict_var = F)
      
      RMSE <- sqrt(mean((gp_model_pred$mu-y_test)^2))
      vec_RMSE <- rep(0,length(y_test)/(n/n_t_steps))
      for (i in 1:length(vec_RMSE)) {
        start <- (i-1)*n/n_t_steps + 1 
        end <- start + n/n_t_steps -1
        vec_RMSE[i] <- sqrt(mean((gp_model_pred$mu[start:end]-y_test[start:end])^2))
      }
      matrix_vec_corr[ii,] <- c(RMSE,vec_RMSE)
      matrixp_vec_corr[ii,] <- gp_model$get_cov_pars()
      gc()
      gp_model <- fitGPModel(gp_coords = coords_train, cov_function = "space_time_gneiting", 
                             y = y_train,likelihood = "gaussian",gp_approx = "fitc",num_ind_points = 500,ind_points_selection = "kmeans++",
                             params = list(trace = T,estimate_cov_par_index=c(1,1,1,1,1,1,1,1),
                                           optimizer_cov = "lbfgs"))
      
      gp_model_pred <- predict(gp_model,gp_coords_pred = coords_test, 
                               y = y_train, predict_var = F)
      
      RMSE <- sqrt(mean((gp_model_pred$mu-y_test)^2))
      vec_RMSE_3 <- rep(0,length(y_test)/(n/n_t_steps))
      for (i in 1:length(vec_RMSE)) {
        start <- (i-1)*n/n_t_steps + 1 
        end <- start + n/n_t_steps -1
        vec_RMSE_3[i] <- sqrt(mean((gp_model_pred$mu[start:end]-y_test[start:end])^2))
      }
      matrix_fitc[ii,] <- c(RMSE,vec_RMSE_3)
      matrixp_fitc[ii,] <- gp_model$get_cov_pars()
      gc()
      gp_model <- fitGPModel(gp_coords = coords_train, cov_function = "space_time_gneiting", 
                             y = y_train,likelihood = "gaussian",gp_approx = "fitc",num_ind_points = 520,ind_points_selection = "space_time_kmeans++",
                             params = list(estimate_cov_par_index=c(1,1,1,1,1,1,1,1),
                                           optimizer_cov = "lbfgs"))
      
      gp_model_pred <- predict(gp_model,gp_coords_pred = coords_test, 
                               y = y_train, predict_var = F)
      
      RMSE <- sqrt(mean((gp_model_pred$mu-y_test)^2))
      vec_RMSE_3 <- rep(0,length(y_test)/(n/n_t_steps))
      for (i in 1:length(vec_RMSE)) {
        start <- (i-1)*n/n_t_steps + 1 
        end <- start + n/n_t_steps -1
        vec_RMSE_3[i] <- sqrt(mean((gp_model_pred$mu[start:end]-y_test[start:end])^2))
      }
      matrix_fitc_gird[ii,] <- c(RMSE,vec_RMSE_3)
      matrixp_fitc_grid[ii,] <- gp_model$get_cov_pars()
      gc()
      gp_model <- fitGPModel(gp_coords = coords_train, cov_function = "space_time_gneiting", 
                             y = y_train,likelihood = "gaussian",gp_approx = "vif_correlation_based",
                             num_ind_points = 520,ind_points_selection = "space_time_kmeans++",num_neighbors = 30,
                             params = list(estimate_cov_par_index=c(1,1,1,1,1,1,1,1),
                                           optimizer_cov = "lbfgs"))
      
      gp_model_pred <- predict(gp_model,gp_coords_pred = coords_test, 
                               y = y_train, predict_var = F)
      
      RMSE <- sqrt(mean((gp_model_pred$mu-y_test)^2))
      vec_RMSE_3 <- rep(0,length(y_test)/(n/n_t_steps))
      for (i in 1:length(vec_RMSE)) {
        start <- (i-1)*n/n_t_steps + 1 
        end <- start + n/n_t_steps -1
        vec_RMSE_3[i] <- sqrt(mean((gp_model_pred$mu[start:end]-y_test[start:end])^2))
      }
      matrix_fsva[ii,] <- c(RMSE,vec_RMSE_3)
      matrixp_fsva[ii,] <- gp_model$get_cov_pars()
      
}
