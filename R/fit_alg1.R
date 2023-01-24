#' Fit Algorithm 1
#' Estimate B and z parameters through block coordinate descent to maximize the log likelihood function, details given in Algorithm 1.
#'
#' @param formula_rhs The right hand side of a formula specifying which covariates to include in the model, must be used with the \code{covariate_data} parameter or replaced by the \code{X} parameter.
#' @param Y An outcome matrix with n rows (for samples) and J columns (for KOs) containing coverage data.
#' @param X Design matrix with n rows (for samples) and p columns (for covariates), should have a leading intercept column of \code{1}s, can be replaced by \code{formula_rhs} and \code{covariate_data}.
#' @param covariate_data A data frame including all covariates given in \code{formula_rhs}, can be replaced by design matrix \code{X}.
#' @param B Optional initial parameter estimate for B matrix 
#' @param tolerance The tolerance used to stop the algorithm when log likelihood values are within \code{tolerance} of each other.
#' @param maxit The maximum number of iterations of the coordinate descent algorithm.
#' @param constraint_fn A constraint function to make the B matrix identifiable.
#' @param maxit_glm The maximum number of iterations when running the glm to update the block of Bj parameters in the coordinate descent algorithm.
#' @param ncores The desired number of cores to optimize block of B parameters in parallel. If not provided, an appropriate number will be chosen for your machine.
#'
#' @return A list including values of the log likelihood, the B matrix, and the z vector at each iteration.
#'
#' @examples
#' X <- cbind(1, rep(c(0, 1), each = 20))
#' z <- rnorm(40) + 8
#' b0 <- rnorm(10)
#' b1 <- 1:10
#' b <- rbind(b0, b1)
#' Y <- matrix(NA, ncol = 10, nrow = 40)
#' 
#' for (i in 1:40) {
#'  for (j in 1:10) {
#'    temp_mean <- exp(X[i, , drop = FALSE] %*% b[, j, drop = FALSE] + z[i])
#'    Y[i,j] <- rpois(1, lambda = temp_mean)
#'  }
#' }
#' 
#' res <- fit_alg1(Y = Y, X = X, constraint_fn = function(x) {median(x)}, ncores = 2)
#' 
#' @export
fit_alg1 <- function(formula_rhs = NULL,
                      Y,
                      X = NULL,
                      covariate_data = NULL,
                      B = NULL,
                      tolerance = 1e-1,
                      maxit = 100,
                      constraint_fn,
                      maxit_glm = 100,
                      ncores = NULL) {
  
  # check that data has been given with either X matrix or formula and covariate data
  if(!is.null(formula_rhs)){
    if(is.null(covariate_data)){
      stop("If formula_rhs is provided, covariates named in formula must be
           provided inside covariate_data.")
    }
    if(!is.data.frame(covariate_data)){
      stop("Argument covariate_data must be a data frame.")
    }
    if(!inherits(formula_rhs,"formula")){
      formula_rhs <- stats::as.formula(formula_rhs)
    }
    X <- stats::model.matrix(formula_rhs,covariate_data)
  }
  
  # set up initial hyperparameters and parameters 
  p <- ncol(X)
  if (p < 2) {
    stop("X must contain an intercept column and at least one other (linearly
          independent) column")
  }
  J <- ncol(Y)
  n <- nrow(X)
  
  # set B values to 0 to start if not pre-specified
  if (is.null(B)) {
    B <- matrix(0, nrow = p, ncol = J)
  }
  
  # find j*, KO to use for initial constraint 
  j_star <- which.max(colSums(Y > 0))
  
  # set temporary identifiability constraint as B^(j*(0)) = 0_p 
  for(k in 1:p){
    B[k, ] <- B[k, ] - B[k, j_star]
  }
  
  # set z_i^(0) for each sample i 
  Y_cross <- pmax(1, Y[, j_star])
  # is this the right B*? If so, won't it always be exp(0) = 1? 
  z0 <- log(Y_cross) - X %*% B[, j_star]
  
  # get log likelihood for initial parameter values 
  f0 <- compute_loglik(Y, X, B, z0)
  
  # update B and z parameters through iteration 
  z <- z0
  t <- 1 
  f_old <- -Inf
  f_new <- f0
  
  # save log likelihood values and beta values 
  lik_vec <- rep(NA, maxit)
  B_array <- array(NA, dim = c(nrow(B), ncol(B), maxit))
  z_array <- array(NA, dim = c(length(z), maxit))
  
  while ((f_new - f_old) > tolerance & t < maxit) {
    
    # update each B vector in parallel using Poisson regression 
    base_theta <- list(Y = Y, X = X, B = B, z = z, maxit_glm = maxit_glm)
    thetas <- list()
    for (j in 1:J) {
      thetas[[j]] <- base_theta
      thetas[[j]]$j <- j
    }
    
    if (is.null(ncores)) {
      cores <- parallel::detectCores() - 1
    } else {
      cores <- ncores
    }
    
    B_res <- parallel::mclapply(X = thetas, FUN = update_Bj_list, mc.cores = cores)
    
    for (j in 1:J) {
      B[, j] <- B_res[[j]]
    }
    print(B)
    
    # enforce identifiability constraint 
    for (k in 1:p) {
      B[k, ] <- B[k, ] - constraint_fn(B[k, ])
    }
    print(B)
    
    # update z values 
    z <- update_z(Y, X, B)
    print(z)
    
    # update t and likelihood value
    lik_vec[t] <- compute_loglik(Y, X, B, z)
    print(lik_vec[t])
    B_array[, , t] <- B
    z_array[, t] <- z
    f_old <- f_new
    f_new <- lik_vec[t]
    t <- t + 1
  }
  
  # get rid of NA values if algorithm finished before maxit
  if (t - 1 < maxit) {
    lik_vec <- lik_vec[1:(t - 1)]
    B_array <- B_array[, , 1:(t - 1)]
    z_array <- z_array[, 1:(t - 1)]
  }
  return(list(likelihood = lik_vec, B = B_array, z = z_array))
}
