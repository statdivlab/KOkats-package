#' Estimate under the null model with first category constraint
#' Estimate B and z parameters through block coordinate descent to maximize the unconstrained log likelihood function under a simple null hypothesis and first category constraint.
#'
#' @param formula_rhs The right hand side of a formula specifying which covariates to include in the model, must be used with the \code{covariate_data} parameter or replaced by the \code{X} parameter.
#' @param Y An outcome matrix with n rows (for samples) and J columns (for KOs) containing coverage data.
#' @param X Design matrix with n rows (for samples) and p columns (for covariates), should have a leading intercept column of \code{1}s, can be replaced by \code{formula_rhs} and \code{covariate_data}.
#' @param covariate_data A data frame including all covariates given in \code{formula_rhs}, can be replaced by design matrix \code{X}.
#' @param z Optional initial parameter estimate for z vector.
#' @param null_k Coefficient of the covariate in design matrix that is set to \code{0} under the null hypothesis.
#' @param null_j Category for which covariate \code{k} is set to \code{0} under the null hypothesis.
#' @param tolerance The tolerance used to stop the algorithm when log likelihood values are within \code{tolerance} of each other.
#' @param tolerance_nr The tolerance used to stop the Newton-Raphon method.
#' @param use_tolerance If \code{FALSE}, will run until \code{maxit} regardless of convergence.
#' @param maxit The maximum number of iterations of the coordinate descent algorithm.
#' @param maxit_nr The maximum number of iterations of the Newton-Raphson algorithm within the coordinate descent.
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
#' # don't run this function is deprecated 
#' #res <- fit_bcd_null_fc_constraint(Y = Y, X = X, ncores = 2, null_k = 2, null_j = 2)
#' 
#' @export
fit_bcd_null_fc_constraint <- function(formula_rhs = NULL,
                                  Y,
                                  X = NULL,
                                  covariate_data = NULL,
                                  z = NULL,
                                  null_k,
                                  null_j, 
                                  tolerance = 1e-1,
                                  tolerance_nr = 1e-1,
                                  use_tolerance = TRUE, 
                                  maxit = 100,
                                  maxit_nr = 100,
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
  if (null_j == 1) {
    stop("You cannot test a null hypothesis about category 1 with a category 1 constraint.")
  }

  # set z values to 0 to start if not pre-specified 
  if (is.null(z)) {
    z <- rep(0, n)
  } 
  
  # initialize B matrix
  B <- matrix(0, nrow = p, ncol = J)
  
  # # set B values as solutions to Poisson regression with z as offsets 
  # # update each B vector in parallel using Poisson regression 
  # base_theta <- list(Y = Y, X = X, z = z, maxit_glm = maxit_glm)
  # thetas <- list()
  # for (j in 1:J) {
  #   thetas[[j]] <- base_theta
  #   thetas[[j]]$j <- j
  # }
  # 
  # if (is.null(ncores)) {
  #   cores <- parallel::detectCores() - 1
  # } else {
  #   cores <- ncores
  # }
  # 
  # B_res <- parallel::mclapply(X = thetas, FUN = update_Bj_list, mc.cores = cores)
  # 
  # for (j in 1:J) {
  #   B[, j] <- B_res[[j]]
  # }
  
  # initialize lambda vector 
  # lambda <- update_lambda(Y, X, B, z, j = 1)
  lambda <- rep(1, p)
  lambda_null <- 1
  # lambda_null <- update_lambda_null(Y, X, B, z, null_j, null_k)
  
  # get log likelihood with constraint for initial parameter values 
  constraint_fn <- function(x) {x[, 1]}
  f0 <- compute_loglik_constrained(Y, X, B, z, lambda, lambda_null,
                                   constraint_fn, null_j, null_k)
  
  # update B and z parameters through iteration 
  t <- 1 
  f_old <- -Inf
  f_new <- f0
  
  # save log likelihood values and beta values 
  lik_vec <- rep(NA, maxit)
  B_array <- array(NA, dim = c(nrow(B), ncol(B), maxit))
  z_array <- array(NA, dim = c(length(z), maxit))
  score_mat <- matrix(NA, nrow = (p*J + n + p + 1), ncol = maxit)
  
  if (use_tolerance) {
    condition <- t == 1 || (f_new - f_old)/f_old > tolerance & t < maxit
  } else {
    condition <- t < maxit
  }
  while (condition) {
    
    # update first category with Newton-Raphson 
    B_old <- B
    B_new <- B
    t_nr <- 1
    while (t_nr == 1 | (sum(((B_new[, 1] - B_old[, 1])/B_old[, 1])^2) > tolerance_nr 
                        & t_nr < maxit_nr)) {
      F_B <- update_lambda(Y, X, B_new, z, j = 1) - lambda
      J_B <- compute_J(X, B_new, z, j = 1)
      B_old <- B_new
      B_new[, 1] <- B_old[, 1] - solve(J_B) %*% F_B
      t_nr <- t_nr + 1
    }
    B <- B_new
    
    # update null category with Newton-Raphson
    B_old <- B
    B_new <- B
    t_nr <- 1
    while (t_nr == 1 | (sum(((B_new[, null_j] - B_old[, null_j])/
                             B_old[, null_j])^2) > tolerance_nr 
                        & t_nr < maxit_nr)) {
      e <- rep(0, p)
      e[null_j] <- 1
      F_B <- update_lambda(Y, X, B_new, z, j = null_j) - lambda_null * e
      J_B <- compute_J(X, B_new, z, j = null_j)
      B_old <- B_new
      B_new[, null_j] <- B_old[, null_j] - solve(J_B) %*% F_B
      t_nr <- t_nr + 1
    }
    B <- B_new
    
    # update all other B vectors in parallel using Poisson regression 
    base_theta <- list(Y = Y, X = X, B = B, z = z, maxit_glm = maxit_glm)
    thetas <- list()
    j_set <- (1:J)[-c(1, null_j)]
    for (j in 1:length(j_set)) {
      thetas[[j]] <- base_theta
      thetas[[j]]$j <- j_set[j]
    }
    
    if (is.null(ncores)) {
      cores <- parallel::detectCores() - 1
    } else {
      cores <- ncores
    }
    
    B_res <- parallel::mclapply(X = thetas, FUN = update_Bj_list, mc.cores = cores)
    
    for (j in 1:length(j_set)) {
      B[, j_set[j]] <- B_res[[j]]
    }
    
    # update z values 
    z <- update_z(Y, X, B)
    
    # update lambda vector 
    lambda <- update_lambda(Y, X, B, z, j = 1)
    lambda_null <- update_lambda_null(Y, X, B, z, null_j, null_k)
    
    # update t and likelihood value
    lik_vec[t] <- compute_loglik_constrained(Y, X, B, z, lambda, lambda_null,
                                             constraint_fn, null_j, null_k)
    B_array[, , t] <- B
    z_array[, t] <- z
    f_old <- f_new
    f_new <- lik_vec[t]
    score <- compute_scores_constraint(X, Y, B, z, constraint_fn, null_j, null_k)
    score_mat[, t] <- score
    print(t)
    print(score[55:63, t])
    t <- t + 1
    
    if (use_tolerance) {
      #condition <- (f_new - f_old) > tolerance & t < maxit
      condition <- sum((score)^2) > tolerance & t < maxit
    } else {
      condition <- t < maxit
    }
  }
  
  final_B <- B_array[, , (t-1)]
  final_z <- update_z(Y, X, final_B)
  final_lambda <- update_lambda(Y, X, final_B, final_z, j = 1)
  
  # get rid of NA values if algorithm finished before maxit
  if (t - 1 < maxit) {
    lik_vec <- lik_vec[1:(t - 1)]
    B_array <- B_array[, , 1:(t - 1)]
    z_array <- z_array[, 1:(t - 1)]
  }
  return(list(likelihood = lik_vec, B = B_array, z = z_array,
              final_B = final_B, final_z = final_z, final_lambda = final_lambda,
              score_mat = score_mat))
}
