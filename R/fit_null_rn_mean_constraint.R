#' Estimate under the null model with mean constraint
#' Estimate B and z parameters through Newton-Raphson to maximize the unconstrained log likelihood function under a simple null hypothesis and mean constraint.
#'
#' @param formula_rhs The right hand side of a formula specifying which covariates to include in the model, must be used with the \code{covariate_data} parameter or replaced by the \code{X} parameter.
#' @param Y An outcome matrix with n rows (for samples) and J columns (for KOs) containing coverage data.
#' @param X Design matrix with n rows (for samples) and p columns (for covariates), should have a leading intercept column of \code{1}s, can be replaced by \code{formula_rhs} and \code{covariate_data}.
#' @param covariate_data A data frame including all covariates given in \code{formula_rhs}, can be replaced by design matrix \code{X}.
#' @param B Optional initial parameter estimate for B matrix. 
#' @param z Optional initial parameter estimate for z vector.
#' @param null_k Coefficient of the covariate in design matrix that is set to \code{0} under the null hypothesis.
#' @param null_j Category for which covariate \code{k} is set to \code{0} under the null hypothesis.
#' @param tolerance The tolerance used to stop the algorithm when log likelihood values are within \code{tolerance} of each other.
#' @param use_tolerance If \code{FALSE}, will run until \code{maxit} regardless of convergence.
#' @param maxit The maximum number of iterations of the coordinate descent algorithm.
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
#' res_mean <- fit_null_nr_mean_constraint(Y = Y, X = X, ncores = 2, null_k = 2, null_j = 2)
#' 
#' @export
fit_null_nr_mean_constraint <- function(formula_rhs = NULL,
                                      Y,
                                      X = NULL,
                                      covariate_data = NULL,
                                      B = NULL, 
                                      z = NULL,
                                      null_k = NULL,
                                      null_j = NULL, 
                                      tolerance = 1e-1,
                                      use_tolerance = TRUE, 
                                      maxit = 100,
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
  
  if (is.null(null_j) | is.null(null_k)) {
    h0 <- FALSE
  } else {
    h0 <- TRUE
  }
  
  # set up initial hyperparameters and parameters 
  p <- ncol(X)
  if (p < 2) {
    stop("X must contain an intercept column and at least one other (linearly
          independent) column")
  }
  J <- ncol(Y)
  n <- nrow(X)
  
  # find optimal parameters without null hypothesis constraint 
  if (is.null(B)) {
    initial_constr <- function(x) {mean(x)}
    res <- fit_bcd_unconstrained(Y = Y, X = X, tolerance = tolerance, maxit = maxit,
                                 maxit_glm = maxit_glm, ncores = ncores, 
                                 constraint_fn = initial_constr)
    z <- res$final_z
    B <- res$final_B
  } else {
    z <- update_z(Y, X, B)
  }
  # set null hypothesis constraint 
  if (h0) {
    B[null_k, null_j] <- 0
  }
  
  # get log likelihood with constraint for initial parameter values 
  f0 <- compute_loglik(Y, X, B, z)
  
  # update B and z parameters through iteration 
  t <- 1 
  f_old <- -Inf
  f_new <- f0
  
  # save log likelihood values and beta values 
  lik_vec <- rep(NA, maxit)
  B_array <- array(NA, dim = c(nrow(B), ncol(B), maxit))
  z_array <- array(NA, dim = c(length(z), maxit))
  score_mat <- matrix(NA, nrow = (p*J + n), ncol = maxit)
  
  if (use_tolerance) {
    condition <- t == 1 || (f_new - f_old)/f_old > tolerance & t < maxit
  } else {
    condition <- t < maxit
  }
  while (condition) {
    
    # get score equations 
    score_vec <- compute_scores_mean_constr(X, Y, B, z)
    
    # augmented design matrix 
    X_prime <- generate_X_prime(X, J)
    A_prime <- generate_A_prime(X_prime, J)
    
    # transform parameters to get theta and W 
    B_prime <- as.vector(B)
    theta <- generate_theta(B_prime, z)
    
    # add extra components on info matrix 
    info_mean <- matrix(0, nrow = (p*J + n), ncol = (p*J + n))
    i_vec <- c(rep(0, p*J), 1:n)
    j_vec <- c(rep(1:J, each = p), rep(0, n))
    k_vec <- c(rep(1:p, J), rep(0, n))
    for (r in (p + 1):(p*J + n)) {
      for (c in 1:r) {
        if (r <= p*(J)) {
          if (j_vec[r] == j_vec[c]) {
            val1 <- sum(X[, k_vec[r]]*X[, k_vec[c]]*
                          exp(X %*% B[, j_vec[r]] + z))
            val2 <- sum(X[, k_vec[r]]*X[, k_vec[c]]*
                          exp(X %*% B + z))
            info_mean[r, c] <- -(val1 + val2)
            info_mean[c, r] <- -(val1 + val2)
          } else {
            val <- sum(X[, k_vec[r]]*X[, k_vec[c]]*
                         exp(X %*% B[, 1] + z))
            info_mean[r, c] <- -val
            info_mean[c, r] <- -val 
          }
        }
        else {
          if (c <= p*(J)) {
            val1 <- X[i_vec[r], k_vec[c]]*
              exp(X[i_vec[r], ] %*% B[, j_vec[c]] + z[i_vec[r]])
            val2 <- -X[i_vec[r], k_vec[c]]*
              exp(X[i_vec[r], ] %*% B[, 1] + z[i_vec[r]])
            info_mean[r, c] <- -(val1 + val2)
            info_mean[c, r] <- -(val1 + val2)
          } else if (r == c) {
            val <- sum(exp(X[i_vec[r], ] %*% B + z[i_vec[r]]))
            info_mean[r, c] <- -val
            info_mean[c, r] <- -val
          }
        }
      }
    }
    
    # just take part of theta vector, score, and info matrix for update 
    if (h0) {
      ind <- (1:(p*J + n))[-c(get_theta_ind(j = 1, k = 1:p, p = p),
                              get_theta_ind(j = null_j, k = null_k, p = p))]
    } else {
      ind <- (1:(p*J + n))[-c(get_theta_ind(j = 1, k = 1:p, p = p))]
    }
    
    theta_sm <- theta[ind]
    score_sm <- score_vec[ind]
    info_sm <- info_mean[ind, ind]
    
    # Newton-Raphson update 
    theta_new <- theta_sm - solve(info_sm) %*% score_sm 
    theta[ind] <- theta_new
    B <- matrix(theta[1:(p*J)], nrow = p, ncol = J)
    B[, 1] <- -rowSums(B[, 2:J])
    z <- theta[(p*J + 1):(p*J + n)]
    
    # save info 
    lik_vec[t] <- compute_loglik(Y, X, B, z)
    B_array[, , t] <- B
    z_array[, t] <- z
    f_old <- f_new
    f_new <- lik_vec[t]
    score <- compute_scores(X, Y, B, z)
    score_mat[, t] <- score
    t <- t + 1
    
    if (use_tolerance) {
      condition <- (f_new - f_old) > tolerance & t < maxit
      #condition <- sum((score)^2) > tolerance & t < maxit
    } else {
      condition <- t < maxit
    }
  }
  
  final_B <- B_array[, , (t-1)]
  final_z <- z_array[, (t-1)]
  
  # get rid of NA values if algorithm finished before maxit
  if (t - 1 < maxit) {
    lik_vec <- lik_vec[1:(t - 1)]
    B_array <- B_array[, , 1:(t - 1)]
    z_array <- z_array[, 1:(t - 1)]
    score_mat <- score_mat[, 1:(t - 1)]
  }
  return(list(likelihood = lik_vec, B = B_array, z = z_array,
              final_B = final_B, final_z = final_z, score_mat = score_mat))
}
