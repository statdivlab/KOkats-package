#' Estimate under the null model with single category constraint
#' Estimate B and z parameters through Newton-Raphson to maximize the unconstrained log likelihood function under a simple null hypothesis and single category constraint.
#'
#' @param formula_rhs The right hand side of a formula specifying which covariates to include in the model, must be used with the \code{covariate_data} parameter or replaced by the \code{X} parameter.
#' @param Y An outcome matrix with n rows (for samples) and J columns (for KOs) containing coverage data.
#' @param X Design matrix with n rows (for samples) and p columns (for covariates), should have a leading intercept column of \code{1}s, can be replaced by \code{formula_rhs} and \code{covariate_data}.
#' @param covariate_data A data frame including all covariates given in \code{formula_rhs}, can be replaced by design matrix \code{X}.
#' @param B Optional initial parameter estimate for B matrix. 
#' @param z Optional initial parameter estimate for z vector.
#' @param constraint_cat Category to constrain coefficients to equal \code{0}.
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
#' res_first <- fit_null_nr_fc_constraint(Y = Y, X = X, ncores = 2, null_k = 2, null_j = 2, tolerance = 1e-6)
#' 
#' @export
fit_null_nr_fc_constraint <- function(formula_rhs = NULL,
                                       Y,
                                       X = NULL,
                                       covariate_data = NULL,
                                       B = NULL,
                                       z = NULL,
                                       constraint_cat = 1, 
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
    if (null_j == constraint_cat) {
      stop("You cannot test a null hypothesis about category 1 with a category 1 constraint.")
    }
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
    initial_constr <- function(x) {x[constraint_cat]}
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
    score_vec <- compute_scores(X, Y, B, z)
    
    # augmented design matrix 
    X_prime <- generate_X_prime(X, J)
    A_prime <- generate_A_prime(X_prime, J)
    
    # transform parameters to get theta and W 
    B_prime <- as.vector(B)
    theta <- generate_theta(B_prime, z)
    W <- generate_W(A_prime, theta)
    
    # get negative information matrix 
    W_half <- sqrt(W)
    info_left <- Matrix::crossprod(A_prime, W_half) 
    info_right <- W_half %*% A_prime
    info <- info_left %*% info_right
    neg_info <- -info
    
    # just take part of theta vector, score, and info matrix for update 
    if (h0) {
      ind <- (1:(p*J + n))[-c(get_theta_ind(j = constraint_cat, k = 1:p, p = p),
                              get_theta_ind(j = null_j, k = null_k, p = p))]
    } else {
      ind <- (1:(p*J + n))[-c(get_theta_ind(j = constraint_cat, k = 1:p, p = p))]
    }
    
    theta_sm <- theta[ind]
    score_sm <- score_vec[ind]
    info_sm <- neg_info[ind, ind]
    
    # Newton-Raphson update 
    theta_new <- theta_sm - solve(info_sm) %*% score_sm 
    theta[ind] <- theta_new
    B <- matrix(theta[1:(p*J)], nrow = p, ncol = J)
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

fc_via_optim <- function(Y, X) {
  par <- rep(0, 57)
  optim(par, fn_to_min, control = list(maxit = 500000))
}

fn_to_min <- function(theta_sm) {
  theta <- c(rep(0, 2), theta_sm[1], 0, theta_sm[2:57])
  B <- matrix(theta[1:20], nrow = 2)
  z <- theta[21:60]
  return(-compute_loglik(Y, X, B, z))
}
