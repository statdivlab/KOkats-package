#' Fit Algorithm 2
#' Estimate B and z parameters through block coordinate descent to maximize the Firth penalized log likelihood function, details given in Algorithm 2.
#'
#' @param formula_rhs The right hand side of a formula specifying which covariates to include in the model, must be used with the \code{covariate_data} parameter or replaced by the \code{X} parameter.
#' @param Y An outcome matrix with n rows (for samples) and J columns (for KOs) containing coverage data.
#' @param X Design matrix with n rows (for samples) and p columns (for covariates), should have a leading intercept column of \code{1}s, can be replaced by \code{formula_rhs} and \code{covariate_data}.
#' @param covariate_data A data frame including all covariates given in \code{formula_rhs}, can be replaced by design matrix \code{X}.
#' @param z Optional initial parameter estimate for z vector.
#' @param tolerance The tolerance used to stop the algorithm when log likelihood values are within \code{tolerance} of each other.
#' @param maxit The maximum number of iterations of the coordinate descent algorithm.
#' @param constraint_fn A constraint function to make the B matrix identifiable.
#' @param maxit_glm The maximum number of iterations when running the glm to update the block of Bj parameters in the coordinate descent algorithm.
#' @param ncores The desired number of cores to optimize block of B parameters in parallel. If not provided, an appropriate number will be chosen for your machine.
#'
#' @return A list including values of the log likelihood, the B matrix, and the z vector at each iteration.
#'
#' @examples
#' n <- 100
#' J <- 100
#' X <- cbind(1, rep(c(0, 1), each = n/2))
#' z <- rnorm(n) + 5
#' b0 <- rnorm(J)
#' b1 <- rnorm(J)
#' b <- rbind(b0, b1)
#' Y <- matrix(NA, ncol = J, nrow = n)
#' 
#' for (i in 1:n) {
#'  for (j in 1:J) {
#'    temp_mean <- exp(X[i, , drop = FALSE] %*% b[, j, drop = FALSE] + z[i])
#'    Y[i,j] <- rpois(1, lambda = temp_mean)
#'  }
#' }
#' 
#' start <- proc.time()
#' res <- fit_penalized_bcd_unconstrained_fast(Y = Y, X = X, 
#'           constraint_fn = function(x) {mean(x)}, ncores = 2)
#' proc.time() - start
#' 
#' @export
fit_penalized_bcd_unconstrained_fast <- function(formula_rhs = NULL,
                                            Y,
                                            X = NULL,
                                            covariate_data = NULL,
                                            z = NULL,
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
  
  # get X_star and Y_tilde 
  X_prime <- generate_X_prime(X, J)
  A_prime <- generate_A_prime(X_prime, J)
  Y_tilde <- generate_Y_tilde(Y)
  
  # set B values to 0 
  B <- matrix(0, nrow = p, ncol = J)
  
  # set z values to 0 if not prespecified
  if (is.null(z)) {
    z <- rep(0, n)
  } 
  
  # transform parameters to get theta and W 
  B_prime <- as.vector(B)
  theta <- generate_theta(B_prime, z)
  W <- generate_W(A_prime, theta)
  
  # get information matrix 
  W_half <- sqrt(W)
  info_left <- Matrix::crossprod(X_prime, W_half) 
  info_right <- W_half %*% X_prime
  info <- info_left %*% info_right
  
  # update B and z parameters through iteration 
  t <- 1 
  B_old <- matrix(-Inf, nrow = p, ncol = J)
  B_new <- B
  
  # save log likelihood values and beta values 
  lik_vec <- rep(NA, maxit)
  score_vec <- rep(NA, maxit)
  B_array <- array(NA, dim = c(nrow(B), ncol(B), maxit))
  z_array <- array(NA, dim = c(length(z), maxit))
  
  while (t == 1 | 
         (sum(((B_new - B_old)/B_old)^2) > tolerance & t < maxit)) {
    
    # figure out number of cores for parallel actions 
    if (is.null(ncores)) {
      cores <- parallel::detectCores() - 1
    } else {
      cores <- ncores
    }
    
    # compute Y+
    info_inv_blocks <- parallel::mclapply(X = 1:J, FUN = invert_block,
                                          p = p, info = info,
                                          mc.cores = cores)
    info_inv <- Matrix::bdiag(info_inv_blocks)
    nJ <- nrow(W)
    
    aug_res <- info_right %*% info_inv %*% info_left
    aug_res_T <- as(aug_res, "TsparseMatrix")
    inds <- which(aug_res_T@i == aug_res_T@j)
    aug_vec <- aug_res@x[inds]
    
    Y_tilde_plus <- Y_tilde + aug_vec/2
    Y_plus <- generate_Y_plus(Y_tilde_plus, J)
    
    # update each B vector in parallel using Poisson regression 
    base_theta <- list(Y = Y_plus, X = X, B = B, z = z, maxit_glm = maxit_glm)
    thetas <- list()
    for (j in 1:J) {
      thetas[[j]] <- base_theta
      thetas[[j]]$j <- j
    }
    
    B_res <- parallel::mclapply(X = thetas, FUN = update_Bj_list, mc.cores = cores)
    
    for (j in 1:J) {
      B[, j] <- B_res[[j]]
    }
    
    # update z values 
    z <- update_z(Y_plus, X, B)
    
    # update t and likelihood value
    B_prime <- as.vector(B)
    theta <- generate_theta(B_prime, z)
    W <- generate_W(A_prime, theta)
    W_half <- sqrt(W)
    info_left <- Matrix::crossprod(X_prime, W_half) 
    info_right <- W_half %*% X_prime
    info <- info_left %*% info_right
    lik_vec[t] <- compute_loglik(Y, X, B, z) + 
      0.5*Matrix::determinant(info)$modulus
    score_vec[t] <- sum(compute_scores(X, Y_plus, B, z)^2)
    B_array[, , t] <- B
    z_array[, t] <- z
    B_old <- B_new
    B_new <- B
    t <- t + 1
  }
  
  final_B <- B_array[, , (t-1)]
  # enforce identifiability constraint 
  for (k in 1:p) {
    final_B[k, ] <- final_B[k, ] - constraint_fn(final_B[k, ])
  }
  final_z <- update_z(Y, X, final_B)
  
  # get rid of NA values if algorithm finished before maxit
  if (t - 1 < maxit) {
    lik_vec <- lik_vec[1:(t - 1)]
    score_vec <- score_vec[1:(t - 1)]
    B_array <- B_array[, , 1:(t - 1)]
    z_array <- z_array[, 1:(t - 1)]
  }
  return(list(likelihood = lik_vec, score = score_vec, 
              B = B_array, z = z_array,
              final_B = final_B, final_z = final_z))
}
