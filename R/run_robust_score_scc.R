#' Run a robust score test with single category constraint
#' Calculate test statistic and p-value for score test under single category constraint.
#'
#' @param formula_rhs The right hand side of a formula specifying which covariates to include in the model, must be used with the \code{covariate_data} parameter or replaced by the \code{X} parameter.
#' @param Y An outcome matrix with n rows (for samples) and J columns (for KOs) containing coverage data.
#' @param X Design matrix with n rows (for samples) and p columns (for covariates), should have a leading intercept column of \code{1}s, can be replaced by \code{formula_rhs} and \code{covariate_data}.
#' @param covariate_data A data frame including all covariates given in \code{formula_rhs}, can be replaced by design matrix \code{X}.
#' @param B Optional initial parameter estimate for B matrix.
#' @param constraint_cat Category to constrain coefficients to equal \code{0}.
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
#' res <- run_robust_score_scc(Y = Y, X = X, ncores = 2, null_k = 2, null_j = 2)
#' 
#' @export
run_robust_score_scc <- function(formula_rhs = NULL,
                          Y,
                          X = NULL,
                          covariate_data = NULL,
                          B = NULL,
                          constraint_cat = 1,
                          null_k = NULL,
                          null_j = NULL, 
                          tolerance = 1e-10,
                          tolerance_nr = 1e-10,
                          use_tolerance = TRUE, 
                          maxit = 100,
                          maxit_nr = 100,
                          maxit_glm = 100,
                          ncores = NULL) {
  
  # hyperparameters 
  n <- nrow(X)
  p <- ncol(X)
  J <- ncol(Y)
  upd_ind <- (1:(p*J + n))[-c(get_theta_ind(j = constraint_cat, k = 1:p, p = p),
                              get_theta_ind(j = null_j, k = null_k, p = p))]
  
  # get optimal values under null hypothesis 
  null_res <- fit_null_bcd_scc_alt(formula_rhs = formula_rhs, Y = Y, X = X, 
                                   covariate_data = covariate_data, B = B, 
                                   constraint_cat = constraint_cat, null_k = null_k,
                                   null_j = null_j, tolerance = tolerance, 
                                   use_tolerance = use_tolerance, maxit = maxit,
                                   maxit_nr = maxit_nr, ncores = ncores)
  B_mle <- null_res$final_B
  z_mle <- null_res$final_z
  
  # get scores using mles under null 
  scores <- compute_scores(X, Y, B_mle, z_mle)
  
  # get information matrix using mles under null 
  B_prime <- as.vector(B_mle)
  theta <- generate_theta(B_prime, z_mle)
  X_prime <- generate_X_prime(X, J)
  A_prime <- generate_A_prime(X_prime, J)
  W <- generate_W(A_prime, theta)
  W_half <- sqrt(W)
  info_left <- Matrix::crossprod(A_prime, W_half) 
  info_right <- W_half %*% A_prime
  info <- info_left %*% info_right
  
  # get covariance of score using mles under null
  js <- (1:J)[-constraint_cat]
  j_vec <- c(rep(js, each = p), rep(0, n))
  k_vec <- c(rep(1:p, J-1), rep(0, n))
  D <- matrix(0, nrow = p*(J - 1) + n, ncol = p*(J - 1) + n)
  for (i in 1:n) {
    for (j in js) {
      score_first <- rep(0, p*(J - 1) + n)
      for (ind in 1:(p*(J - 1))) {
        if (j_vec[ind] == j) {
          score_first[ind] <- Y[i, j] * X[i, k_vec[ind]] - 
            X[i, k_vec[ind]] * exp(X[i, ] %*% B_mle[, j] + z_mle[i])
        }
      }
      for (ind in 1:n) {
        if (ind == i) {
          score_first[ind + p*(J - 1)] <- 1/(J - 1) * (Y[i, 1] - exp(z_mle[i])) + 
            Y[i, j] - exp(X[i, ] %*% B_mle[, j] + z_mle[i])
        }
      }
      D <- D + score_first %*% t(score_first)
    }
  }
  full_D <- matrix(0, nrow = p*J + n, ncol = p*J + n)
  constraint_ind <- get_theta_ind(constraint_cat, 1:p, p)
  full_D[-constraint_ind, -constraint_ind] <- D
  
  # calculate score statistic 
  null_ind <- get_theta_ind(null_j, null_k, p)
  # get variance of score 
  D_11 <- full_D[null_ind, null_ind]
  I_12 <- info[null_ind, upd_ind]
  I_22_inv <- chol2inv(chol(info[upd_ind, upd_ind]))
  D_12 <- full_D[null_ind, upd_ind]
  D_22 <- full_D[upd_ind, upd_ind]
  score_var <- D_11 - I_12 %*% I_22_inv %*% D_12 - 
    D_12 %*% I_22_inv %*% I_12 + 
    I_12 %*% I_22_inv %*% D_22 %*% I_22_inv %*% I_12
  # score stat 
  test_stat <- scores[null_ind] %*% chol2inv(chol(score_var)) %*% scores[null_ind]
  p_val <- 1 - pchisq(test_stat, 1)
  
  return(list(null_estimates = null_res, test_stat = test_stat, p_val = p_val,
              D = full_D))
}
