#' Run a score test with mean constraint
#' Calculate test statistic and p-value for score test under mean constraint.
#'
#' @param formula_rhs The right hand side of a formula specifying which covariates to include in the model, must be used with the \code{covariate_data} parameter or replaced by the \code{X} parameter.
#' @param Y An outcome matrix with n rows (for samples) and J columns (for KOs) containing coverage data.
#' @param X Design matrix with n rows (for samples) and p columns (for covariates), should have a leading intercept column of \code{1}s, can be replaced by \code{formula_rhs} and \code{covariate_data}.
#' @param covariate_data A data frame including all covariates given in \code{formula_rhs}, can be replaced by design matrix \code{X}.
#' @param B Optional initial parameter estimate for B matrix.
#' @param constraint_cat Category to constrain coefficients to equal the negative sum of all other categories.
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
#' res <- run_score_mc(Y = Y, X = X, ncores = 2, null_k = 2, null_j = 2)
#' 
#' @export
run_score_mc <- function(formula_rhs = NULL,
                          Y,
                          X = NULL,
                          covariate_data = NULL,
                          B = NULL,
                          constraint_cat = 1,
                          null_k = NULL,
                          null_j = NULL, 
                          tolerance = 1e-10,
                          use_tolerance = TRUE, 
                          maxit = 100,
                          maxit_glm = 100,
                          ncores = NULL) {
  
  # hyperparameters 
  n <- nrow(X)
  p <- ncol(X)
  J <- ncol(Y)
  upd_B_alt_ind <- (1:(p*J))[-c(get_theta_ind(j = constraint_cat, k = 1:p, p = p))]
  upd_B_ind <- (1:(p*J))[-c(get_theta_ind(j = constraint_cat, k = 1:p, p = p),
                                get_theta_ind(j = null_j, k = null_k, p = p))]
  upd_ind <- c(upd_B_ind, (p*J + 1):(p*J + n))
  
  # get optimal values under null hypothesis 
  null_res <- fit_null_bcd_mc(formula_rhs = formula_rhs, Y = Y, X = X, 
                              covariate_data = covariate_data, B = B, 
                               constraint_cat = constraint_cat, null_k = null_k, 
                               null_j = null_j, tolerance = tolerance,
                               use_tolerance = use_tolerance, maxit = maxit,
                              maxit_glm = maxit_glm, ncores = ncores)
  B_mle <- null_res$final_B
  z_mle <- null_res$final_z
  
  # get scores using mles under null 
  scores <- compute_scores(X, Y, B_mle, z_mle)
  
  # get information matrix using mles under null 
  info <- matrix(0, nrow = p*J + n, ncol = p*J + n)
  # info block B by B
  J_B <- matrix(NA, nrow = p*(J - 1), ncol = p*(J - 1))
  j_vec <- rep(1:J, each = p)[upd_B_alt_ind]
  k_vec <- rep(1:p, J)[upd_B_alt_ind]
  for (r in 1:(p*(J - 1))) {
    for (c in 1:r) {
      val <- sum(-X[, k_vec[r]]*X[, k_vec[c]]*exp(X %*% B_mle[, 1] + z_mle))
      if (j_vec[r] == j_vec[c]) {
        val <- val + sum(-X[, k_vec[r]]*X[, k_vec[c]]*
                           exp(X %*% B_mle[, j_vec[r]] + z_mle))
      }
      J_B[r, c] <- val
      J_B[c, r] <- val
    }
  }
  info[upd_B_alt_ind, upd_B_alt_ind] <- -J_B
  # info blocks B by z and z by B
  # for (r in 1:(p*(J - 1) - 1)) {
  for (r in 1:(p*(J - 1))) {
    for (c in 1:n) {
      val <- -X[c, k_vec[r]]*exp(X[c, ] %*% B_mle[, 1] + z_mle[c]) + 
        X[c, k_vec[r]]*exp(X[c, ] %*% B_mle[, j_vec[r]] + z_mle[c])
      #info[upd_ind[r], (p*J + c)] <- val
      info[upd_B_alt_ind[r], (p*J + c)] <- val
      #info[(p*J + c), upd_ind[r]] <- val
      info[(p*J + c), upd_B_alt_ind[r]] <- val
    }
  }
  # info block z by z 
  for (i in 1:n) {
    info[p*J + i, p*J + i] <- sum(exp(X[i, ] %*% B_mle + z_mle[i]))
  }
  
  # calculate score statistic 
  null_ind <- get_theta_ind(null_j, null_k, p)
  # inner <- info[null_ind, upd_ind] %*% solve(info[upd_ind, upd_ind]) %*% info[upd_ind, null_ind]
  inner <- info[null_ind, upd_ind] %*% chol2inv(chol(info[upd_ind, upd_ind])) %*% info[upd_ind, null_ind]
  # test_stat <- scores[null_ind] %*% solve(inner) %*% scores[null_ind]
  test_stat <- scores[null_ind] %*% chol2inv(chol(inner)) %*% scores[null_ind]
  p_val <- 1 - pchisq(test_stat, 1)
  
  return(list(null_estimates = null_res, test_stat = test_stat, p_val = p_val))
}
