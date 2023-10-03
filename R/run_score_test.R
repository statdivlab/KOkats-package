#' Run a robust score test 
#' Calculate test statistic and p-value for robust score test under given constraint.
#'
#' @param formula_rhs The right hand side of a formula specifying which covariates to include in the model, must be used with the \code{covariate_data} parameter or replaced by the \code{X} parameter.
#' @param covariate_data A data frame including all covariates given in \code{formula_rhs}, can be replaced by design matrix \code{X}.
#' @param Y An outcome matrix with n rows (for samples) and J columns (for categories) containing abundance data.
#' @param X Design matrix with n rows (for samples) and p columns (for covariates).
#' @param B Optional initial parameter estimate for B matrix.
#' @param constraint Should be "scc" for single category constraint, "mc" for mean constraint, or "msc" for mean 
#' over a subset constraint.
#' @param constraint_cat Category to constrain coefficients to equal the negative sum of all other categories.
#' @param subset_j Indices of categories to include in constraint for the mean over a subset constraint.
#' @param null_k Coefficient of the covariate in design matrix that is set to \code{0} under the null hypothesis.
#' @param null_j Category for which covariate \code{k} is set to \code{0} under the null hypothesis.
#' @param robust If TRUE runs a robust score test, if FALSE runs a model-based score test.
#' @param tolerance The tolerance used to stop the algorithm when log likelihood values are within \code{tolerance} of each other.
#' @param tolerance_nr The tolerance used to stop the Newton-Raphon method.
#' @param use_tolerance If \code{FALSE}, will run until \code{maxit} regardless of convergence.
#' @param maxit The maximum number of iterations of the coordinate descent algorithm.
#' @param maxit_nr The maximum number of iterations of the Newton-Raphson algorithm within the coordinate descent.
#' @param ncores The desired number of cores to optimize block of B parameters in parallel. If not provided, an appropriate number will be chosen for your machine.
#' @param solve_option way to solve in newton raphson 
#' @param opt_option "bcd" for block coordinate descent or "nr" for newton-raphson
#'
#' @return A list including values of the log likelihood, the B matrix, and the z vector at each iteration.
#'
#' @examples
#' X <- cbind(1, rep(c(0, 1), each = 20))
#' z <- rnorm(40) + 8
#' b0 <- rnorm(10)
#' b1 <- rnorm(10)
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
#' res <- run_score_test(Y = Y, X = X, ncores = 2, null_k = 2, null_j = 2,
#'                       constraint = "scc", solve_option = "solve")
#' 
#' @export
run_score_test <- function(formula_rhs = NULL,
                           Y,
                           X = NULL,
                           covariate_data = NULL,
                           B = NULL,
                           constraint, 
                           constraint_cat = 1,
                           subset_j = NULL,
                           null_k,
                           null_j, 
                           robust = TRUE, 
                           tolerance = 1e-10,
                           tolerance_nr = 1e-10,
                           use_tolerance = TRUE, 
                           maxit = 1000,
                           maxit_nr = 1000,
                           ncores = NULL,
                           solve_option,
                           opt_option = "nr") {
  
  # check for valid constraint
  if (!(constraint %in% c("scc", "mc", "msc"))) {
    stop("Please provide a valid constraint, either 'scc', 'mc', or 'msc'.")
  }
  # check for requirements for mean over a subset constraint
  if (constraint == "msc") {
    if (is.null(subset_j)) {
      stop("If using the 'msc' constraint, please provide 'subset_j'.")
    } else {
      if (!(constraint_cat %in% subset_j)) {
        stop("Please choose a constraint category that is part of 'subset_j'.")
      }
    }
  }
  
  # hyperparameters 
  n <- nrow(X)
  p <- ncol(X)
  J <- ncol(Y)
  upd_ind <- (1:(p*J + n))[-c(get_theta_ind(j = constraint_cat, k = 1:p, p = p),
                              get_theta_ind(j = null_j, k = null_k, p = p))]
  
  # get optimal values under null hypothesis 
  if (opt_option == "bcd") {
    null_res <- fit_null_mle(formula_rhs = formula_rhs, Y = Y, X = X, 
                             covariate_data = covariate_data, B = B, constraint = constraint,
                             constraint_cat = constraint_cat, subset_j = subset_j, null_k = null_k,
                             null_j = null_j, tolerance = tolerance, 
                             tolerance_nr = tolerance_nr, 
                             use_tolerance = use_tolerance, maxit = maxit,
                             maxit_nr = maxit_nr, ncores = ncores,
                             solve_option = solve_option)
  } else {
    null_res <- fit_null_mle_nr(formula_rhs = formula_rhs, Y = Y, X = X, 
                             covariate_data = covariate_data, B = B, constraint = constraint,
                             constraint_cat = constraint_cat, subset_j = subset_j, null_k = null_k,
                             null_j = null_j, tolerance = tolerance, 
                             tolerance_nr = tolerance_nr, 
                             use_tolerance = use_tolerance, maxit = maxit,
                             maxit_nr = maxit_nr, ncores = ncores,
                             solve_option = solve_option)
  }
  B_mle <- null_res$final_B
  z_mle <- null_res$final_z
  
  # get scores using mles under null 
  scores <- compute_scores_cstr(X = X, Y = Y, B = B_mle, z = z_mle, constraint = constraint, 
                                constraint_cat = constraint_cat, subset_j = subset_j)
  
  # get information matrix using mles under null 
  info <- compute_info_cstr(X = X, B = B_mle, z = z_mle, constraint = constraint,
                            constraint_cat = constraint_cat, subset_j = subset_j)
  
  # get covariance of score using mles under null
  null_ind <- get_theta_ind(null_j, null_k, p)
  if (robust) {
    D <- compute_score_var_cstr(Y = Y, X = X, B = B_mle, z = z_mle, constraint = constraint,
                                constraint_cat = constraint_cat, subset_j = subset_j)
    score_var <- compute_var_single_score(null_ind = null_ind, upd_ind = upd_ind, 
                                          p = p, D = D, info = info, robust = TRUE)
  } else {
    score_var <- compute_var_single_score(null_ind = null_ind, upd_ind = upd_ind,
                                          p = p, D = NULL, info = info, robust = FALSE)
  }

  # calculate score statistic 
  test_stat <- scores[null_ind] %*% chol2inv(chol(score_var)) %*% scores[null_ind]
  p_val <- 1 - stats::pchisq(test_stat, 1)
  
  return(list(null_estimates = null_res, test_stat = test_stat, p_val = p_val,
              scores = scores, info = info, D = D))
}
