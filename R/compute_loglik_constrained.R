#' Compute log likelihood with constraints
#' Compute the log likelihood for given values of B and z and constraints.
#' 
#' @param Y An outcome matrix with n rows (for samples) and J columns (for KOs) containing coverage data.
#' @param X Design matrix with n rows (for samples) and p columns (for covariates).
#' @param B Beta parameter matrix with p rows (for covariates) and J columns (for KOs).
#' @param z z parameter vector of length n.
#' @param lambda Vector of Lagrange multipliers for identifiability constraints. 
#' @param lambda_null Lagrange multiplier for null hypothesis constraint.
#' @param constraint_fn A constraint function to make the B matrix identifiable.#' 
#' @param null_j Category for which covariate \code{k} is set to \code{0} under the null hypothesis.
#' @param null_k Coefficient of the covariate in design matrix that is set to \code{0} under the null hypothesis.
#' 
#' @return The value of the log likelihood.
#' 
#' @export
compute_loglik_constrained <- function(Y, X, B, z, lambda, lambda_null, constraint_fn,
                           null_j, null_k) {
  J <- ncol(Y)
  log_means <- X %*% B + matrix(z, ncol = 1) %*% matrix(1, ncol = J, nrow = 1)
  ll <- sum((Y*log_means - exp(log_means))) - 
    lambda %*% constraint_fn(B) - lambda_null * B[null_k, null_j]
  return(ll)
}
