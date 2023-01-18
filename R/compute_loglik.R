#' Compute log likelihood 
#' Compute the log likelihood for given values of B and z.
#' 
#' @param Y An outcome matrix with n rows (for samples) and J columns (for KOs) containing coverage data.
#' @param X Design matrix with n rows (for samples) and p columns (for covariates).
#' @param B Beta parameter matrix with p rows (for covariates) and J columns (for KOs).
#' @param z z parameter vector of length n.
#' 
#' @return The value of the log likelihood.
#' 
#' @export
compute_loglik <- function(Y, X, B, z) {
  J <- ncol(Y)
  log_means <- X %*% B + matrix(z, ncol = 1) %*% matrix(1, ncol = J, nrow = 1)
  ll <- sum((Y*log_means - exp(log_means)))
  return(ll)
}