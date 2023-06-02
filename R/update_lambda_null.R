#' Update lambda values for null hypothesis
#' Update lambda values used for constrained estimation under the null hypothesis.
#' 
#' @param Y An outcome matrix with n rows (for samples) and J columns (for KOs) containing coverage data.
#' @param X Design matrix with n rows (for samples) and p columns (for covariates).
#' @param B Beta parameter matrix with p rows (for covariates) and J columns (for KOs).
#' @param z z parameter vector of length n.
#' @param null_j Category tested in null hypothesis.
#' @param null_k Coefficient index tested in null hypothesis.
#' 
#' @return A vector of length p of lambda values.
#' 
#' @export 
update_lambda_null <- function(Y, X, B, z, null_j, null_k) {
  val <- sum(Y[, null_j] * X[, null_k] - X[, null_k] * exp(X %*% B[, null_j] + z))
  return(val)
}
