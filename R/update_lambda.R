#' Update lambda values for constraint
#' Update lambda values used for constrained estimation under the null hypothesis.
#' 
#' @param Y An outcome matrix with n rows (for samples) and J columns (for KOs) containing coverage data.
#' @param X Design matrix with n rows (for samples) and p columns (for covariates).
#' @param B Beta parameter matrix with p rows (for covariates) and J columns (for KOs).
#' @param z z parameter vector of length n.
#' @param j Category used for constraint. 
#' 
#' @return A vector of length p of lambda values.
#' 
#' @export 
update_lambda <- function(Y, X, B, z, j) {
  p <- ncol(X)
  vec <- rep(0, p)
  for (i in 1:nrow(X)) {
    inner <- Y[i, j] - exp(X[i, ] %*% B[, j] + z[i])
    vec <- vec + X[i, ] * inner 
  }
  return(vec)
}
