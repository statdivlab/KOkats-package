#' Compute the Jacobian for Newton-Raphson for updating B vectors.
#' 
#' @param X Design matrix with n rows (for samples) and p columns (for covariates).
#' @param B Beta parameter matrix with p rows (for covariates) and J columns (for KOs).
#' @param z z parameter vector of length n.
#' @param j Category being updated.
#' 
#' @return \code{k} by \code{k} matrix J. 
#' 
#' @export
compute_J <- function(X, B, z, j) {
  p <- ncol(X)
  J <- matrix(NA, nrow = p, ncol = p) 
  for (r in 1:p) {
    for (c in 1:r) {
      val <- sum(X[, r] * X[, c] * exp(X %*% B[, j] + z))
      J[r, c] <- val
      J[c, r] <- val
    }
  }
  return(J)
}
