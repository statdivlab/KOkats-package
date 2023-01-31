#' Compute score vector
#' Compute the score equations for all parameters for given values of B and z.
#' 
#' @param Y An outcome matrix with n rows (for samples) and J columns (for KOs) containing coverage data.
#' @param X Design matrix with n rows (for samples) and p columns (for covariates).
#' @param B Beta parameter matrix with p rows (for covariates) and J columns (for KOs).
#' @param z z parameter vector of length n.
#' 
#' @return The value of each score equation, evaluated at the values of B and z.
#' 
#' @export
compute_scores <- function(X, Y, B, z) {
  n <- nrow(X)
  p <- ncol(X)
  J <- ncol(Y)
  scores <- rep(NA, p*J + n) 
  z_mat <- matrix(z, ncol = 1) %*% matrix(1, ncol = J, nrow = 1)
  log_means <- X %*% B + z_mat
  for (j in 1:J) {
    for (k in 1:p) {
      tmp <- 0 
      for (i in 1:n) {
        tmp <- tmp + Y[i, j]*X[i, k] - X[i, k]*exp(X[i, ] %*% B[, j] + z[i])
      }
      scores[(j - 1)*p + k] <- tmp
    }
  }
  for (i in 1:n) {
    scores[p*J + i] <- sum(Y[i, ] - exp(log_means[i, ]))
  }
  return(scores)
}

