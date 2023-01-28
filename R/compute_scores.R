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

