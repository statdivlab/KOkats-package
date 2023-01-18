#' Generate X tilde
#' Create expanded design matrix X tilde such that log expected value of Y tilde is equal to X tilde times theta.
#' 
#' @param X Matrix of X values, n rows by p columns 
#' @param J The number of KOs. 
#' 
#' @return X_tilde, a matrix with \code{nJ} rows and \code{pJ + n} columns.
#' 
#' @export 
generate_X_tilde <- function(X, J) {
  n <- nrow(X)
  p <- ncol(X)
  Ei <- diag(nrow = n)
  Ej <- diag(nrow = J) 
  X_tilde <- matrix(0, nrow = n*J, ncol = (p*J + n))
  for (i in 1:n) {
    for (j in 1:J) {
      k <- (i - 1)*J + j 
      X_tilde[k, ] <- c(kronecker(X[i, ], Ej[, j]), Ei[i, ])
    }
  }
  return(X_tilde)
}
