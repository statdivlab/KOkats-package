#' Generate X prime
#' Create expanded design matrix X prime.
#' 
#' @param X Matrix of X values, n rows by p columns 
#' @param J The number of KOs. 
#' 
#' @return X_prime, a matrix with \code{nJ} rows and \code{pJ} columns.
#' 
#' @export 
generate_X_prime <- function(X, J) {
  n <- nrow(X)
  p <- ncol(X) 
  i <- rep(unlist(lapply(1:J, function(x) {(1:n - 1)*J + x})), each = p)
  j <- unlist(lapply(1:J, function(x) {rep(((x - 1)*p + 1):((x - 1)*p + p), n)}))
  x <- rep(as.vector(t(X)), J)
  X_prime <- Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(n*J, p*J))
  return(X_prime)
}
