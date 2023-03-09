#' Generate A prime
#' Create expanded design matrix A prime.
#' 
#' @param X_prime Expanded matrix of X values, nJ rows by pJ columns 
#' @param J The number of KOs. 
#' 
#' @return A_prime, a matrix with \code{nJ} rows and \code{pJ + n} columns.
#' 
#' @export 
generate_A_prime <- function(X_prime, J) {
  n <- nrow(X_prime)/J
  i <- 1:(n*J) 
  j <- rep(1:n, each = J)
  add_cols <- Matrix::sparseMatrix(i = i, j = j, x = 1, dims = c(n*J, n))
  return(cbind(X_prime, add_cols))
}
