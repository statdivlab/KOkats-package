#' Generate W
#' Create diagonal matrix with expected value of Y tilde on diagonals. 
#' 
#' @param X_tilde An expanded version of the design matrix X. 
#' @param theta A vector of all B and z parameters 
#' 
#' @return W, a diagonal matrix with \code{nJ} rows and columns. 
#' 
#' @export 
generate_W <- function(X_tilde, theta) {
  W <- Matrix::Diagonal(x = as.vector(exp(X_tilde %*% theta)))
  return(W)
}
