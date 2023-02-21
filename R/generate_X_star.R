#' Generate X star
#' Create expanded design matrix X star (the first pJ columns of X_tilde).
#' 
#' @param X Matrix of X values, n rows by p columns 
#' @param J The number of KOs. 
#' @param transpose If true, will generate the transpose of X_star. Set to \code{FALSE}.
#' 
#' @return X_star, a matrix with \code{nJ} rows and \code{pJ} columns.
#' 
#' @export 
generate_X_star <- function(X, J, transpose = FALSE) {
  n <- nrow(X)
  p <- ncol(X) 
  if (!transpose) {
    i <- rep(unlist(lapply(1:J, function(x) {(1:n - 1)*J + x})), each = p)
    j <- unlist(lapply(1:J, function(x) {rep((1:p - 1)*J + x, n)}))
    x <- rep(as.vector(t(X)), J)
    X_star <- Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(n*J, p*J))
    return(X_star)
  } else {
    i <- rep(unlist(lapply(1:J, function(x) {(1:n - 1)*J + x})), each = p)
    j <- unlist(lapply(1:J, function(x) {rep((1:p - 1)*J + x, n)}))
    x <- rep(as.vector(t(X)), J)
    X_star <- Matrix::sparseMatrix(i = j, j = i, x = x, dims = c(p*J, n*J))
    return(X_star)
  }
}