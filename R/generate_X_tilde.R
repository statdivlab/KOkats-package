#' Generate X tilde
#' Create expanded design matrix X tilde such that log expected value of Y tilde is equal to X tilde times theta.
#' 
#' @param X Matrix of X values, n rows by p columns 
#' @param J The number of KOs. 
#' @param transpose If true, will generate the transpose of X_tilde. Set to \code{FALSE}.
#' 
#' @return X_tilde, a matrix with \code{nJ} rows and \code{pJ + n} columns.
#' 
#' @export 
generate_X_tilde <- function(X, J, transpose = FALSE) {
  n <- nrow(X)
  p <- ncol(X) 
  if (!transpose) {
    i <- rep(unlist(lapply(1:J, function(x) {(1:n - 1)*J + x})), each = p)
    j <- unlist(lapply(1:J, function(x) {rep((1:p - 1)*J + x, n)}))
    x <- rep(as.vector(t(X)), J)
    X_tilde <- Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(n*J, p*J))
    i <- 1:(n*J) 
    j <- rep(1:n, each = J)
    add_cols <- Matrix::sparseMatrix(i = i, j = j, x = 1, dims = c(n*J, n))
    return(cbind(X_tilde, add_cols))
  } else {
    i <- rep(unlist(lapply(1:J, function(x) {(1:n - 1)*J + x})), each = p)
    j <- unlist(lapply(1:J, function(x) {rep((1:p - 1)*J + x, n)}))
    x <- rep(as.vector(t(X)), J)
    X_tilde <- Matrix::sparseMatrix(i = j, j = i, x = x, dims = c(p*J, n*J))
    i <- 1:(n*J) 
    j <- rep(1:n, each = J)
    add_cols <- Matrix::sparseMatrix(i = j, j = i, x = 1, dims = c(n, n*J))
    return(rbind(X_tilde, add_cols))
  }
}

