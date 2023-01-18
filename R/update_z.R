#' Update z values 
#' Use the closed form solution to the score equation for z to update the block of z parameters in the block coordinate descent algorithm.
#' 
#' @param Y An outcome matrix with n rows (for samples) and J columns (for KOs) containing coverage data.
#' @param X Design matrix with n rows (for samples) and p columns (for covariates).
#' @param B Beta parameter matrix with p rows (for covariates) and J columns (for KOs).
#'
#' @return A vector of new values of z. 
#' 
#' @export 
update_z <- function(Y, X, B) {
  z <- log(Matrix::rowSums(Y)) - log(Matrix::rowSums(exp(X %*% B)))
  return(z)
}