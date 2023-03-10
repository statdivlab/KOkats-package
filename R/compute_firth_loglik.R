#' Compute Firth penalized log likelihood 
#' Compute the Firth penalized log likelihood for given values of B and z.
#' 
#' @param Y An outcome matrix with n rows (for samples) and J columns (for KOs) containing coverage data.
#' @param X Design matrix with n rows (for samples) and p columns (for covariates).
#' @param B Beta parameter matrix with p rows (for covariates) and J columns (for KOs).
#' @param z z parameter vector of length n.
#' @param X_aug An expanded version of the design matrix.
#' @param W A diagonal matrix with expected values of Y tilde on the diagonal.
#' 
#' @return The value of the log likelihood.
#' 
#' @export
compute_firth_loglik <- function(Y, X, B, z, X_aug, W) {
  # get log likelihood 
  loglik <- compute_loglik(Y, X, B, z)
  # compute information matrix using fact that I = X_tilde^T W X_tilde
  # info <- Matrix::crossprod(X_aug, W) %*% X_aug
  # eigs <- eigen(info)$values
  # non_zero_eigs <- eigs[eigs > 0]
  # log_non_zero_eigs <- log(non_zero_eigs)
  # firth_portion <- 1/2*sum(log_non_zero_eigs)
  
  firth_portion <- 0.5*Matrix::determinant(
    Matrix::crossprod(X_aug, W) %*% X_aug)$modulus
  
  return(loglik + firth_portion)
}
