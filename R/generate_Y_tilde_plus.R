#' Generate Y tilde plus 
#' Generate an augmented outcome vector to use with Firth penalized likelihood.
#' 
#' @param Y_tilde An outcome vector of length \code{nJ}.
#' @param X_tilde An expanded version of the design matrix.
#' @param W A diagonal matrix with expected values of Y tilde on the diagonal.
#'  
#' @return Y_tilde_plus, the augmented outcome vector.
#' 
#' @export
generate_Y_tilde_plus <- function(Y_tilde, X_tilde, W, X_tilde_trans) {
  info <- X_tilde_trans %*% W %*% X_tilde
  info_inv <- MASS::ginv(as.matrix(info))
  info_inv <- Matrix::Matrix(info_inv, sparse = TRUE)
  W_half <- sqrt(W)
  # augmented portion 
  nJ <- nrow(W)
  left <- W_half %*% X_tilde
  right <- X_tilde_trans %*% W_half 
  aug_vec <- vector(0, length = nJ)
  for (i in 1:nJ) {
    aug_vec <- left[i, ] %*% info_inv %*% right[, i]
  }
  #aug_mat <- W_half %*% X_tilde %*% info_inv %*% X_tilde_trans %*% W_half
  #aug <- spam::diag(aug_mat)
  return(Y_tilde + aug_vec/2)
}
