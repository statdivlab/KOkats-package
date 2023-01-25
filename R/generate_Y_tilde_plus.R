#' Generate Y tilde plus 
#' Generate an augmented outcome vector to use with Firth penalized likelihood.
#' 
#' @param Y_tilde An outcome vector of length \code{nJ}.
#' @param X_tilde An expanded version of the design matrix.
#' @param W A diagonal matrix with expected values of Y tilde on the diagonal.
#' @param X_tilde_trans The transpose of the \code{X_tilde} matrix.  
#' @param cores The number of cores to use in parallel.  
#'  
#' @return Y_tilde_plus, the augmented outcome vector.
#' 
#' @export
generate_Y_tilde_plus <- function(Y_tilde, X_tilde, W, X_tilde_trans, cores) {
  info <- X_tilde_trans %*% W %*% X_tilde
  info_inv <- MASS::ginv(as.matrix(info))
  info_inv <- Matrix::Matrix(info_inv, sparse = TRUE)
  W_half <- sqrt(W)
  # augmented portion 
  nJ <- nrow(W)
  left <- W_half %*% X_tilde
  right <- X_tilde_trans %*% W_half 
  #aug_vec <- rep(0, nJ)
  #for (i in 1:nJ) {
  #  aug_vec[i] <- left[i, ] %*% info_inv %*% right[, i]
  #}
  #aug_mat <- W_half %*% X_tilde %*% info_inv %*% X_tilde_trans %*% W_half
  #aug <- spam::diag(aug_mat)
  aug_res <- parallel::mclapply(X = 1:nJ, FUN = get_augment_diag, mc.cores = cores)
  aug_vec <- as.vector(unlist(aug_res))
  return(Y_tilde + aug_vec/2)
}

