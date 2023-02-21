#' Generate Y tilde plus 
#' Generate an augmented outcome vector to use with Firth penalized likelihood.
#' 
#' @param Y_tilde An outcome vector of length \code{nJ}.
#' @param X_aug An expanded version of the design matrix.
#' @param W A diagonal matrix with expected values of Y tilde on the diagonal.
#' @param cores The number of cores to use in parallel.  
#'  
#' @return Y_tilde_plus, the augmented outcome vector.
#' 
#' @export
generate_Y_tilde_plus <- function(Y_tilde, X_aug, W, cores) {
  W_half <- sqrt(W)
  info_left <- Matrix::crossprod(X_aug, W_half)
  info_right <- W_half %*% X_aug
  info <- info_left %*% info_right
  info_inv <- MASS::ginv(as.matrix(info))
  info_inv <- Matrix::Matrix(info_inv, sparse = TRUE)
  # augmented portion 
  nJ <- nrow(W)
  aug_res <- parallel::mclapply(X = 1:nJ, FUN = get_augment_diag, 
                                left = info_right, right = info_left, 
                                info_inv = info_inv, mc.cores = cores)
  aug_vec <- as.vector(unlist(aug_res))
  return(Y_tilde + aug_vec/2)
}

