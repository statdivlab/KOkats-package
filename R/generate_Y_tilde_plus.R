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
generate_Y_tilde_plus <- function(Y_tilde, X_tilde, W) {
  info <- t(X_tilde) %*% W %*% X_tilde
  #info_chol <- Matrix::chol(info) 
  # grabbing directly from radEmu - come back to think more about this later!!!
  # info_chol <- suppressWarnings(
  #   try(Matrix::chol(info),silent = TRUE)
  # )
  # if(!inherits(info_chol,"Matrix")){
  #   perturb <- 1e-8
  #   while(!("matrix" %in% class(info_chol))){
  #     message("Reattempting Cholesky decomposition")
  #     info_chol <- suppressWarnings(
  #       try(Matrix::chol(info +
  #                          Matrix::Diagonal(
  #                            x = rep(perturb*max(info),
  #                                    nrow(info)),
  #                          )),silent = TRUE)
  #     )
  #     perturb <- perturb*10
  #   }
  # }
  # chol_inv <- Matrix::solve(info_chol)
  info_inv <- MASS::ginv(info)
  W_half <- sqrt(W)
  # augmented portion 
  aug <- diag(W_half %*% X_tilde %*% 
                info_inv %*% 
                t(X_tilde) %*% W_half)
  return(Y_tilde + aug/2)
}
