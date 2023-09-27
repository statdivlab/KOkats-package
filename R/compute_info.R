#' Compute information
#' Calculate the information matrix with MLEs of parameters.
#'
#' @param X Design matrix with n rows (for samples) and p columns (for covariates), should have a leading intercept column of \code{1}s, can be replaced by \code{formula_rhs} and \code{covariate_data}.
#' @param B B parameter matrix.
#' @param z z parameter vector. 
#' @param constraint Should be "scc" for single category constraint or "mc" for mean constraint.
#' @param constraint_cat Category to constrain coefficients (for single category the B parameters for this
#' category will equal 0, for mean the B parameters for this category will equal the negative sum of the others.
#'
#' @return The information matrix for a constrained likelihood. 
#'
#' @export
compute_info <- function(X, B, z, constraint, constraint_cat) {
  J <- ncol(B)
  if (constraint == "scc") {
    B_prime <- as.vector(B)
    theta <- generate_theta(B_prime, z)
    X_prime <- generate_X_prime(X, J)
    A_prime <- generate_A_prime(X_prime, J)
    W <- generate_W(A_prime, theta)
    W_half <- sqrt(W)
    info_left <- Matrix::crossprod(A_prime, W_half) 
    info_right <- W_half %*% A_prime
    info <- info_left %*% info_right
  } else {
    n <- nrow(X)
    p <- nrow(B)
    upd_B_alt_ind <- (1:(p*J))[-c(get_theta_ind(j = constraint_cat, k = 1:p, p = p))]
    info <- matrix(0, nrow = p*J + n, ncol = p*J + n)
    # info block B by B
    J_B <- matrix(NA, nrow = p*(J - 1), ncol = p*(J - 1))
    j_vec <- rep(1:J, each = p)[upd_B_alt_ind]
    k_vec <- rep(1:p, J)[upd_B_alt_ind]
    for (r in 1:(p*(J - 1))) {
      for (c in 1:r) {
        val <- sum(-X[, k_vec[r]]*X[, k_vec[c]]*exp(X %*% B[, constraint_cat] + z))
        if (j_vec[r] == j_vec[c]) {
          val <- val + sum(-X[, k_vec[r]]*X[, k_vec[c]]*
                             exp(X %*% B[, j_vec[r]] + z))
        }
        J_B[r, c] <- val
        J_B[c, r] <- val
      }
    }
    info[upd_B_alt_ind, upd_B_alt_ind] <- -J_B
    # info blocks B by z and z by B
    for (r in 1:(p*(J - 1))) {
      for (c in 1:n) {
        val <- -X[c, k_vec[r]]*exp(X[c, ] %*% B[, constraint_cat] + z[c]) + 
          X[c, k_vec[r]]*exp(X[c, ] %*% B[, j_vec[r]] + z[c])
        info[upd_B_alt_ind[r], (p*J + c)] <- val
        info[(p*J + c), upd_B_alt_ind[r]] <- val
      }
    }
    # info block z by z 
    for (i in 1:n) {
      info[p*J + i, p*J + i] <- sum(exp(X[i, ] %*% B + z[i]))
    }
  }
  return(info)
}
