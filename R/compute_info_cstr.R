#' Compute information matrix with constraint
#' Compute the information matrix for all parameters for given values of B and z with constraint.
#' 
#' @param X Design matrix with n rows (for samples) and p columns (for covariates).
#' @param B Beta parameter matrix with p rows (for covariates) and J columns (for categories).
#' @param z z parameter vector of length n.
#' @param constraint Should be "scc" for single category constraint, "mc" for mean constraint, or "msc" for mean 
#' over a subset constraint.
#' @param constraint_cat Category to constrain coefficients to equal the negative sum of all other categories.
#' @param subset_j Indices of categories to include in constraint for the mean over a subset constraint.
#' 
#' @return The square \code{p*J + n} information matrix, evaluated at the values of B and z with the constraint.
#' 
#' @export
compute_info_cstr <- function(X, B, z, constraint, constraint_cat, subset_j) {
  
  J <- ncol(B)
  n <- nrow(X)
  p <- nrow(B)
  # single category constraint
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
  # mean constraint
  } else if (constraint == "mc") {
    upd_B_alt_ind <- (1:(p*J))[-c(get_theta_ind(j = constraint_cat, k = 1:p, p = p))]
    z_mat <- matrix(z, ncol = 1) %*% matrix(1, ncol = J, nrow = 1)
    log_means <- X %*% B + z_mat
    info <- matrix(0, nrow = p*J + n, ncol = p*J + n)
    # info block B by B
    J_B <- matrix(NA, nrow = p*(J - 1), ncol = p*(J - 1))
    j_vec <- rep(1:J, each = p)[upd_B_alt_ind]
    k_vec <- rep(1:p, J)[upd_B_alt_ind]
    for (r in 1:(p*(J - 1))) {
      for (c in 1:r) {
        val <- sum(-X[, k_vec[r]]*X[, k_vec[c]]*exp(log_means[, constraint_cat]))
        if (j_vec[r] == j_vec[c]) {
          val <- val + sum(-X[, k_vec[r]]*X[, k_vec[c]] *
                             exp(log_means[, j_vec[r]]))
        }
        J_B[r, c] <- val
        J_B[c, r] <- val
      }
    }
    info[upd_B_alt_ind, upd_B_alt_ind] <- -J_B
    # info blocks B by z and z by B
    for (r in 1:(p*(J - 1))) {
      for (c in 1:n) {
        val <- -X[c, k_vec[r]]*exp(log_means[c, constraint_cat]) + 
          X[c, k_vec[r]]*exp(log_means[c, j_vec[r]])
        info[upd_B_alt_ind[r], (p*J + c)] <- val
        info[(p*J + c), upd_B_alt_ind[r]] <- val
      }
    }
    # info block z by z 
    for (i in 1:n) {
      info[p*J + i, p*J + i] <- sum(exp(log_means[i, ]))
    }
  # mean over a subset constraint
  } else {
    upd_B_alt_ind <- (1:(p*J))[-c(get_theta_ind(j = constraint_cat, k = 1:p, p = p))]
    z_mat <- matrix(z, ncol = 1) %*% matrix(1, ncol = J, nrow = 1)
    log_means <- X %*% B + z_mat
    info <- matrix(0, nrow = p*J + n, ncol = p*J + n)
    # info block B by B
    J_B <- matrix(NA, nrow = p*(J - 1), ncol = p*(J - 1))
    j_vec <- rep(1:J, each = p)[upd_B_alt_ind]
    k_vec <- rep(1:p, J)[upd_B_alt_ind]
    for (r in 1:(p*(J - 1))) {
      for (c in 1:r) {
        val <- 0
        if (j_vec[r] %in% subset_j & j_vec[c] %in% subset_j) {
          val <- val + sum(-X[, k_vec[r]]*X[, k_vec[c]] * exp(log_means[, constraint_cat]))
        }
        if (j_vec[r] == j_vec[c]) {
          val <- val + sum(-X[, k_vec[r]]*X[, k_vec[c]] * exp(log_means[, j_vec[r]]))
        }
        J_B[r, c] <- val
        J_B[c, r] <- val
      }
    }
    info[upd_B_alt_ind, upd_B_alt_ind] <- -J_B
    # info blocks B by z and z by B
    for (r in 1:(p*(J - 1))) {
      for (c in 1:n) {
        val <-  X[c, k_vec[r]]*exp(log_means[c, j_vec[r]])
        if (j_vec[r] %in% subset_j) {
          val <- val - X[c, k_vec[r]]*exp(log_means[c, constraint_cat])
        }
        info[upd_B_alt_ind[r], (p*J + c)] <- val
        info[(p*J + c), upd_B_alt_ind[r]] <- val
      }
    }
    # info block z by z 
    for (i in 1:n) {
      info[p*J + i, p*J + i] <- sum(exp(log_means[i, ]))
    }
  }
  
  # set constraint category part of info matrix to NA since these are not valid parameters
  info[constraint_cat, ] <- NA
  info[, constraint_cat] <- NA
  return(info)
}
