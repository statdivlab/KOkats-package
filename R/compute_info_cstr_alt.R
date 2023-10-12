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
compute_info_cstr_alt <- function(X, B, z, constraint, constraint_cat, subset_j) {
  
  # hyperparameters 
  J <- ncol(B)
  n <- nrow(X)
  p <- nrow(B)
  # observations and parameters without the constraint category 
  obs_j_ind <- rep(1:J, n)
  obs_cstr_ind <- obs_j_ind != constraint_cat
  param_j_ind <- sort(rep(1:J, p))
  param_cstr_ind <- c(param_j_ind != constraint_cat, rep(FALSE, n))
  # compute info components 
  B_prime <- as.vector(B)
  theta <- generate_theta(B_prime, z)
  X_prime <- generate_X_prime(X, J)
  A_prime <- generate_A_prime(X_prime, J)
  W <- generate_W(A_prime, theta)
  Ws <- generate_W_cstr(W = W, X = X, B = B, z = z, constraint = constraint, 
                        constraint_cat = constraint_cat, subset_j = subset_j)
  W_BB <- Ws$W_BB
  W_Bz <- Ws$W_Bz
  # generate blocks of info matrix
  info_BB <- Matrix::crossprod(A_prime[obs_cstr_ind, param_cstr_ind], W_BB) %*%
    A_prime[obs_cstr_ind, param_cstr_ind]
  info_Bz <- Matrix::crossprod(A_prime[obs_cstr_ind, param_cstr_ind], W_Bz) %*%
    A_prime[obs_cstr_ind, (J*p + 1):(J*p + n)]
  info_zz <- Matrix::crossprod(A_prime[, (J*p + 1):(J*p + n)], W) %*%
    A_prime[, (J*p + 1):(J*p + n)]
  info_zB <- MatrixExtra::t_shallow(info_Bz)
  info <- cbind(rbind(info_BB, info_zB), rbind(info_Bz, info_zz))
  # set constraint category part of info matrix to NA since these are not valid parameters
  constraint_ind <- get_theta_ind(constraint_cat, 1:p, p)
  info_na <- cbind(info[, 1:nrow(info) < constraint_ind[1]],
                   matrix(NA, nrow = nrow(info), ncol = p),
                   info[, 1:nrow(info) >= constraint_ind[1]])
  info_na <- rbind(info_na[1:nrow(info_na) < constraint_ind[1], ],
                   matrix(NA, nrow = p, ncol = ncol(info_na)),
                   info_na[1:nrow(info_na) >= constraint_ind[1], ])
  return(info_na)
}
