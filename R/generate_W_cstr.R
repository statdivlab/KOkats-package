#' Variance of score
#' Calculate the empirical variance of the score for the constrained likelihood.
#'
#' @param W W matrix from unconstrained likelihood, diagonal matrix with entry i, j giving expected Y_ij based on mean model.
#' @param X Design matrix with n rows (for samples) and p columns (for covariates).
#' @param B Beta parameter matrix with p rows (for covariates) and J columns (for categories).
#' @param z z parameter vector of length n.
#' @param constraint Should be "scc" for single category constraint, "mc" for mean constraint, or "msc" for mean 
#' over a subset constraint.
#' @param constraint_cat Category to constrain coefficients to equal the negative sum of all other categories.
#' @param subset_j Indices of categories to include in constraint for the mean over a subset constraint.
#'
#' @return A matrix giving the empirical variance of the score vector.
#' 
#' @export
generate_W_cstr <- function(W, X, B, z, constraint, constraint_cat, subset_j) {
  
  # hyperparameters
  n <- nrow(X)
  J <- ncol(B)
  obs_j_ind <- rep(1:J, n)
  obs_cstr_ind <- obs_j_ind != constraint_cat
  W_cstr <- W[obs_cstr_ind, obs_cstr_ind]
  
  if (constraint == "scc") {
    W_list <- list(W_BB = W_cstr, W_Bz = W_cstr)
  } else if (constraint == "mc") {
    # portion for beta by beta block of W matrix 
    cstr_cat_vals <- exp(X %*% B[, constraint_cat] + z)
    diag_list <- list()
    for (i in 1:n) {
      diag_list[[i]] <- matrix(data = cstr_cat_vals[i], nrow = J - 1, ncol = J - 1)
    }
    W_BB <- bdiag_m(diag_list) + W_cstr
    #vals <- lapply(1:n, function(x) {c(rep(0, (J - 1)*(x - 1)), 
    #                                   rep(cstr_cat_vals[x], (J - 1)), 
    #                                   rep(0, (J - 1)*(n - x)))})
    #new_vals <- lapply(vals, function(x) {rep(x, J - 1)})
    #W_BB <- matrix(unlist(new_vals), nrow(W_cstr), ncol(W_cstr), byrow = TRUE) + 
    #  W_cstr
    # portion for beta by z block of W matrix
    W_Bz <- W_cstr - diag(rep(cstr_cat_vals, each = (J - 1)))
    W_list <- list(W_BB = W_BB, W_Bz = W_Bz)
  } else {
    j_ind <- (1:J)[-constraint_cat]
    subset_j_ind <- j_ind %in% subset_j
    subset_row <- matrix(subset_j_ind, J - 1, J - 1, byrow = TRUE)
    subset_col <- matrix(subset_j_ind, J - 1, J - 1, byrow = FALSE)
    subset_mat <- subset_row * subset_col 
    # portion for beta by beta block of W matrix 
    cstr_cat_vals <- exp(X %*% B[, constraint_cat] + z)
    diag_list <- list()
    for (i in 1:n) {
      diag_list[[i]] <- cstr_cat_vals[i] * subset_mat 
    }
    W_BB <- bdiag_m(diag_list) + W_cstr
    # portion for beta by z block of W matrix
    W_Bz <- W_cstr - diag(rep(cstr_cat_vals, each = (J - 1)) * subset_j_ind)
    W_list <- list(W_BB = W_BB, W_Bz = W_Bz)
  }
  return(W_list)
}
