#' Variance of score
#' Calculate the empirical variance of the score for the constrained likelihood.
#'
#' @param Y An outcome matrix with n rows (for samples) and J columns (for categories) containing abundance data.
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
compute_score_var_cstr <- function(Y, X, B, z, constraint, constraint_cat, subset_j) {
  
  # get hyperparameters 
  n <- nrow(X)
  p <- ncol(X)
  J <- ncol(Y)
  
  js <- (1:J)[-constraint_cat]
  j_vec <- c(rep(js, each = p), rep(0, n))
  k_vec <- c(rep(1:p, J-1), rep(0, n))
  z_mat <- matrix(z, ncol = 1) %*% matrix(1, ncol = J, nrow = 1)
  log_means <- X %*% B + z_mat
  D <- matrix(0, nrow = p*(J - 1) + n, ncol = p*(J - 1) + n)
  for (i in 1:n) {
    score <- rep(0, p*(J - 1) + n)
    for (ind in 1:(p*(J - 1))) {
      j <- j_vec[ind]
      if (constraint == "scc") {
        score[ind] <- X[i, k_vec[ind]] * (Y[i, j] - exp(log_means[i, j]))
      } else if (constraint == "mc") {
        score[ind] <- X[i, k_vec[ind]] * (-Y[i, constraint_cat] + Y[i, j] +
                                            exp(log_means[i, constraint_cat]) - 
                                            exp(log_means[i, j]))
      } else {
        val <- X[i, k_vec[ind]] * (Y[i, j] - exp(log_means[i, j]))
        if (j %in% subset_j) {
          val <- val + X[i, k_vec[ind]] * (-Y[i, constraint_cat] + 
                                             exp(log_means[i, constraint_cat]))
        }
        score[ind] <- val
      }
    }
    for (ind in 1:n) {
      if (ind == i) {
        score[ind + p*(J - 1)] <- sum(Y[i, ] - exp(X[i, ] %*% B + z[i]))
      }
    }
    D <- D + score %*% t(score)
  }
  
  full_D <- matrix(NA, nrow = p*J + n, ncol = p*J + n)
  constraint_ind <- get_theta_ind(constraint_cat, 1:p, p)
  full_D[-constraint_ind, -constraint_ind] <- D
  
  # return empirical score variance
  return(full_D)
}

