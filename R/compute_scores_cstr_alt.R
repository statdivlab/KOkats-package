#' Compute score vector with constraint
#' Compute the score equations for all parameters for given values of B and z with constraint.
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
#' @return The value of each score equation, evaluated at the values of B and z with the constraint.
#' 
#' @export
compute_scores_cstr_alt <- function(X, Y, B, z, constraint, constraint_cat, subset_j) {
  n <- nrow(X)
  p <- ncol(X)
  J <- ncol(Y)
  scores <- rep(NA, p*J + n) 
  z_mat <- matrix(z, ncol = 1) %*% matrix(1, ncol = J, nrow = 1)
  log_means <- X %*% B + z_mat
  js <- (1:J)[-constraint_cat]
  for (j in js) {
    for (k in 1:p) {
      if (constraint == "scc") {
        #scores[(j - 1)*p + k] <-  sum(X[, k] %*% (Y[, j] - exp(log_means[, j])))
        scores[(j - 1)*p + k] <-  sum(X[, k] %*% (Y[, j]) - 
                                        sum(X[, k] %*% exp(log_means[, j])))
      } else if (constraint == "mc") {
        # scores[(j - 1)*p + k] <- sum(X[, k] %*% (Y[, j] - exp(log_means[, j]) -
        #                                            Y[, constraint_cat] +
        #                                            exp(log_means[, constraint_cat])))
        scores[(j - 1)*p + k] <- sum(X[, k] %*% Y[, j]) - 
          sum(X[, k] %*% exp(log_means[, j])) - sum(X[, k] %*% Y[, constraint_cat]) +
          sum(X[, k] %*% exp(log_means[, constraint_cat]))
      } else {
        if (j %in% subset_j) {
          # scores[(j - 1)*p + k] <- sum(X[, k] %*% (Y[, j] - exp(log_means[, j]) -
          #                                            Y[, constraint_cat] +
          #                                            exp(log_means[, constraint_cat])))
          scores[(j - 1)*p + k] <- sum(X[, k] %*% Y[, j]) - 
            sum(X[, k] %*% exp(log_means[, j])) - sum(X[, k] %*% Y[, constraint_cat]) + 
            sum(X[, k] %*% exp(log_means[, constraint_cat]))
        } else {
          #scores[(j - 1)*p + k] <-  sum(X[, k] %*% (Y[, j] - exp(log_means[, j])))
          scores[(j - 1)*p + k] <-  sum(X[, k] %*% Y[, j]) - 
            sum(X[, k] %*% exp(log_means[, j]))
        }
      }
    }
  }
  for (i in 1:n) {
    scores[p*J + i] <- sum(Y[i, ]) - sum(exp(log_means[i, ])) 
  }
  return(scores)
}

