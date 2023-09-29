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
#' 
#' @return The value of each score equation, evaluated at the values of B and z with the constraint.
#' 
#' @export
compute_scores_cstr <- function(X, Y, B, z, constraint, constraint_cat = 1) {
  if (!(constraint %in% c("scc", "mc", "msc"))) {
    stop("Please provide a valid constraint, either 'scc', 'mc', or 'msc'.")
  }
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
          scores[(j - 1)*p + k] <-  sum(X[, k] %*% (Y[, j] - exp(log_means[, j])))
        } else if (constraint == "mc") {
          scores[(j - 1)*p + k] <- sum(X[, k] %*% (Y[, j] - exp(log_means[, j]) -
                                                     Y[, constraint_cat] +
                                                     exp(log_means[, constraint_cat])))
        }
    }
  }
  for (i in 1:n) {
    scores[p*J + i] <- sum(Y[i, ] - exp(log_means[i, ])) 
  }
  return(scores)
}

