#' Variance of score
#' Calculate the empirical variance of the score for the constrained likelihood.
#'
#' @param Y An outcome matrix with n rows (for samples) and J columns (for KOs) containing coverage data.
#' @param X Design matrix with n rows (for samples) and p columns (for covariates), should have a leading intercept column of \code{1}s, can be replaced by \code{formula_rhs} and \code{covariate_data}.
#' @param B B parameter matrix.
#' @param z z parameter vector. 
#' @param constraint Should be "scc" for single category constraint or "mc" for mean constraint.
#' @param constraint_cat Category to constrain coefficients (for single category the B parameters for this
#' category will equal 0, for mean the B parameters for this category will equal the negative sum of the others.
#'
#' @return The empirical variance of the score.
#'
#' @examples
#' X <- cbind(1, rep(c(0, 1), each = 20))
#' z <- rnorm(40) + 8
#' b0 <- rnorm(10)
#' b1 <- 1:10
#' b <- rbind(b0, b1)
#' Y <- matrix(NA, ncol = 10, nrow = 40)
#' 
#' for (i in 1:40) {
#'  for (j in 1:10) {
#'    temp_mean <- exp(X[i, , drop = FALSE] %*% b[, j, drop = FALSE] + z[i])
#'    Y[i,j] <- rpois(1, lambda = temp_mean)
#'  }
#' }
#' 
#' res <- null_score_var(Y = Y, X = X, B = b, z = z, 
#'                       constraint = "scc", constraint_cat = 1)
#' 
#' @export
null_score_var <- function(Y, X, B, z, constraint, constraint_cat) {
  
  # check for valid constraint 
  if (!(constraint %in% c("scc", "mc"))) {
    stop("Please enter 'scc' for a single category constraint or 'mc' for a mean constraint.")
  }
  
  # get hyperparameters 
  n <- nrow(X)
  p <- ncol(X)
  J <- ncol(Y)
  
  js <- (1:J)[-constraint_cat]
  j_vec <- c(rep(js, each = p), rep(0, n))
  k_vec <- c(rep(1:p, J-1), rep(0, n))
  D <- matrix(0, nrow = p*(J - 1) + n, ncol = p*(J - 1) + n)
  for (i in 1:n) {
    score <- rep(0, p*(J - 1) + n)
    for (ind in 1:(p*(J - 1))) {
      j <- j_vec[ind]
      if (constraint == "sc") {
        score[ind] <- X[i, k_vec[ind]] * (Y[i, j] - exp(X[i, ] %*% B[, j] + z[i]))
      } else {
        score[ind] <- X[i, k_vec[ind]] * (-Y[i, 1] + Y[i, j] +
          exp(X[i, ] %*% B[, 1] + z[i]) - exp(X[i, ] %*% B[, j] + z[i]))
      }
    }
    for (ind in 1:n) {
      if (ind == i) {
        if (constraint == "scc") {
          score[ind + p*(J - 1)] <- Y[i, 1] - exp(z[i]) + 
            Y[i, j] - sum(exp(X[i, ] %*% B[, js] + z[i]))
        } else {
          score[ind + p*(J - 1)] <- Y[i, 1] - exp(X[i, ] %*% B[, 1] + z[i]) + 
            Y[i, j] - sum(exp(X[i, ] %*% B[, js] + z[i]))
        }
      }
    }
    D <- D + score %*% t(score)
  }
  
  full_D <- matrix(0, nrow = p*J + n, ncol = p*J + n)
  constraint_ind <- get_theta_ind(constraint_cat, 1:p, p)
  full_D[-constraint_ind, -constraint_ind] <- D
  
  # return empirical score variance
  return(D)
}

