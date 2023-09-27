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
#' res <- null_score_var_alt(Y = Y, X = X, B = b, z = z, 
#'                       constraint = "scc", constraint_cat = 1)
#' 
#' @export
null_score_var_alt <- function(Y, X, B, z, constraint, constraint_cat) {
  
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
      if (constraint == "scc") {
        score[ind] <- X[i, k_vec[ind]] * (Y[i, j] - exp(X[i, ] %*% B[, j] + z[i]))
      } else {
        score[ind] <- X[i, k_vec[ind]] * (-Y[i, constraint_cat] + Y[i, j] +
                                            exp(X[i, ] %*% B[, constraint_cat] + z[i]) - exp(X[i, ] %*% B[, j] + z[i]))
      }
    }
    for (ind in 1:n) {
      if (ind == i) {
        if (constraint == "scc") {
          score[ind + p*(J - 1)] <- Y[i, constraint_cat] - exp(z[i]) + 
            Y[i, j] - sum(exp(X[i, ] %*% B[, js] + z[i]))
        } else {
          score[ind + p*(J - 1)] <- Y[i, constraint_cat] - exp(X[i, ] %*% B[, constraint_cat] + z[i]) + 
            Y[i, j] - sum(exp(X[i, ] %*% B[, js] + z[i]))
        }
      }
    }
    D <- D + score %*% t(score)
  }
  
  full_D <- matrix(0, nrow = p*J + n, ncol = p*J + n)
  constraint_ind <- get_theta_ind(constraint_cat, 1:p, p)
  full_D[-constraint_ind, -constraint_ind] <- D
  
  # calculate in original way
  org_D <- matrix(0, nrow = p*(J - 1) + n, ncol = p*(J - 1) + n)
  for (i in 1:n) {
    for (j in js) {
      score <- rep(0, p*(J - 1) + n)
      for (ind in 1:(p*(J - 1))) {
        if (j_vec[ind] == j) {
          if (constraint == "scc") {
            score[ind] <- Y[i, j] * X[i, k_vec[ind]] - 
              X[i, k_vec[ind]] * exp(X[i, ] %*% B[, j] + z[i])
          } else {
            #score[ind] <- -Y[i, 1] * X[i, k_vec[ind]] / (J - 1) + 
            score[ind] <- -Y[i, constraint_cat] * X[i, k_vec[ind]] + 
              1/(J - 1) * X[i, k_vec[ind]] * exp(X[i, ] %*% B[, constraint_cat] + z[i]) +
              Y[i, j] * X[i, k_vec[ind]] - X[i, k_vec[ind]] * 
              exp(X[i, ] %*% B[, j] + z[i]) 
          }
        } else {
          if (constraint == "mc") {
            #score[ind] <- -Y[i, 1] * X[i, k_vec[ind]] / (J - 1) + 
            #  1/(J - 1) * X[i, k_vec[ind]] * exp(X[i, ] %*% B[, 1] + z[i])
            score[ind] <- 1/(J - 1) * X[i, k_vec[ind]] * 
              exp(X[i, ] %*% B[, constraint_cat] + z[i])
          }
        }
      }
      for (ind in 1:n) {
        if (ind == i) {
          if (constraint == "scc") {
            score[ind + p*(J - 1)] <- 1/(J - 1) * (Y[i, constraint_cat] - exp(z[i])) + 
              Y[i, j] - exp(X[i, ] %*% B[, j] + z[i])
          } else {
            score[ind + p*(J - 1)] <- 1/(J - 1) * 
              (Y[i, constraint_cat] - exp(X[i, ] %*% B[, constraint_cat] + z[i])) +
              Y[i, j] - exp(X[i, ] %*% B[, j] + z[i])
          }
        }
      }
      org_D <- org_D + score %*% t(score)
    }
  }
  org_full_D <- matrix(0, nrow = p*J + n, ncol = p*J + n)
  constraint_ind <- get_theta_ind(constraint_cat, 1:p, p)
  org_full_D[-constraint_ind, -constraint_ind] <- org_D
  
  # return empirical score variance
  return(list(D = full_D, orig_D = org_full_D))
}