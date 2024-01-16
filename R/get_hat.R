#' Get hat matrix 
#' Calculate diagonals of hat matrix for Poisson log linear model 
#'
#' @param X Design matrix with n rows (for samples) and p columns (for covariates), should have a leading intercept column of \code{1}s.
#' @param B Beta matrix 
#' @param z z vector 
#' @param cc Constraint category 
#'
#' @return The diagonals of the hat matrix for the Poisson log linear model with design matrix X, parameters B, offsets z, and cc category used as the baseline category.
#'
#' @examples
#' n <- 50
#' J <- 10
#' cc <- 10
#' #X <- cbind(1, rnorm(n))
#' X <- cbind(1, rep(0:1, n/2))
#' b0 <- c(rnorm(J - 1), 0)
#' b1 <- c(rnorm(J - 1), 0)
#' B <- rbind(b0, b1)
#' z_init <- rnorm(n) + 8
#' Y <- matrix(NA, ncol = J, nrow = n)
#' 
#' for (i in 1:n) {
#'  for (j in 1:J) {
#'    temp_mean <- exp(X[i, , drop = FALSE] %*% B[, j, drop = FALSE] + z_init[i])
#'    Y[i,j] <- rpois(1, lambda = temp_mean)
#'  }
#' }
#' z <- update_z(Y, X, B)
#' get_hat(X, B, z, cc)
#' 
#' @export
get_hat <- function(X, B, z, cc = NULL) {
  
  J <- ncol(B)
  n <- nrow(X)
  p <- ncol(X)
  
  if (is.null(cc)) {
    # set the constraint category as the final category if not specified 
    cc <- J
  }
  
  # make matrix to save values in 
  hat_diags <- matrix(data = NA, nrow = n, ncol = J)
  
  # get W diagonal matrix 
  X_prime <- generate_X_prime(X, J)
  A_prime <- generate_A_prime(X_prime, J)
  B_prime <- as.vector(B)
  theta <- generate_theta(B_prime, z)
  W <- generate_W(A_prime, theta)
  
  # get information matrix (ignore cc)
  param_j_ind <- c(sort(rep(1:J, p)), rep(0, n))
  param_cstr_ind <- param_j_ind != cc
  full_info <- Matrix::crossprod(A_prime, W) %*% A_prime 
  info <- full_info[param_cstr_ind, param_cstr_ind]
  
  # get useful quantities
  ident <- diag(J - 1)
  one_vec <- rep(1, J - 1)
  
  # loop over samples 
  for (i in 1:n) {
    # construct G_i matrix 
    G_i <- kronecker(ident, t(X[i, ]))
    
    # construct Z_i matrix 
    e_i <- rep(0, n)
    e_i[i] <- 1 
    partA <- cbind(G_i, one_vec %x% t(e_i))
    partB <- c(rep(0, p*(J - 1)), t(e_i))
    if (cc == 1) {
      Z_i <- rbind(partB, partA)
    } else if (cc == J) {
      Z_i <- rbind(partA, partB)
    } else {
      Z_i <- rbind(partA[1:(cc - 1),], partB, partA[cc:(J - 1), ])
    }
    
    # note! Z_i is equivalent to A_prime[i_ind, -cc param ind] 
    # remove redundancy from code later on!
    
    # construct H_i matrix 
    i_ind <- (1 + (i - 1)*J):(J + (i - 1)*J)
    H_i <- Z_i %*% solve(info) %*% t(Z_i) %*% W[i_ind, i_ind]
    
    # save H_i diagonals 
    hat_diags[i, ] <- diag(as.matrix(H_i))
  }
  
  # return values 
  return(hat_diags)
  
}
