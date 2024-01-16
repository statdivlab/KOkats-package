#' Get reduced bias B estimates and augmentation values 
#' Implement algorithm from Kosmidis & Firth (2011) to get reduced bias B estimates and augmentation values
#'
#' @param X Design matrix with n rows (for samples) and p columns (for covariates), should have a leading intercept column of \code{1}s.
#' @param B Beta matrix 
#' @param z z vector 
#' @param cc Constraint category 
#' @param Y Data matrix 
#' @param tolerance The tolerance used to stop the algorithm when log likelihood values are within \code{tolerance} of each other.
#' @param maxit The maximum number of iterations of the coordinate descent algorithm.
#'
#' @return B estimates and augmentation values 
#'
#' @examples
#' n <- 50
#' J <- 10
#' X <- cbind(1, rnorm(n))
#' #X <- cbind(1, rep(0:1, n/2))
#' b0 <- c(rnorm(J - 1), 0)
#' b1 <- c(rnorm(J - 1), 0)
#' B <- rbind(b0, b1)
#' z_init <- rnorm(n) 
#' Y <- matrix(NA, ncol = J, nrow = n)
#' 
#' set.seed(2)
#' for (i in 1:n) {
#'  for (j in 1:J) {
#'    temp_mean <- exp(X[i, , drop = FALSE] %*% B[, j, drop = FALSE] + z_init[i])
#'    Y[i,j] <- rpois(1, lambda = temp_mean)
#'  }
#' }
#' 
#' 
#' res <- reduced_bias_B_alg(X, cc = NULL, Y, tolerance = 1e-6, maxit = 10000)
#' 
#' @export
reduced_bias_B_alg <- function(X, cc, Y, tolerance = 1e-2, maxit = 1000) {
  
  # hyperparameters 
  n <- nrow(X)
  p <- ncol(X)
  J <- ncol(Y)
  if (is.null(cc)) {
    # set the constraint category as the final category if not specified 
    cc <- J
  }
  
  # initialization values 
  B_old <- matrix(-Inf, nrow = p, ncol = J)
  B_new <- matrix(0, nrow = p, ncol = J)
  t <- 1
  cores <- parallel::detectCores() - 1
  
  # until tolerance 
  while (t == 1 | 
         (sum(((B_new[, -cc] - B_old[, -cc])/B_old[, -cc])^2) > tolerance & t < maxit)) {
   
    # get restricted z's (under restriction that sum_j Y_ij = sum_j mu_ij for all i)
    z_res <- update_z(Y, X, B_new)
    
    # use current B's and restricted z's to calculate H matrix 
    augs <- get_hat(X, B_new, z_res, cc)
    Y_new <- Y + augs/2
    
    # fit model with ML using adjusted responses Y_ij + h_ijj/2 to get new B and z estimates 
    res <- fit_bcd_unconstrained(Y = Y_new, X = X, tolerance = tolerance, maxit = maxit,
                                 maxit_glm = NULL, ncores = cores, 
                                 constraint_fn = function(x) {x[cc]}) 
    
    # iterate 
    B_old <- B_new
    B_new <- res$final_B
    t <- t + 1
  }
  
  # return final B and last h_ijj/2 values
  return(list(B = B_new, augs = augs/2, t = t))
  
}
