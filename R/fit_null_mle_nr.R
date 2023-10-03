#' Estimate under the null model
#' Estimate B and z parameters through block coordinate descent to maximize the unconstrained log likelihood function under a simple null hypothesis and given constraint.
#'
#' @param formula_rhs The right hand side of a formula specifying which covariates to include in the model, must be used with the \code{covariate_data} parameter or replaced by the \code{X} parameter.
#' @param covariate_data A data frame including all covariates given in \code{formula_rhs}, can be replaced by design matrix \code{X}.
#' @param Y An outcome matrix with n rows (for samples) and J columns (for categories) containing abundance data.
#' @param X Design matrix with n rows (for samples) and p columns (for covariates).
#' @param B Optional initial parameter estimate for B matrix.
#' @param constraint Should be "scc" for single category constraint, "mc" for mean constraint, or "msc" for mean 
#' over a subset constraint.
#' @param constraint_cat Category to constrain coefficients to equal the negative sum of all other categories.
#' @param subset_j Indices of categories to include in constraint for the mean over a subset constraint.
#' @param null_k Coefficient of the covariate in design matrix that is set to \code{0} under the null hypothesis.
#' @param null_j Category for which covariate \code{k} is set to \code{0} under the null hypothesis.
#' @param tolerance The tolerance used to stop the algorithm when log likelihood values are within \code{tolerance} of each other.
#' @param tolerance_nr The tolerance used to stop the Newton-Raphon method.
#' @param use_tolerance If \code{FALSE}, will run until \code{maxit} regardless of convergence.
#' @param maxit The maximum number of iterations of the coordinate descent algorithm.
#' @param maxit_nr The maximum number of iterations of the Newton-Raphson algorithm within the coordinate descent.
#' @param ncores The desired number of cores to optimize block of B parameters in parallel. If not provided, an appropriate number will be chosen for your machine.
#' @param solve_option
#'
#' @return A list including values of the log likelihood, the B matrix, and the z vector at each iteration.
#'
#' @examples
#' J <- 50
#' n <- 40
#' X <- cbind(1, rep(c(0, 1), each = n/2))
#' z <- rnorm(n) + 8
#' b0 <- rnorm(J)
#' b1 <- rnorm(J)
#' b <- rbind(b0, b1)
#' Y <- matrix(NA, ncol = J, nrow = n)
#' 
#' for (i in 1:n) {
#'  for (j in 1:J) {
#'    temp_mean <- exp(X[i, , drop = FALSE] %*% b[, j, drop = FALSE] + z[i])
#'    Y[i,j] <- rpois(1, lambda = temp_mean)
#'  }
#' }
#' 
#' null_mle <- fit_null_mle(Y = Y, X = X, ncores = 2, null_k = 2, null_j = 2,
#'                           constraint = "scc", solve_option = "solve")
#' 
#' @export
fit_null_mle_nr <- function(formula_rhs = NULL, Y, X = NULL, covariate_data = NULL, B = NULL,
                            constraint, constraint_cat = 1, subset_j = NULL, null_k = NULL, null_j = NULL,  
                            tolerance = 1e-10, tolerance_nr = 1e-10, use_tolerance = TRUE,  
                            maxit = 1000, maxit_nr = 1000, ncores = NULL, solve_option) {
  
  # check for valid constraint
  if (!(constraint %in% c("scc", "mc", "msc"))) {
    stop("Please provide a valid constraint, either 'scc', 'mc', or 'msc'.")
  }
  # check for requirements for mean over a subset constraint
  if (constraint == "msc") {
    if (is.null(subset_j)) {
      stop("If using the 'msc' constraint, please provide 'subset_j'.")
    } else {
      if (!(constraint_cat %in% subset_j)) {
        stop("Please choose a constraint category that is part of 'subset_j'.")
      }
    }
  }
  
  # check that data has been given with either X matrix or formula and covariate data
  if(!is.null(formula_rhs)){
    if(is.null(covariate_data)){
      stop("If formula_rhs is provided, covariates named in formula must be
           provided inside covariate_data.")
    }
    if(!is.data.frame(covariate_data)){
      stop("Argument covariate_data must be a data frame.")
    }
    if(!inherits(formula_rhs,"formula")){
      formula_rhs <- stats::as.formula(formula_rhs)
    }
    X <- stats::model.matrix(formula_rhs,covariate_data)
  }
  
  if (is.null(null_j) | is.null(null_k)) {
    h0 <- FALSE
  } else {
    h0 <- TRUE
  }
  
  # set up initial hyperparameters and parameters 
  p <- ncol(X)
  if (p < 2) {
    stop("X must contain an intercept column and at least one other (linearly
          independent) column")
  }
  J <- ncol(Y)
  n <- nrow(X)
  
  # find optimal parameters without null hypothesis constraint 
  if (is.null(B)) {
    if (is.null(B)) {
      if (constraint == "scc") {
        initial_constr <- function(x) {x[constraint_cat]}
      } else if (constraint == "mc") {
        initial_constr <- function(x) {mean(x)}
      } else {
        initial_constr <- function(x) {mean(x[subset_j])}
      }
      res <- fit_bcd_unconstrained(Y = Y, X = X, tolerance = tolerance, maxit = maxit,
                                   maxit_glm = NULL, ncores = ncores, 
                                   constraint_fn = initial_constr)
      z <- res$final_z
      B <- res$final_B
    } else {
      if (constraint == "scc") {
        B <- B - B[, constraint_cat]
      } else if (constraint == "mc") {
        B <- B - rowMeans(B)
      } else {
        B <- B - rowMeans(B[, subset_j])
      }
      z <- update_z(Y, X, B)
    }
  }
  
  # set null hypothesis constraint 
  if (h0) {
    B[null_k, null_j] <- 0
    if (constraint == "mc") {
      null_sum <- sum(B[null_k, ])
      B[null_k, -null_j] <- B[null_k, -null_j] - null_sum/(J - 1)
    } else if (constraint == "msc") {
      if (null_j %in% subset_j) {
        null_sum <- sum(B[null_k, subset_j])
        B[null_k, -null_j] <- B[null_k, -null_j] - null_sum/(length(subset_j) - 1)
      }
    }
  }
  
  z <- update_z(Y, X, B)
  f0 <- compute_loglik(Y, X, B, z)
  first_lik <- f0 
  
  # update B and z parameters through iteration 
  t <- 1 
  f_old <- -Inf
  f_new <- f0
  
  # save log likelihood values and beta values 
  lik_vec <- rep(NA, maxit)
  B_array <- array(NA, dim = c(nrow(B), ncol(B), maxit))
  z_array <- array(NA, dim = c(length(z), maxit))
  score_mat <- matrix(NA, nrow = (p*(J - 1) + n - 1), ncol = maxit)
  
  B_upd_ind <- (1:(p*J))[-c(get_theta_ind(j = constraint_cat, k = 1:p, p = p),
                            get_theta_ind(j = null_j, k = null_k, p = p))]
  z_ind <- (p*J + 1):(p*J + n)
  upd_ind <- c(B_upd_ind, z_ind) 
  
  # update parameter values via Newton-Raphson update 
  B_old <- as.vector(B)[B_upd_ind]
  B_new <- as.vector(B)[B_upd_ind]
  z_old <- z
  z_new <- z 
  t <- 1
  if (constraint == "msc") {
    subset_j_no_c <- subset_j[subset_j != constraint_cat]
  }
  
  # set initial parameter vector 
  B_full <- rep(0, p*J)
  # fill in B values that were updated 
  B_full[B_upd_ind] <- B_new
  B_mat <- matrix(B_full, nrow = p, ncol = J)
  # fill in B values that are predetermined by constraint
  if (constraint == "scc") {
    B_mat[, constraint_cat] <- 0
  } else if (constraint == "mc") {
    B_mat[, constraint_cat] <- -rowSums(B_mat[, -constraint_cat])
  } else {
    B_mat[, constraint_cat] <- -rowSums(B_mat[, subset_j_no_c]) 
  }
  #z_new <- update_z(Y, X, B_mat)
  # set null B to 0 (Sarah check - this should never change from 0)
  B_mat[null_k, null_j] <- 0 
  
  if (use_tolerance) {
    condition <- t == 1 || (f_new - f_old)/f_old > tolerance & t < maxit
  } else {
    condition <- t < maxit
  }
  while (condition) {
    
    score <- compute_scores_cstr(X = X, Y = Y, B = B_mat, z = z_new, 
                                 constraint = constraint, 
                                 constraint_cat = constraint_cat,
                                 subset_j = subset_j)[upd_ind]
    score_deriv <- -compute_info_cstr(X = X, B = B_mat, z = z_new, 
                                      constraint = constraint,
                                      constraint_cat = constraint_cat,
                                      subset_j = subset_j)[upd_ind, upd_ind]
    B_old <- B_new
    z_old <- z_new 
    theta_old <- c(B_old, z_old)
    if (solve_option == "chol") {
      theta_new <- theta_old + chol2inv(chol(-score_deriv)) %*% score
    } else if (solve_option == "solve") {
      if (constraint == "scc") {
        tmp <- Matrix::solve(score_deriv, -score, sparse = TRUE)
      } else {
        tmp <- Matrix::solve(score_deriv, -score)
      }
      theta_new <- tmp + theta_old
    } else if (solve_option == "limSolve") {
      if (constraint == "scc") {
        # to get this to work need to get this into the correct form 
        # could be faster - especially when solve fails because there are
        # multiple possible solutions
        # tmp <- limSolve::Solve.block(score_deriv, -score)
        tmp <- limSolve::Solve(as.matrix(score_deriv), -score)
      } else {
        tmp <- limSolve::Solve(score_deriv, -score)
      }
      
      theta_new <- tmp + theta_old
    } else {
      stop("Your option is not a valid solving option.")
    }
    B_new <- theta_new[upd_ind %in% B_upd_ind]
    z_new <- theta_new[upd_ind %in% z_ind]
    
    # put B back into matrix
    B_full[B_upd_ind] <- B_new
    B_mat <- matrix(B_full, nrow = p, ncol = J)
    # fill in B values that are predetermined by constraint
    if (constraint == "scc") {
      B_mat[, constraint_cat] <- 0
    } else if (constraint == "mc") {
      B_mat[, constraint_cat] <- -rowSums(B_mat[, -constraint_cat])
    } else {
      B_mat[, constraint_cat] <- -rowSums(B_mat[, subset_j_no_c]) 
    }
    #z_new <- update_z(Y, X, B_mat)
    # set null B to 0 (Sarah check - this should never change from 0)
    #B_mat[null_k, null_j] <- 0 
    
    # update t and likelihood value
    lik_vec[t] <- compute_loglik(Y, X, B_mat, z_new)
    scores <- compute_scores_cstr(X = X, Y = Y, B = B_mat, z = z_new, 
                                  constraint = constraint, 
                                  constraint_cat = constraint_cat,
                                  subset_j = subset_j)
    score_mat[, t] <- scores[upd_ind]
    B_array[, , t] <- B_mat
    z_array[, t] <- z_new
    f_old <- f_new
    f_new <- lik_vec[t]
    t <- t + 1
    
    if (use_tolerance) {
      condition <- (f_new - f_old)/f_old > tolerance & t < maxit
    } else {
      condition <- t < maxit
    }
  }
  
  final_B <- B_array[, , (t-1)]
  final_z <- update_z(Y, X, final_B)
  final_lik <- compute_loglik(Y, X, final_B, final_z)
  
  # get rid of NA values if algorithm finished before maxit
  converged <- FALSE
  if (t - 1 < maxit) {
    lik_vec <- lik_vec[1:(t - 1)]
    B_array <- B_array[, , 1:(t - 1)]
    z_array <- z_array[, 1:(t - 1)]
    score_mat <- score_mat[, 1:(t - 1)]
    converged <- TRUE
  }
  return(list(likelihood = lik_vec, B = B_array, z = z_array,
              final_B = final_B, final_z = final_z, 
              first_lik = first_lik, score_mat = score_mat,
              converged = converged))
}
