#' Estimate variance of single score value
#' Calculate the empirical variance of the score for the parameter tested.
#'
#' @param null_k Coefficient of the covariate in design matrix that is set to \code{0} under the null hypothesis.
#' @param null_j Category for which covariate \code{k} is set to \code{0} under the null hypothesis.
#' @param upd_ind Indices of parameter vector, not including parameters from constraint row or parameter implicated in null hypothesis.
#' @param p Number of columns in design matrix. 
#' @param D Estimated empirical variance of score. 
#' @param info Information matrix. 
#'
#' @return The empirical variance of the score for a single value.
#'
#' @export
var_single_score <- function(null_j, null_k, upd_ind, p, D, info) {
  # calculate score statistic 
  null_ind <- get_theta_ind(null_j, null_k, p)
  # get variance of score 
  D_11 <- D[null_ind, null_ind]
  I_12 <- info[null_ind, upd_ind]
  I_22_inv <- chol2inv(chol(info[upd_ind, upd_ind]))
  D_12 <- D[null_ind, upd_ind]
  D_22 <- D[upd_ind, upd_ind]
  score_var <- D_11 - I_12 %*% I_22_inv %*% D_12 - 
    D_12 %*% I_22_inv %*% I_12 + 
    I_12 %*% I_22_inv %*% D_22 %*% I_22_inv %*% I_12
  return(score_var)
}
