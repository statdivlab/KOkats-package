#' Estimate variance of single score value
#' Calculate the empirical variance of the score for the parameter tested.
#'
#' @param null_ind Index of parameter vector corresponding with parameter being tested.
#' @param upd_ind Indices of parameter vector, not including parameters from constraint row or parameter implicated in null hypothesis.
#' @param p Number of columns in design matrix. 
#' @param D Estimated empirical variance of score. 
#' @param info Information matrix. 
#' @param robust If TRUE uses the D matrix and creates a robust test, if FALSE uses only the information matrix for a model-based test.
#'
#' @return The empirical variance of the score for a single value.
#'
#' @export
compute_var_single_score <- function(null_ind, upd_ind, p, D = NULL, info, robust = TRUE) {
  
  if (is.null(D) & robust) {
    stop("If doing a robust test, please provide D matrix.")
  }
  
  # get variance of score 
  if (!(robust)) {
    score_var <- info[null_ind, null_ind] - 
      info[null_ind, upd_ind] %*% chol2inv(chol(info[upd_ind, upd_ind])) %*% info[upd_ind, null_ind]
  } else {
    D_11 <- D[null_ind, null_ind]
    I_12 <- info[null_ind, upd_ind]
    I_22_inv <- chol2inv(chol(info[upd_ind, upd_ind]))
    D_12 <- D[null_ind, upd_ind]
    D_22 <- D[upd_ind, upd_ind]
    score_var <- D_11 - I_12 %*% I_22_inv %*% D_12 - 
      D_12 %*% I_22_inv %*% I_12 + 
      I_12 %*% I_22_inv %*% D_22 %*% I_22_inv %*% I_12
  }
  
  return(score_var)
}
