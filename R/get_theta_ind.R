#' Switch indices
#' Go from index in Beta matrix to theta vector. 
#' 
#' @param j Category or vector of categories.
#' @param k Covariate index or vector of covariate indices. 
#' @param p Number of covariates. 
#' 
#' @return The index of Beta_k^j in the theta vector. 
#' 
#' @export
get_theta_ind <- function(j, k, p) {
  return((j - 1)*p + k)
}
