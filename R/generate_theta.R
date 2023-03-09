#' Generate theta
#' Put all B and z parameters together into one vector. 
#' 
#' @param B_vec vector of B values of length p*J. 
#' @param z vector of z values of length n.
#' 
#' @return theta, a vector of length p*J + n
#' 
#' @export 
generate_theta <- function(B_vec, z) {
  return(c(B_vec, z))
}