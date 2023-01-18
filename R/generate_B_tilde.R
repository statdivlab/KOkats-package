#' Generate B tilde
#' Compress B matrix into vector B tilde.
#' 
#' @param B Matrix of B values, p rows by J columns 
#' 
#' @return B_tilde, a vector of length p*J
#' 
#' @export 
generate_B_tilde <- function(B) {
  return(as.vector(t(B)))
}