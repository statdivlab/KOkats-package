#' Generate Y tilde
#' Compress Y matrix into vector Y tilde.
#' 
#' @param Y Matrix of Y values, n rows by J columns 
#' 
#' @return Y_tilde, a vector of length n*J
#' 
#' @export 
generate_Y_tilde <- function(Y) {
  return(as.vector(t(Y)))
}