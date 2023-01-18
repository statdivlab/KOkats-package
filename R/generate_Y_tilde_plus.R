#' Generate Y plus 
#' Generate an augmented outcome matrix to use with Firth penalized likelihood.
#' 
#' @param Y_tilde_plus An augmented outcome vector of length \code{nJ}.
#' @param J The number of KOs. 
#'  
#' @return Y_plus, the augmented outcome matrix.
#' 
#' @export
generate_Y_plus <- function(Y_tilde_plus, J) {
  n <- length(Y_tilde_plus)/J
  Y_plus <- matrix(Y_tilde_plus, nrow = n, ncol = J, byrow = TRUE) 
  return(Y_plus)
}