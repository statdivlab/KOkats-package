#' Generate values to augment Y
#' Generate the diagonal values of the augmentation matrix, a helper function for \code{generate_Y_tilde_plus}
#' 
#' @param i The index of the diagonal. 
#'  
#' @return The value of the augmentation matrix for the ith diagonal.
#' 
#' @export
get_augment_diag <- function(i) {
  return(as.numeric(left[i, ] %*% info_inv %*% right[, i]))
}