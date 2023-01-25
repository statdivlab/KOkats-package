#' Generate values to augment Y
#' Generate the diagonal values of the augmentation matrix, a helper function for \code{generate_Y_tilde_plus}
#' 
#' @param i The index of the diagonal. 
#' @param left Left matrix to multiply.  
#' @param right Right matrix to multiply. 
#' @param info_inv Inverse of info matrix, middle matrix to multiply.
#'  
#' @return The value of the augmentation matrix for the ith diagonal.
#' 
#' @export
get_augment_diag <- function(i, left, right, info_inv) {
  return(as.numeric(left[i, ] %*% info_inv %*% right[, i]))
}
