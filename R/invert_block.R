#' Invert matrix
#' Invert a block of the block diagonal information matrix for the beta parameters.
#' 
#' @param j Category j. 
#' @param p Number of covariates in the model. 
#' @param info Full sparse information matrix for beta.
#' 
#' @return The inverse of the block of the info matrix corresponding with category j.
#' 
#' @export 
invert_block <- function(j, p, info) {
  inds <- ((j - 1)*p + 1):((j - 1)*p + p)
  return(chol2inv(chol(info[inds, inds])))
}
