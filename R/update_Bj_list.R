#' Update Bj values 
#' Use Poisson regression to update the Bj block in the block coordinate descent algorithm.
#' 
#' @param theta A list containing all necessary parameters for this function. This list includes \code{j}, the current KO to be updated, \code{Y}, the outcome matrix, \code{X}, the design matrix, \code{B}, the current value of the B matrix, and \code{z}, the current value of the z vector.
#' 
#' @return A vector of new values of B for KO j, the coefficients from fitting a glm.
#' 
#' @export 
update_Bj_list <- function(theta) {
  j <- theta$j
  Y <- theta$Y
  X <- theta$X
  z <- theta$z
  maxit_glm <- theta$maxit_glm
  
  Bj_update <- try(fastglm::fastglm(x = X,
                                    offset = z,
                                    family = "poisson",
                                    control = list(
                                      maxit = maxit_glm),
                                    y = Y[, j])
  )
  
  if(is.null(Bj_update)){
    stop("glm update failed!")
  }
  
  return(Bj_update$coefficients)
}
