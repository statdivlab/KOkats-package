# update B values using Poisson regression
#' @export 
update_Bj_list <- function(theta) {
  j <- theta$j
  Y <- theta$Y
  X <- theta$X
  B <- theta$B
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
