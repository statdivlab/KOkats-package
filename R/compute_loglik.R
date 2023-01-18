# compute log likelihood 
#' @export
compute_loglik <- function(Y, X, B, z) {
  J <- ncol(Y)
  log_means <- X %*% B + matrix(z, ncol = 1) %*% matrix(1, ncol = J, nrow = 1)
  ll <- sum((Y*log_means - exp(log_means)))
  return(ll)
}