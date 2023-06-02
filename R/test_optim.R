generate_data <- function(n, J, seed = 1) {
  set.seed(seed)
  X <- cbind(1, rep(c(0, 1), each = n/2))
  z <- rnorm(n) + 8
  b0 <- rnorm(J)
  b1 <- 1:J
  b <- rbind(b0, b1)
  Y <- matrix(NA, ncol = J, nrow = n)

  for (i in 1:n) {
   for (j in 1:J) {
     temp_mean <- exp(X[i, , drop = FALSE] %*% b[, j, drop = FALSE] + z[i])
     Y[i,j] <- rpois(1, lambda = temp_mean)
   }
  }
  return(list(X = X, Y = Y))
}

fc_via_optim <- function(Y, X, constraint) {
  J <- ncol(Y)
  p <- ncol(X)
  n <- nrow(X)
  par <- rep(0, p*(J-1) + n - 1)
  if (constraint == "first") {
    optim(par, fn_to_min, control = list(maxit = 500000))
  } else {
    optim(par, fn_to_min_mean, control = list(maxit = 500000))
  }
}

fn_to_min <- function(theta_sm) {
  theta <- c(rep(0, p), theta_sm[1], 0, theta_sm[2:(p*(J - 1) + n - 1)])
  B <- matrix(theta[1:(p*J)], nrow = 2)
  z <- theta[(p*J + 1):length(theta)]
  return(-compute_loglik(Y, X, B, z))
}

fn_to_min_mean <- function(theta_sm) {
  theta <- c(rep(0, p), theta_sm[1], 0, theta_sm[2:(p*(J - 1) + n - 1)])
  B <- matrix(theta[1:(p*J)], nrow = 2)
  B[, 1] <- -rowSums(B[, 2:J])
  z <- theta[(p*J + 1):length(theta)]
  return(-compute_loglik(Y, X, B, z))
}

# to run:
# dat <- generate_data(40, 10) 
# p <- 2
# J <- 10
# n <- 40 
# optim_res <- fc_via_optim(dat$Y, dat$X, "first")
# optim_res_mean <- fc_via_optim(dat$Y, dat$X, "mean")