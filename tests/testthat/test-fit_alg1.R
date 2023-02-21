test_that("last log likelihood value is highest", {
  X <- cbind(1, rep(c(0, 1), each = 20))
  z <- rnorm(40) + 8
  b0 <- rnorm(10)
  b1 <- 1:10
  b <- rbind(b0, b1)
  Y <- matrix(NA, ncol = 10, nrow = 40)

  for (i in 1:40) {
   for (j in 1:10) {
     temp_mean <- exp(X[i, , drop = FALSE] %*% b[, j, drop = FALSE] + z[i])
     Y[i,j] <- rpois(1, lambda = temp_mean)
   }
  }

  res <- fit_alg1(Y = Y, X = X, constraint_fn = function(x) {median(x)}, ncores = 2)
  expect_true(res$likelihood[length(res$likelihood)] == max(res$likelihood))
})

test_that("fast alg gives same result as slow alg", {
  set.seed(1)
  n <- 10
  J <- 5
  X <- cbind(1, rep(0:1, n/2))
  #X <- cbind(1, rnorm(n))
  z <- rnorm(n) + 8
  b0 <- rnorm(J, 20)
  b1 <- rnorm(J, 2)
  b <- rbind(b0, b1)
  Y <- matrix(NA, ncol = J, nrow = n)

  for (i in 1:n) {
   for (j in 1:J) {
     temp_mean <- exp(X[i, , drop = FALSE] %*% b[, j, drop = FALSE] + z[i])
     Y[i,j] <- rpois(1, lambda = temp_mean)
   }
  }
  Y[1, 1] <- 0

  #res1 <- fit_alg1(Y = Y, X = X, constraint_fn = function(x) {mean(x)}, ncores = 2)
  res1 <- fit_alg1(Y = Y, X = X, constraint_fn = function(x) {mean(x)}, ncores = 2)
  res2 <- fit_alg1_alt(Y = Y, X = X, constraint_fn = function(x) {mean(x)}, maxit = 1000, ncores = 2)
  ll1 <- compute_loglik(Y, X, res1$final_B, res1$final_z)
  ll2 <- compute_loglik(Y, X, res2$final_B, res2$final_z)
  ll1 - ll2
  (res1$final_B - res2$final_B)/res1$final_B
  (res1$final_z - res2$final_z)/res1$final_z
  ll1 - max(res2$likelihood)
  res2_length <- length(res2$likelihood)
  (res2$likelihood[res2_length] - res2$likelihood[res2_length - 1])/res2$likelihood[res2_length]
  (res2$B[,,res2_length] - res2$B[,,res2_length - 1])/res2$B[,,res2_length]
})
