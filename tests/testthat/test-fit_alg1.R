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
