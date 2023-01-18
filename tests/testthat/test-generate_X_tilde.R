test_that("the XB + z and X_tilde*theta elements match as expected", {
  n <- 10
  J <- 3
  X <- cbind(1, rep(c(0, 1), each = n/2))
  z <- rnorm(n) + 8
  b0 <- rnorm(J)
  b1 <- 1:J
  b <- rbind(b0, b1)
  theta <- generate_theta(generate_B_tilde(b), z)
  X_tilde <- generate_X_tilde(X, J)
  i <- sample(1:n, 1)
  j <- sample(1:J, 1)
  k <- (i - 1)*J + j
  expect_equal(as.numeric(X[i, ] %*% b[, j] + z[i]), (X_tilde %*% theta)[k])
})
