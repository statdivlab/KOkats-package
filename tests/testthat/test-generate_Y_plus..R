test_that("the Y_plus and Y_tilde_plus elements match as expected", {
  n <- 5
  J <- 2
  p <- 2
  Y <- matrix(rpois(n*J, 100), nrow = n, ncol = J)
  Y_tilde <- generate_Y_tilde(Y)
  X <- cbind(rep(1, n), rnorm(n))
  X_tilde <- generate_X_tilde(X, J)
  theta <- rnorm(p*J + n)
  W <- generate_W(X_tilde, theta)
  Y_tilde_plus <- generate_Y_tilde_plus(Y_tilde, X_tilde, W, cores = 2)
  Y_plus <- generate_Y_plus(Y_tilde_plus, J)
  i <- sample(1:n, 1)
  j <- sample(1:J, 1)
  k <- (i - 1)*J + j
  expect_equal(Y_plus[i, j], Y_tilde_plus[k])
})

