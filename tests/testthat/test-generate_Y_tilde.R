test_that("the Y and Y_tilde elements match as expected", {
  Y <- matrix(rpois(10, 100), nrow = 5, ncol = 2)
  Y_tilde <- generate_Y_tilde(Y)
  i <- sample(1:5, 1)
  j <- sample(1:2, 1)
  k <- (i - 1)*ncol(Y) + j
  expect_equal(Y[i, j], Y_tilde[k])
})
