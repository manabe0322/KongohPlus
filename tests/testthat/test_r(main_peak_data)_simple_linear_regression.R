test_that("simple_linear_regression", {
  # Condition
  x <- c(13, 14, 15, 16, 17)
  y <- c(120, 124, 128, 132, 136)

  # Run
  param <- simple_linear_regression(x, y)

  # Test
  expect_equal(length(param), 2)
  expect_equal(as.numeric(param[1]), 4)
  expect_equal(as.numeric(param[2]), 68)
})
