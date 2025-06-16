test_that("log_sum_exp, pattern 1", {
  # Condition
  x <- c(log(1), log(2), log(3))

  # Run
  val <- log_sum_exp(x)

  # Test
  expect_equal(round(val, 6), round(log(6), 6))
})

test_that("log_sum_exp, pattern 2", {
  # Condition
  x <- c(log(10^-300), log(10^-301), log(10^-302))

  # Run
  val <- log_sum_exp(x)

  # Test
  expect_equal(round(val, 6), round(log(10^-300) + log(1.11), 6))
})
