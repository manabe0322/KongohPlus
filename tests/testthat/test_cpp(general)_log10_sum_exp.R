test_that("log10_sum_exp", {
  # Condition
  x <- c(log(0.001), log(0.0001), log(0.00001))

  # Run
  val <- log10_sum_exp(x)

  # Test
  expect_equal(round(val, 6), round(log10(0.001 + 0.0001 + 0.00001), 6))
})

test_that("log10_sum_exp", {
  # Condition
  x <- c(-100, -200, -300)

  # Run
  val <- log10_sum_exp(x)

  # Test
  expect_equal(round(val, 6), round(log10(exp(-100) + exp(-200) + exp(-300)), 6))
})

