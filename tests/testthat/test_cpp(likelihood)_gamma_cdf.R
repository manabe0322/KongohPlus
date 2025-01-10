test_that("gamma_cdf", {
  # Condition
  threshold <- 2
  shape <- 2
  scale <- 2

  # Run
  cdf <- gamma_cdf(threshold, shape, scale)

  # Test
  expect_equal(cdf, pgamma(threshold, shape = shape, scale = scale))
})
