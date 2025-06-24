test_that("gamma_cdf, pattern 1", {
  # Condition
  threshold <- 2
  shape <- 2
  scale <- 2

  # Run
  cdf <- gamma_cdf(threshold, shape, scale)

  # Test
  expect_equal(cdf, pgamma(threshold, shape = shape, scale = scale))
})

test_that("gamma_cdf, pattern 2", {
  # Condition
  threshold <- 0
  shape <- 2
  scale <- 2

  # Run
  cdf <- gamma_cdf(threshold, shape, scale)

  # Test
  expect_equal(cdf, pgamma(threshold, shape = shape, scale = scale))
})

test_that("gamma_cdf, pattern 3", {
  # Condition
  threshold <- -1
  shape <- 2
  scale <- 2

  # Run
  cdf <- gamma_cdf(threshold, shape, scale)

  # Test
  expect_equal(cdf, 0)
})

test_that("gamma_cdf, pattern 4", {
  # Condition
  threshold <- 2
  shape <- -1
  scale <- 2

  # Run
  cdf <- gamma_cdf(threshold, shape, scale)

  # Test
  expect_equal(cdf, 0)
})

test_that("gamma_cdf, pattern 5", {
  # Condition
  threshold <- 2
  shape <- 2
  scale <- -1

  # Run
  cdf <- gamma_cdf(threshold, shape, scale)

  # Test
  expect_equal(cdf, 0)
})

test_that("gamma_cdf, pattern 6", {
  # Condition
  threshold <- 2
  shape <- Inf
  scale <- 2

  # Run
  cdf <- gamma_cdf(threshold, shape, scale)

  # Test
  expect_equal(cdf, 0)
})

test_that("gamma_cdf, pattern 7", {
  # Condition
  threshold <- 2
  shape <- 2
  scale <- Inf

  # Run
  cdf <- gamma_cdf(threshold, shape, scale)

  # Test
  expect_equal(cdf, 0)
})
