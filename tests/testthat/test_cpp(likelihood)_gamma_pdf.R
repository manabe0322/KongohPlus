test_that("gamma_pdf, pattern 1", {
  # Condition
  value <- 2
  shape <- 2
  scale <- 2

  # Run
  pdf <- gamma_pdf(value, shape, scale)

  # Test
  expect_equal(pdf, dgamma(value, shape = shape, scale = scale))
})

test_that("gamma_pdf, pattern 2", {
  # Condition
  value <- 0
  shape <- 2
  scale <- 2

  # Run
  pdf <- gamma_pdf(value, shape, scale)

  # Test
  expect_equal(pdf, dgamma(value, shape = shape, scale = scale))
})

test_that("gamma_pdf, pattern 3", {
  # Condition
  value <- -1
  shape <- 2
  scale <- 2

  # Run
  pdf <- gamma_pdf(value, shape, scale)

  # Test
  expect_equal(pdf, 0)
})

test_that("gamma_pdf, pattern 4", {
  # Condition
  value <- 2
  shape <- -1
  scale <- 2

  # Run
  pdf <- gamma_pdf(value, shape, scale)

  # Test
  expect_equal(pdf, 0)
})

test_that("gamma_pdf, pattern 5", {
  # Condition
  value <- 2
  shape <- 2
  scale <- -1

  # Run
  pdf <- gamma_pdf(value, shape, scale)

  # Test
  expect_equal(pdf, 0)
})

test_that("gamma_pdf, pattern 6", {
  # Condition
  value <- 2
  shape <- Inf
  scale <- 2

  # Run
  pdf <- gamma_pdf(value, shape, scale)

  # Test
  expect_equal(pdf, 0)
})

test_that("gamma_pdf, pattern 6", {
  # Condition
  value <- 2
  shape <- 2
  scale <- Inf

  # Run
  pdf <- gamma_pdf(value, shape, scale)

  # Test
  expect_equal(pdf, 0)
})
