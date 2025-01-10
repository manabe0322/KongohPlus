test_that("gamma_pdf", {
  # Condition
  value <- 2
  shape <- 2
  scale <- 2

  # Run
  pdf <- gamma_pdf(value, shape, scale)

  # Test
  expect_equal(pdf, dgamma(value, shape = shape, scale = scale))
})
