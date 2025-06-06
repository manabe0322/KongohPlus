test_that("softmax", {
  # Condition
  x <- c(1, 2, 3)

  # Run
  val <- softmax(x)

  # Test
  expect_equal(round(val[1], 4), 0.0900)
  expect_equal(round(val[2], 4), 0.2447)
  expect_equal(round(val[3], 4), 0.6652)
})
