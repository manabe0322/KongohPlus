test_that("calc_log_sum_10, pattern 1", {
  # Condition
  log10_vals <- c(-1, -2, -3)

  # Run
  result <- calc_log_sum_10(log10_vals)

  # Test
  expect_equal(result, log10(0.111))
})

test_that("calc_log_sum_10, pattern 2", {
  # Condition
  log10_vals <- c(-321, -322, -323)

  # Run
  result <- calc_log_sum_10(log10_vals)

  # Test
  expect_equal(result, -320 + log10(0.111))
})

test_that("calc_log_sum_10, pattern 3", {
  # Condition
  log10_vals <- c(-1, -3, -5)

  # Run
  result <- calc_log_sum_10(log10_vals)

  # Test
  expect_equal(result, log10(0.10101))
})

test_that("calc_log_sum_10, pattern 4", {
  # Condition
  log10_vals <- c(-1000, -998, -996)

  # Run
  result <- calc_log_sum_10(log10_vals)

  # Test
  expect_equal(result, -995 + log10(0.10101))
})

test_that("calc_log_sum_10, pattern 5", {
  # Condition
  log10_vals <- c(-1000, -980, -960)

  # Run
  result <- calc_log_sum_10(log10_vals)

  # Test
  expect_equal(result, -960)
})

