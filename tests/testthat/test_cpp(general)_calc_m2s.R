test_that("calc_m2s", {
  # Condition
  al <- c(15, 15.1, 15.2, 15.3)

  # Run
  m2s <- calc_m2s(al)

  # Test
  expect_equal(length(m2s), 4)
  expect_equal(m2s[1], 14.2)
  expect_equal(m2s[2], 14.3)
  expect_equal(m2s[3], 15)
  expect_equal(m2s[4], 15.1)
})
