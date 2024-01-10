test_that("calc_bs", {
  # Condition
  al <- c(15, 15.1, 15.2, 15.3)

  # Run
  bs <- calc_bs(al)

  # Test
  expect_equal(length(bs), 4)
  expect_equal(bs[1], 14)
  expect_equal(bs[2], 14.1)
  expect_equal(bs[3], 14.2)
  expect_equal(bs[4], 14.3)
})
