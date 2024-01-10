test_that("calc_fs", {
  # Condition
  al <- c(15, 15.1, 15.2, 15.3)

  # Run
  fs <- calc_fs(al)

  # Test
  expect_equal(length(fs), 4)
  expect_equal(fs[1], 16)
  expect_equal(fs[2], 16.1)
  expect_equal(fs[3], 16.2)
  expect_equal(fs[4], 16.3)
})
