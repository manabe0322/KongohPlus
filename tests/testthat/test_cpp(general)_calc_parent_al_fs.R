test_that("calc_parent_al_fs", {
  # Condition
  fs <- c(15, 15.1, 15.2, 15.3)

  # Run
  al <- calc_parent_al_fs(fs)

  # Test
  expect_equal(length(al), 4)
  expect_equal(al[1], 14)
  expect_equal(al[2], 14.1)
  expect_equal(al[3], 14.2)
  expect_equal(al[4], 14.3)
})
