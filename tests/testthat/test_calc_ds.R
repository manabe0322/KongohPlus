test_that("calc_ds", {
  # Condition
  al <- c(15, 15.1, 15.2, 15.3)

  # Run
  ds <- calc_ds(al)

  # Test
  expect_equal(length(ds), 4)
  expect_equal(ds[1], 13)
  expect_equal(ds[2], 13.1)
  expect_equal(ds[3], 13.2)
  expect_equal(ds[4], 13.3)
})
