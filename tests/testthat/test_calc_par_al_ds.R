test_that("calc_par_al_ds", {
  # Condition
  ds <- c(15, 15.1, 15.2, 15.3)

  # Run
  al <- calc_par_al_ds(ds)

  # Test
  expect_equal(length(al), 4)
  expect_equal(al[1], 17)
  expect_equal(al[2], 17.1)
  expect_equal(al[3], 17.2)
  expect_equal(al[4], 17.3)
})
