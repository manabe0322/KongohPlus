test_that("calc_weighted_lus", {
  # Condition
  lus_one_al <- c(12, 13, 12, 14)
  count_one_al <- c(99, 400, 1, 300)

  # Run
  wlus <- calc_weighted_lus(lus_one_al, count_one_al)

  # Test
  expect_equal(wlus, 13.25)
})
