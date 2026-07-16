test_that("calc_parent_al_m2s", {
  # Condition
  m2s <- c(15, 15.1, 15.2, 15.3)

  # Run
  al <- calc_parent_al_m2s(m2s)

  # Test
  expect_equal(length(al), 4)
  expect_equal(al[1], 15.2)
  expect_equal(al[2], 15.3)
  expect_equal(al[3], 16)
  expect_equal(al[4], 16.1)
})
