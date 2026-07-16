test_that("calc_parent_al_bs", {
  # Condition
  bs <- c(15, 15.1, 15.2, 15.3)

  # Run
  al <- calc_parent_al_bs(bs)

  # Test
  expect_equal(length(al), 4)
  expect_equal(al[1], 16)
  expect_equal(al[2], 16.1)
  expect_equal(al[3], 16.2)
  expect_equal(al[4], 16.3)
})
