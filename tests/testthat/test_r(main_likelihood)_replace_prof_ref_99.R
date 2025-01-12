test_that("replace_prof_ref_99", {
  # Condition
  prof_ref <- c(11, 12, 13, 14)
  peaks <- c(11, 12, 13, 99)

  # Run
  prof_ref_99 <- replace_prof_ref_99(prof_ref, peaks)

  # Test
  expect_equal(length(prof_ref_99), 4)
  expect_equal(prof_ref_99[1], 11)
  expect_equal(prof_ref_99[2], 12)
  expect_equal(prof_ref_99[3], 13)
  expect_equal(prof_ref_99[4], 99)
})
