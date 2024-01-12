test_that("make_cand_peak pattern1", {

  # Conditions
  gt_comb <- c(10, 15, 17, 19.2)
  sr_consider <- rep(TRUE, 4)

  # Run
  cand_peak <- make_cand_peak(gt_comb, sr_consider)

  # Test
  expect_equal(setequal(cand_peak, c(8, 9, 9.2, 10, 11, 13, 14, 14.2, 15, 16, 16.2, 17, 17.2, 18, 18.2, 19, 19.2, 20.2)), TRUE)
})

test_that("make_cand_peak pattern2", {

  # Conditions
  gt_comb <- c(10, 15, 17, 19.2)
  sr_consider <- c(TRUE, TRUE, TRUE, FALSE)

  # Run
  cand_peak <- make_cand_peak(gt_comb, sr_consider)

  # Test
  expect_equal(setequal(cand_peak, c(8, 9, 10, 11, 13, 14, 15, 16, 17, 17.2, 18, 18.2, 19.2, 20.2)), TRUE)
})

test_that("make_cand_peak pattern3", {

  # Conditions
  gt_comb <- c(10, 15, 17, 19.2)
  sr_consider <- c(TRUE, TRUE, FALSE, FALSE)

  # Run
  cand_peak <- make_cand_peak(gt_comb, sr_consider)

  # Test
  expect_equal(setequal(cand_peak, c(9, 10, 11, 14, 15, 16, 17, 18, 18.2, 19.2, 20.2)), TRUE)
})

test_that("make_cand_peak pattern4", {

  # Conditions
  gt_comb <- c(10, 15, 17, 19.2)
  sr_consider <- c(TRUE, FALSE, FALSE, FALSE)

  # Run
  cand_peak <- make_cand_peak(gt_comb, sr_consider)

  # Test
  expect_equal(setequal(cand_peak, c(9, 10, 14, 15, 16, 17, 18.2, 19.2)), TRUE)
})

test_that("make_cand_peak pattern5", {

  # Conditions
  gt_comb <- c(10, 15, 17, 19.2)
  sr_consider <- rep(FALSE, 4)

  # Run
  cand_peak <- make_cand_peak(gt_comb, sr_consider)

  # Test
  expect_equal(setequal(cand_peak, c(10, 15, 17, 19.2)), TRUE)
})
