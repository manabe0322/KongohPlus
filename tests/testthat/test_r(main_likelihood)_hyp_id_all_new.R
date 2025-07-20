test_that("select_hyp_id_assumed_contrib", {
  # Condition
  hyp_id_all <- make_hyp_id_all(3, c("A", "B", "C"), c("A", "B", "C"))
  hyp_id_best_assumed_contrib <- c(0, 3, 0)

  # Run
  hyp_id_all_new <- select_hyp_id_assumed_contrib(hyp_id_all, hyp_id_best_assumed_contrib)

  # Test
  expect_equal(nrow(hyp_id_all_new), 2)
  expect_equal(ncol(hyp_id_all_new), 3)
  expect_equal(hyp_id_all_new[1, 1], 1)
  expect_equal(hyp_id_all_new[1, 2], 3)
  expect_equal(hyp_id_all_new[1, 3], 2)
  expect_equal(hyp_id_all_new[2, 1], 2)
  expect_equal(hyp_id_all_new[2, 2], 3)
  expect_equal(hyp_id_all_new[2, 3], 1)
})
