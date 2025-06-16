test_that("select_alternative_hyp, pattern 1", {
  # Condition
  hyp_id_best <- c(0, 1, 0)
  hyp_id_all_alternative <- make_hyp_id_all(3, c("A"), c("A", "B", "C"))
  overall_log10like_all_hyp_alternative <- c(-100, -110, -120)

  # Run
  index_alternative <- select_alternative_hyp(hyp_id_best, hyp_id_all_alternative, overall_log10like_all_hyp_alternative)

  # Test
  expect_equal(index_alternative, 2)
})

test_that("select_alternative_hyp, pattern 2", {
  # Condition
  hyp_id_best <- c(0, 2, 1)
  hyp_id_all_alternative <- make_hyp_id_all(3, c("A", "B"), c("A", "B", "C"))
  overall_log10like_all_hyp_alternative <- c(-100, -110, -120, -130, -140, -150)

  # Run
  index_alternative <- select_alternative_hyp(hyp_id_best, hyp_id_all_alternative, overall_log10like_all_hyp_alternative)

  # Test
  expect_equal(index_alternative, 6)
})

test_that("select_alternative_hyp, pattern 3", {
  # Condition
  hyp_id_best <- c(2, 3, 1)
  hyp_id_all_alternative <- make_hyp_id_all(3, c("A"), c("A", "B", "C"))
  overall_log10like_all_hyp_alternative <- c(-100, -110, -120)

  # Run
  index_alternative <- select_alternative_hyp(hyp_id_best, hyp_id_all_alternative, overall_log10like_all_hyp_alternative)

  # Test
  expect_equal(index_alternative, 3)
})

test_that("select_alternative_hyp, pattern 4", {
  # Condition
  hyp_id_best <- c(3, 0, 1)
  hyp_id_all_alternative <- make_hyp_id_all(3, c("A", "B"), c("A", "B", "C"))
  overall_log10like_all_hyp_alternative <- c(-100, -110, -120, -130, -140, -150)

  # Run
  index_alternative <- select_alternative_hyp(hyp_id_best, hyp_id_all_alternative, overall_log10like_all_hyp_alternative)

  # Test
  expect_equal(index_alternative, 5)
})

test_that("select_alternative_hyp, pattern 5", {
  # Condition
  hyp_id_best <- c(0, 1, 0)
  hyp_id_all_alternative <- make_hyp_id_all(3, c("A", "B"), c("A", "B", "C"))
  overall_log10like_all_hyp_alternative <- c(-100, -110, -120, -130, -140, -150)

  # Run
  index_alternative <- select_alternative_hyp(hyp_id_best, hyp_id_all_alternative, overall_log10like_all_hyp_alternative)

  # Test
  expect_equal(index_alternative, 3)
})
