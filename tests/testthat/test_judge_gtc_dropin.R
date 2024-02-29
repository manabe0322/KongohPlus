test_that("judge_gtc_dropin, condition1", {
  # Condition
  peaks_rm99 <- c(16:19, 21, 22, 23)
  heights_rm99 <- c(209, 2609, 1333, 649, 150, 140, 130)
  gtc <- c(17, 18, 17, 19)
  sr_consider <- c(T, T, T, F)
  n_dropin_fltr <- 4
  h_dropin_fltr <- 300

  # Run
  cand_peak <- make_cand_peak(gtc, sr_consider)
  index_dropin <- search_index_dropin(peaks_rm99, cand_peak)
  adopt_gtc <- judge_gtc_dropin(peaks_rm99, heights_rm99, cand_peak, index_dropin, n_dropin_fltr, h_dropin_fltr)

  # Test
  expect_equal(adopt_gtc, TRUE)
})

test_that("judge_gtc_dropin, condition2", {
  # Condition
  peaks_rm99 <- c(16:19, 21, 22, 23)
  heights_rm99 <- c(209, 2609, 1333, 649, 150, 140, 130)
  gtc <- c(17, 18, 17, 19)
  sr_consider <- c(T, T, T, F)
  n_dropin_fltr <- 2
  h_dropin_fltr <- 300

  # Run
  cand_peak <- make_cand_peak(gtc, sr_consider)
  index_dropin <- search_index_dropin(peaks_rm99, cand_peak)
  adopt_gtc <- judge_gtc_dropin(peaks_rm99, heights_rm99, cand_peak, index_dropin, n_dropin_fltr, h_dropin_fltr)

  # Test
  expect_equal(adopt_gtc, FALSE)
})

test_that("judge_gtc_dropin, condition3", {
  # Condition
  peaks_rm99 <- c(16:19, 21, 22, 23)
  heights_rm99 <- c(209, 2609, 1333, 649, 150, 140, 130)
  gtc <- c(17, 18, 17, 19)
  sr_consider <- c(T, T, T, F)
  n_dropin_fltr <- 4
  h_dropin_fltr <- 135

  # Run
  cand_peak <- make_cand_peak(gtc, sr_consider)
  index_dropin <- search_index_dropin(peaks_rm99, cand_peak)
  adopt_gtc <- judge_gtc_dropin(peaks_rm99, heights_rm99, cand_peak, index_dropin, n_dropin_fltr, h_dropin_fltr)

  # Test
  expect_equal(adopt_gtc, FALSE)
})

test_that("judge_gtc_dropin, condition4", {
  # Condition
  peaks_rm99 <- c(16:19, 20, 22, 23)
  heights_rm99 <- c(209, 2609, 1333, 649, 150, 140, 130)
  gtc <- c(17, 18, 17, 19)
  sr_consider <- c(T, T, T, F)
  n_dropin_fltr <- 4
  h_dropin_fltr <- 145

  # Run
  cand_peak <- make_cand_peak(gtc, sr_consider)
  index_dropin <- search_index_dropin(peaks_rm99, cand_peak)
  adopt_gtc <- judge_gtc_dropin(peaks_rm99, heights_rm99, cand_peak, index_dropin, n_dropin_fltr, h_dropin_fltr)

  # Test
  expect_equal(adopt_gtc, TRUE)
})
