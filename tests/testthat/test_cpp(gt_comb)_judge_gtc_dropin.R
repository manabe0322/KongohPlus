test_that("judge_gtc_dropin, pattern 1 (three drop-ins, n_dropin_fltr = 4)", {
  # Condition
  peaks_rm99 <- c(16:19, 21, 22, 23)
  heights_rm99 <- c(209, 2609, 1333, 649, 150, 140, 130)
  gtc <- c(17, 18, 17, 19)
  sr_consider <- c(T, T, T, F)
  n_dropin_fltr <- 4
  h_dropin_fltr <- 300
  bool_gtc_di <- rep(FALSE, length(peaks_rm99) + 1)
  bool_gtc_di[1] <- TRUE

  # Run
  cand_peak <- make_cand_peak(gtc, sr_consider)
  index_dropin <- search_index_dropin(peaks_rm99, cand_peak)
  bool_gtc_di_new <- judge_gtc_dropin(peaks_rm99, heights_rm99, cand_peak, index_dropin, n_dropin_fltr, h_dropin_fltr, bool_gtc_di)

  # Test
  expect_equal(length(bool_gtc_di_new), length(bool_gtc_di))
  expect_equal(bool_gtc_di_new[1], TRUE)
  expect_equal(bool_gtc_di_new[2], FALSE)
  expect_equal(bool_gtc_di_new[3], FALSE)
  expect_equal(bool_gtc_di_new[4], FALSE)
  expect_equal(bool_gtc_di_new[5], FALSE)
  expect_equal(bool_gtc_di_new[6], TRUE)
  expect_equal(bool_gtc_di_new[7], TRUE)
  expect_equal(bool_gtc_di_new[8], TRUE)
})

test_that("judge_gtc_dropin, pattern 2  (three drop-ins, n_dropin_fltr = 2)", {
  # Condition
  peaks_rm99 <- c(16:19, 21, 22, 23)
  heights_rm99 <- c(209, 2609, 1333, 649, 150, 140, 130)
  gtc <- c(17, 18, 17, 19)
  sr_consider <- c(T, T, T, F)
  n_dropin_fltr <- 2
  h_dropin_fltr <- 300
  bool_gtc_di <- rep(FALSE, length(peaks_rm99) + 1)
  bool_gtc_di[1] <- TRUE

  # Run
  cand_peak <- make_cand_peak(gtc, sr_consider)
  index_dropin <- search_index_dropin(peaks_rm99, cand_peak)
  bool_gtc_di_new <- judge_gtc_dropin(peaks_rm99, heights_rm99, cand_peak, index_dropin, n_dropin_fltr, h_dropin_fltr, bool_gtc_di)

  # Test
  expect_equal(length(bool_gtc_di_new), length(bool_gtc_di))
  expect_equal(bool_gtc_di_new[1], FALSE)
  expect_equal(bool_gtc_di_new[2], FALSE)
  expect_equal(bool_gtc_di_new[3], FALSE)
  expect_equal(bool_gtc_di_new[4], FALSE)
  expect_equal(bool_gtc_di_new[5], FALSE)
  expect_equal(bool_gtc_di_new[6], FALSE)
  expect_equal(bool_gtc_di_new[7], FALSE)
  expect_equal(bool_gtc_di_new[8], FALSE)
})

test_that("judge_gtc_dropin, pattern 3 (h_dropin_fltr)", {
  # Condition
  peaks_rm99 <- c(16:19, 21, 22, 23)
  heights_rm99 <- c(209, 2609, 1333, 649, 150, 140, 130)
  gtc <- c(17, 18, 17, 19)
  sr_consider <- c(T, T, T, F)
  n_dropin_fltr <- 4
  h_dropin_fltr <- 135
  bool_gtc_di <- rep(FALSE, length(peaks_rm99) + 1)
  bool_gtc_di[1] <- TRUE

  # Run
  cand_peak <- make_cand_peak(gtc, sr_consider)
  index_dropin <- search_index_dropin(peaks_rm99, cand_peak)
  bool_gtc_di_new <- judge_gtc_dropin(peaks_rm99, heights_rm99, cand_peak, index_dropin, n_dropin_fltr, h_dropin_fltr, bool_gtc_di)

  # Test
  expect_equal(length(bool_gtc_di_new), length(bool_gtc_di))
  expect_equal(bool_gtc_di_new[1], FALSE)
  expect_equal(bool_gtc_di_new[2], FALSE)
  expect_equal(bool_gtc_di_new[3], FALSE)
  expect_equal(bool_gtc_di_new[4], FALSE)
  expect_equal(bool_gtc_di_new[5], FALSE)
  expect_equal(bool_gtc_di_new[6], FALSE)
  expect_equal(bool_gtc_di_new[7], FALSE)
  expect_equal(bool_gtc_di_new[8], FALSE)
})

test_that("judge_gtc_dropin, pattern 4 (h_dropin_fltr, fs)", {
  # Condition
  peaks_rm99 <- c(16:19, 20, 22, 23)
  heights_rm99 <- c(209, 2609, 1333, 649, 150, 140, 130)
  gtc <- c(17, 18, 17, 19)
  sr_consider <- c(T, T, T, F)
  n_dropin_fltr <- 4
  h_dropin_fltr <- 145
  bool_gtc_di <- rep(FALSE, length(peaks_rm99) + 1)
  bool_gtc_di[1] <- TRUE

  # Run
  cand_peak <- make_cand_peak(gtc, sr_consider)
  index_dropin <- search_index_dropin(peaks_rm99, cand_peak)
  bool_gtc_di_new <- judge_gtc_dropin(peaks_rm99, heights_rm99, cand_peak, index_dropin, n_dropin_fltr, h_dropin_fltr, bool_gtc_di)

  # Test
  expect_equal(length(bool_gtc_di_new), length(bool_gtc_di))
  expect_equal(bool_gtc_di_new[1], TRUE)
  expect_equal(bool_gtc_di_new[2], FALSE)
  expect_equal(bool_gtc_di_new[3], FALSE)
  expect_equal(bool_gtc_di_new[4], FALSE)
  expect_equal(bool_gtc_di_new[5], FALSE)
  expect_equal(bool_gtc_di_new[6], FALSE)
  expect_equal(bool_gtc_di_new[7], TRUE)
  expect_equal(bool_gtc_di_new[8], TRUE)
})
