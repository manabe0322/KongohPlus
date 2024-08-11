test_that("judge_gtc_sr, pattern 1 (general)", {
  # Condition
  gtc <- c(10, 20, 30, 40)
  peaks <- c(8, 9, 9.2, 10, 11, 18, 19, 20, 21, 28, 29, 30, 31, 38, 39, 40, 41, 99)
  peaks_rm99_rmdi <- c(8, 9, 9.2, 10, 11, 18, 19, 20, 21, 28, 29, 30, 31, 38, 39, 40, 41)
  heights_rm99_rmdi <- c(100, 1000, 200, 10000, 500, 100, 1000, 10000, 500, 100, 1000, 10000, 500, 100, 1000, 10000, 500)
  bsrs <- c(0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0)
  fsrs <- c(0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0)
  dsrs <- c(0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0)
  m2srs <- c(0, 0, 0, 0.02, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  bsr_fltr <- 0.5
  fsr_fltr <- 0.5
  dsr_fltr <- 0.5
  m2sr_fltr <- 0.5
  n_dropin_fltr <- 2
  h_dropin_fltr <- 50
  n_dropin <- 0
  bool_gtc_di <- rep(FALSE, length(peaks))
  bool_gtc_di[1] <- TRUE

  # Run
  bool_gtc_di_new <-  judge_gtc_sr(gtc, peaks, peaks_rm99_rmdi, heights_rm99_rmdi,
                                   bsrs, fsrs, dsrs, m2srs,
                                   bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr,
                                   n_dropin_fltr, h_dropin_fltr, n_dropin, bool_gtc_di)

  # Test
  expect_equal(length(bool_gtc_di_new), length(bool_gtc_di))
  expect_equal(bool_gtc_di_new[1], TRUE)
  expect_equal(any(bool_gtc_di_new[-1]), FALSE)
})

test_that("judge_gtc_sr, pattern 2 (bsrs)", {
  # Condition
  gtc <- c(10, 20, 30, 40)
  peaks <- c(8, 9, 9.2, 10, 11, 18, 19, 20, 21, 28, 29, 30, 31, 38, 39, 40, 41, 99)
  peaks_rm99_rmdi <- c(8, 9, 9.2, 10, 11, 18, 19, 20, 21, 28, 29, 30, 31, 38, 39, 40, 41)
  heights_rm99_rmdi <- c(100, 1000, 200, 10000, 500, 100, 1000, 10000, 500, 100, 1000, 10000, 500, 100, 1000, 10000, 500)
  bsrs <- c(0, 0, 0, 0, 0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0)
  fsrs <- c(0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0)
  dsrs <- c(0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0)
  m2srs <- c(0, 0, 0, 0.02, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  bsr_fltr <- 0.5
  fsr_fltr <- 0.5
  dsr_fltr <- 0.5
  m2sr_fltr <- 0.5
  n_dropin_fltr <- 2
  h_dropin_fltr <- 50
  n_dropin <- 0
  bool_gtc_di <- rep(FALSE, length(peaks))
  bool_gtc_di[1] <- TRUE

  # Run
  bool_gtc_di_new <-  judge_gtc_sr(gtc, peaks, peaks_rm99_rmdi, heights_rm99_rmdi,
                                   bsrs, fsrs, dsrs, m2srs,
                                   bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr,
                                   n_dropin_fltr, h_dropin_fltr, n_dropin, bool_gtc_di)

  # Test
  expect_equal(length(bool_gtc_di_new), length(bool_gtc_di))
  expect_equal(bool_gtc_di_new[1], FALSE)
  expect_equal(any(bool_gtc_di_new[-1]), FALSE)
})

test_that("judge_gtc_sr, pattern 3 (fsrs)", {
  # Condition
  gtc <- c(10, 20, 30, 40)
  peaks <- c(8, 9, 9.2, 10, 11, 18, 19, 20, 21, 28, 29, 30, 31, 38, 39, 40, 41, 99)
  peaks_rm99_rmdi <- c(8, 9, 9.2, 10, 11, 18, 19, 20, 21, 28, 29, 30, 31, 38, 39, 40, 41)
  heights_rm99_rmdi <- c(100, 1000, 200, 10000, 500, 100, 1000, 10000, 500, 100, 1000, 10000, 500, 100, 1000, 10000, 500)
  bsrs <- c(0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0)
  fsrs <- c(0, 0, 0, 0, 0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0)
  dsrs <- c(0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0)
  m2srs <- c(0, 0, 0, 0.02, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  bsr_fltr <- 0.5
  fsr_fltr <- 0.5
  dsr_fltr <- 0.5
  m2sr_fltr <- 0.5
  n_dropin_fltr <- 2
  h_dropin_fltr <- 50
  n_dropin <- 0
  bool_gtc_di <- rep(FALSE, length(peaks))
  bool_gtc_di[1] <- TRUE

  # Run
  bool_gtc_di_new <-  judge_gtc_sr(gtc, peaks, peaks_rm99_rmdi, heights_rm99_rmdi,
                                   bsrs, fsrs, dsrs, m2srs,
                                   bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr,
                                   n_dropin_fltr, h_dropin_fltr, n_dropin, bool_gtc_di)

  # Test
  expect_equal(length(bool_gtc_di_new), length(bool_gtc_di))
  expect_equal(bool_gtc_di_new[1], FALSE)
  expect_equal(any(bool_gtc_di_new[-1]), FALSE)
})

test_that("judge_gtc_sr, pattern 4 (dsrs)", {
  # Condition
  gtc <- c(10, 20, 30, 40)
  peaks <- c(8, 9, 9.2, 10, 11, 18, 19, 20, 21, 28, 29, 30, 31, 38, 39, 40, 41, 99)
  peaks_rm99_rmdi <- c(8, 9, 9.2, 10, 11, 18, 19, 20, 21, 28, 29, 30, 31, 38, 39, 40, 41)
  heights_rm99_rmdi <- c(100, 1000, 200, 10000, 500, 100, 1000, 10000, 500, 100, 1000, 10000, 500, 100, 1000, 10000, 500)
  bsrs <- c(0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0)
  fsrs <- c(0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0)
  dsrs <- c(0, 0, 0, 0, 0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0)
  m2srs <- c(0, 0, 0, 0.02, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  bsr_fltr <- 0.5
  fsr_fltr <- 0.5
  dsr_fltr <- 0.5
  m2sr_fltr <- 0.5
  n_dropin_fltr <- 2
  h_dropin_fltr <- 50
  n_dropin <- 0
  bool_gtc_di <- rep(FALSE, length(peaks))
  bool_gtc_di[1] <- TRUE

  # Run
  bool_gtc_di_new <-  judge_gtc_sr(gtc, peaks, peaks_rm99_rmdi, heights_rm99_rmdi,
                                   bsrs, fsrs, dsrs, m2srs,
                                   bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr,
                                   n_dropin_fltr, h_dropin_fltr, n_dropin, bool_gtc_di)

  # Test
  expect_equal(length(bool_gtc_di_new), length(bool_gtc_di))
  expect_equal(bool_gtc_di_new[1], FALSE)
  expect_equal(any(bool_gtc_di_new[-1]), FALSE)
})

test_that("judge_gtc_sr, pattern 5 (m2srs)", {
  # Condition
  gtc <- c(10, 20, 30, 40)
  peaks <- c(8, 9, 9.2, 10, 11, 18, 19, 20, 21, 28, 29, 30, 31, 38, 39, 40, 41, 99)
  peaks_rm99_rmdi <- c(8, 9, 9.2, 10, 11, 18, 19, 20, 21, 28, 29, 30, 31, 38, 39, 40, 41)
  heights_rm99_rmdi <- c(100, 1000, 200, 10000, 500, 100, 1000, 10000, 500, 100, 1000, 10000, 500, 100, 1000, 10000, 500)
  bsrs <- c(0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0)
  fsrs <- c(0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0)
  dsrs <- c(0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0)
  m2srs <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  bsr_fltr <- 0.5
  fsr_fltr <- 0.5
  dsr_fltr <- 0.5
  m2sr_fltr <- 0.5
  n_dropin_fltr <- 2
  h_dropin_fltr <- 50
  n_dropin <- 0
  bool_gtc_di <- rep(FALSE, length(peaks))
  bool_gtc_di[1] <- TRUE

  # Run
  bool_gtc_di_new <-  judge_gtc_sr(gtc, peaks, peaks_rm99_rmdi, heights_rm99_rmdi,
                                   bsrs, fsrs, dsrs, m2srs,
                                   bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr,
                                   n_dropin_fltr, h_dropin_fltr, n_dropin, bool_gtc_di)

  # Test
  expect_equal(length(bool_gtc_di_new), length(bool_gtc_di))
  expect_equal(bool_gtc_di_new[1], FALSE)
  expect_equal(any(bool_gtc_di_new[-1]), FALSE)
})

test_that("judge_gtc_sr, pattern 6 (bsr_fltr)", {
  # Condition
  gtc <- c(10, 20, 30, 40)
  peaks <- c(8, 9, 9.2, 10, 11, 18, 19, 20, 21, 28, 29, 30, 31, 38, 39, 40, 41, 99)
  peaks_rm99_rmdi <- c(8, 9, 9.2, 10, 11, 18, 19, 20, 21, 28, 29, 30, 31, 38, 39, 40, 41)
  heights_rm99_rmdi <- c(100, 1000, 200, 10000, 500, 100, 1000, 10000, 500, 100, 1000, 10000, 500, 100, 1000, 10000, 500)
  bsrs <- c(0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0)
  fsrs <- c(0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0)
  dsrs <- c(0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0)
  m2srs <- c(0, 0, 0, 0.02, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  bsr_fltr <- 0.05
  fsr_fltr <- 0.5
  dsr_fltr <- 0.5
  m2sr_fltr <- 0.5
  n_dropin_fltr <- 2
  h_dropin_fltr <- 50
  n_dropin <- 0
  bool_gtc_di <- rep(FALSE, length(peaks))
  bool_gtc_di[1] <- TRUE

  # Run
  bool_gtc_di_new <-  judge_gtc_sr(gtc, peaks, peaks_rm99_rmdi, heights_rm99_rmdi,
                                   bsrs, fsrs, dsrs, m2srs,
                                   bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr,
                                   n_dropin_fltr, h_dropin_fltr, n_dropin, bool_gtc_di)

  # Test
  expect_equal(length(bool_gtc_di_new), length(bool_gtc_di))
  expect_equal(bool_gtc_di_new[1], FALSE)
  expect_equal(any(bool_gtc_di_new[-1]), FALSE)
})

test_that("judge_gtc_sr, pattern 7 (fsr_fltr)", {
  # Condition
  gtc <- c(10, 20, 30, 40)
  peaks <- c(8, 9, 9.2, 10, 11, 18, 19, 20, 21, 28, 29, 30, 31, 38, 39, 40, 41, 99)
  peaks_rm99_rmdi <- c(8, 9, 9.2, 10, 11, 18, 19, 20, 21, 28, 29, 30, 31, 38, 39, 40, 41)
  heights_rm99_rmdi <- c(100, 1000, 200, 10000, 500, 100, 1000, 10000, 500, 100, 1000, 10000, 500, 100, 1000, 10000, 500)
  bsrs <- c(0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0)
  fsrs <- c(0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0)
  dsrs <- c(0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0)
  m2srs <- c(0, 0, 0, 0.02, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  bsr_fltr <- 0.5
  fsr_fltr <- 0.025
  dsr_fltr <- 0.5
  m2sr_fltr <- 0.5
  n_dropin_fltr <- 2
  h_dropin_fltr <- 50
  n_dropin <- 0
  bool_gtc_di <- rep(FALSE, length(peaks))
  bool_gtc_di[1] <- TRUE

  # Run
  bool_gtc_di_new <-  judge_gtc_sr(gtc, peaks, peaks_rm99_rmdi, heights_rm99_rmdi,
                                   bsrs, fsrs, dsrs, m2srs,
                                   bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr,
                                   n_dropin_fltr, h_dropin_fltr, n_dropin, bool_gtc_di)

  # Test
  expect_equal(length(bool_gtc_di_new), length(bool_gtc_di))
  expect_equal(bool_gtc_di_new[1], FALSE)
  expect_equal(any(bool_gtc_di_new[-1]), FALSE)
})

test_that("judge_gtc_sr, pattern 8 (dsr_fltr)", {
  # Condition
  gtc <- c(10, 20, 30, 40)
  peaks <- c(8, 9, 9.2, 10, 11, 18, 19, 20, 21, 28, 29, 30, 31, 38, 39, 40, 41, 99)
  peaks_rm99_rmdi <- c(8, 9, 9.2, 10, 11, 18, 19, 20, 21, 28, 29, 30, 31, 38, 39, 40, 41)
  heights_rm99_rmdi <- c(100, 1000, 200, 10000, 500, 100, 1000, 10000, 500, 100, 1000, 10000, 500, 100, 1000, 10000, 500)
  bsrs <- c(0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0)
  fsrs <- c(0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0)
  dsrs <- c(0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0)
  m2srs <- c(0, 0, 0, 0.02, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  bsr_fltr <- 0.5
  fsr_fltr <- 0.5
  dsr_fltr <- 0.005
  m2sr_fltr <- 0.5
  n_dropin_fltr <- 2
  h_dropin_fltr <- 50
  n_dropin <- 0
  bool_gtc_di <- rep(FALSE, length(peaks))
  bool_gtc_di[1] <- TRUE

  # Run
  bool_gtc_di_new <-  judge_gtc_sr(gtc, peaks, peaks_rm99_rmdi, heights_rm99_rmdi,
                                   bsrs, fsrs, dsrs, m2srs,
                                   bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr,
                                   n_dropin_fltr, h_dropin_fltr, n_dropin, bool_gtc_di)

  # Test
  expect_equal(length(bool_gtc_di_new), length(bool_gtc_di))
  expect_equal(bool_gtc_di_new[1], FALSE)
  expect_equal(any(bool_gtc_di_new[-1]), FALSE)
})

test_that("judge_gtc_sr, pattern 9 (m2sr_fltr)", {
  # Condition
  gtc <- c(10, 20, 30, 40)
  peaks <- c(8, 9, 9.2, 10, 11, 18, 19, 20, 21, 28, 29, 30, 31, 38, 39, 40, 41, 99)
  peaks_rm99_rmdi <- c(8, 9, 9.2, 10, 11, 18, 19, 20, 21, 28, 29, 30, 31, 38, 39, 40, 41)
  heights_rm99_rmdi <- c(100, 1000, 200, 10000, 500, 100, 1000, 10000, 500, 100, 1000, 10000, 500, 100, 1000, 10000, 500)
  bsrs <- c(0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0)
  fsrs <- c(0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0)
  dsrs <- c(0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0)
  m2srs <- c(0, 0, 0, 0.02, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  bsr_fltr <- 0.5
  fsr_fltr <- 0.5
  dsr_fltr <- 0.5
  m2sr_fltr <- 0.01
  n_dropin_fltr <- 2
  h_dropin_fltr <- 50
  n_dropin <- 0
  bool_gtc_di <- rep(FALSE, length(peaks))
  bool_gtc_di[1] <- TRUE

  # Run
  bool_gtc_di_new <-  judge_gtc_sr(gtc, peaks, peaks_rm99_rmdi, heights_rm99_rmdi,
                                   bsrs, fsrs, dsrs, m2srs,
                                   bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr,
                                   n_dropin_fltr, h_dropin_fltr, n_dropin, bool_gtc_di)

  # Test
  expect_equal(length(bool_gtc_di_new), length(bool_gtc_di))
  expect_equal(bool_gtc_di_new[1], FALSE)
  expect_equal(any(bool_gtc_di_new[-1]), FALSE)
})

test_that("judge_gtc_sr, pattern 10 (bsrs dropin)", {
  # Condition
  gtc <- c(10, 20, 30, 40)
  peaks <- c(8, 9, 9.2, 10, 11, 18, 19, 20, 21, 28, 29, 30, 31, 38, 39, 40, 41, 99)
  peaks_rm99_rmdi <- c(8, 9, 9.2, 10, 11, 18, 19, 20, 21, 28, 29, 30, 31, 38, 39, 40, 41)
  heights_rm99_rmdi <- c(100, 1000, 200, 10000, 500, 100, 1000, 10000, 500, 100, 1000, 10000, 500, 100, 1000, 10000, 500)
  bsrs <- c(0, 0, 0, 0, 0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0)
  fsrs <- c(0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0)
  dsrs <- c(0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0)
  m2srs <- c(0, 0, 0, 0.02, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  bsr_fltr <- 0.5
  fsr_fltr <- 0.5
  dsr_fltr <- 0.5
  m2sr_fltr <- 0.5
  n_dropin_fltr <- 2
  h_dropin_fltr <- 1000
  n_dropin <- 0
  bool_gtc_di <- rep(FALSE, length(peaks))
  bool_gtc_di[1] <- TRUE

  # Run
  bool_gtc_di_new <-  judge_gtc_sr(gtc, peaks, peaks_rm99_rmdi, heights_rm99_rmdi,
                                   bsrs, fsrs, dsrs, m2srs,
                                   bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr,
                                   n_dropin_fltr, h_dropin_fltr, n_dropin, bool_gtc_di)

  # Test
  expect_equal(length(bool_gtc_di_new), length(bool_gtc_di))
  expect_equal(bool_gtc_di_new[1], TRUE)
  expect_equal(bool_gtc_di_new[2], FALSE)
  expect_equal(bool_gtc_di_new[3], TRUE)
  expect_equal(any(bool_gtc_di_new[4:length(bool_gtc_di_new)]), FALSE)
})

test_that("judge_gtc_sr, pattern 11 (fsrs dropin)", {
  # Condition
  gtc <- c(10, 20, 30, 40)
  peaks <- c(8, 9, 9.2, 10, 11, 18, 19, 20, 21, 28, 29, 30, 31, 38, 39, 40, 41, 99)
  peaks_rm99_rmdi <- c(8, 9, 9.2, 10, 11, 18, 19, 20, 21, 28, 29, 30, 31, 38, 39, 40, 41)
  heights_rm99_rmdi <- c(100, 1000, 200, 10000, 500, 100, 1000, 10000, 500, 100, 1000, 10000, 500, 100, 1000, 10000, 500)
  bsrs <- c(0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0)
  fsrs <- c(0, 0, 0, 0, 0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0)
  dsrs <- c(0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0)
  m2srs <- c(0, 0, 0, 0.02, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  bsr_fltr <- 0.5
  fsr_fltr <- 0.5
  dsr_fltr <- 0.5
  m2sr_fltr <- 0.5
  n_dropin_fltr <- 2
  h_dropin_fltr <- 1000
  n_dropin <- 0
  bool_gtc_di <- rep(FALSE, length(peaks))
  bool_gtc_di[1] <- TRUE

  # Run
  bool_gtc_di_new <-  judge_gtc_sr(gtc, peaks, peaks_rm99_rmdi, heights_rm99_rmdi,
                                   bsrs, fsrs, dsrs, m2srs,
                                   bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr,
                                   n_dropin_fltr, h_dropin_fltr, n_dropin, bool_gtc_di)

  # Test
  expect_equal(length(bool_gtc_di_new), length(bool_gtc_di))
  expect_equal(bool_gtc_di_new[1], TRUE)
  expect_equal(bool_gtc_di_new[2], FALSE)
  expect_equal(bool_gtc_di_new[3], FALSE)
  expect_equal(bool_gtc_di_new[4], FALSE)
  expect_equal(bool_gtc_di_new[5], FALSE)
  expect_equal(bool_gtc_di_new[6], TRUE)
  expect_equal(any(bool_gtc_di_new[7:length(bool_gtc_di_new)]), FALSE)
})

test_that("judge_gtc_sr, pattern 12 (dsrs dropin)", {
  # Condition
  gtc <- c(10, 20, 30, 40)
  peaks <- c(8, 9, 9.2, 10, 11, 18, 19, 20, 21, 28, 29, 30, 31, 38, 39, 40, 41, 99)
  peaks_rm99_rmdi <- c(8, 9, 9.2, 10, 11, 18, 19, 20, 21, 28, 29, 30, 31, 38, 39, 40, 41)
  heights_rm99_rmdi <- c(100, 1000, 200, 10000, 500, 100, 1000, 10000, 500, 100, 1000, 10000, 500, 100, 1000, 10000, 500)
  bsrs <- c(0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0)
  fsrs <- c(0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0)
  dsrs <- c(0, 0, 0, 0, 0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0)
  m2srs <- c(0, 0, 0, 0.02, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  bsr_fltr <- 0.5
  fsr_fltr <- 0.5
  dsr_fltr <- 0.5
  m2sr_fltr <- 0.5
  n_dropin_fltr <- 2
  h_dropin_fltr <- 1000
  n_dropin <- 0
  bool_gtc_di <- rep(FALSE, length(peaks))
  bool_gtc_di[1] <- TRUE

  # Run
  bool_gtc_di_new <-  judge_gtc_sr(gtc, peaks, peaks_rm99_rmdi, heights_rm99_rmdi,
                                   bsrs, fsrs, dsrs, m2srs,
                                   bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr,
                                   n_dropin_fltr, h_dropin_fltr, n_dropin, bool_gtc_di)

  # Test
  expect_equal(length(bool_gtc_di_new), length(bool_gtc_di))
  expect_equal(bool_gtc_di_new[1], TRUE)
  expect_equal(bool_gtc_di_new[2], TRUE)
  expect_equal(any(bool_gtc_di_new[3:length(bool_gtc_di_new)]), FALSE)
})

test_that("judge_gtc_sr, pattern 13 (m2srs dropin)", {
  # Condition
  gtc <- c(10, 20, 30, 40)
  peaks <- c(8, 9, 9.2, 10, 11, 18, 19, 20, 21, 28, 29, 30, 31, 38, 39, 40, 41, 99)
  peaks_rm99_rmdi <- c(8, 9, 9.2, 10, 11, 18, 19, 20, 21, 28, 29, 30, 31, 38, 39, 40, 41)
  heights_rm99_rmdi <- c(100, 1000, 200, 10000, 500, 100, 1000, 10000, 500, 100, 1000, 10000, 500, 100, 1000, 10000, 500)
  bsrs <- c(0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0)
  fsrs <- c(0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0)
  dsrs <- c(0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0)
  m2srs <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  bsr_fltr <- 0.5
  fsr_fltr <- 0.5
  dsr_fltr <- 0.5
  m2sr_fltr <- 0.5
  n_dropin_fltr <- 2
  h_dropin_fltr <- 1000
  n_dropin <- 0
  bool_gtc_di <- rep(FALSE, length(peaks))
  bool_gtc_di[1] <- TRUE

  # Run
  bool_gtc_di_new <-  judge_gtc_sr(gtc, peaks, peaks_rm99_rmdi, heights_rm99_rmdi,
                                   bsrs, fsrs, dsrs, m2srs,
                                   bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr,
                                   n_dropin_fltr, h_dropin_fltr, n_dropin, bool_gtc_di)

  # Test
  expect_equal(length(bool_gtc_di_new), length(bool_gtc_di))
  expect_equal(bool_gtc_di_new[1], TRUE)
  expect_equal(bool_gtc_di_new[2], FALSE)
  expect_equal(bool_gtc_di_new[3], FALSE)
  expect_equal(bool_gtc_di_new[4], TRUE)
  expect_equal(any(bool_gtc_di_new[5:length(bool_gtc_di_new)]), FALSE)
})

test_that("judge_gtc_sr, pattern 14 (bsr_fltr dropin)", {
  # Condition
  gtc <- c(10, 20, 30, 40)
  peaks <- c(8, 9, 9.2, 10, 11, 18, 19, 20, 21, 28, 29, 30, 31, 38, 39, 40, 41, 99)
  peaks_rm99_rmdi <- c(8, 9, 9.2, 10, 11, 18, 19, 20, 21, 28, 29, 30, 31, 38, 39, 40, 41)
  heights_rm99_rmdi <- c(100, 6000, 200, 10000, 500, 100, 1000, 10000, 500, 100, 1000, 10000, 500, 100, 1000, 10000, 500)
  bsrs <- c(0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0)
  fsrs <- c(0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0)
  dsrs <- c(0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0)
  m2srs <- c(0, 0, 0, 0.02, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  bsr_fltr <- 0.5
  fsr_fltr <- 0.5
  dsr_fltr <- 0.5
  m2sr_fltr <- 0.5
  n_dropin_fltr <- 2
  h_dropin_fltr <- 10000
  n_dropin <- 0
  bool_gtc_di <- rep(FALSE, length(peaks))
  bool_gtc_di[1] <- TRUE

  # Run
  bool_gtc_di_new <-  judge_gtc_sr(gtc, peaks, peaks_rm99_rmdi, heights_rm99_rmdi,
                                   bsrs, fsrs, dsrs, m2srs,
                                   bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr,
                                   n_dropin_fltr, h_dropin_fltr, n_dropin, bool_gtc_di)

  # Test
  expect_equal(length(bool_gtc_di_new), length(bool_gtc_di))
  expect_equal(bool_gtc_di_new[1], TRUE)
  expect_equal(bool_gtc_di_new[2], FALSE)
  expect_equal(bool_gtc_di_new[3], TRUE)
  expect_equal(any(bool_gtc_di_new[4:length(bool_gtc_di_new)]), FALSE)
})

test_that("judge_gtc_sr, pattern 15 (fsr_fltr dropin)", {
  # Condition
  gtc <- c(10, 20, 30, 40)
  peaks <- c(8, 9, 9.2, 10, 11, 18, 19, 20, 21, 28, 29, 30, 31, 38, 39, 40, 41, 99)
  peaks_rm99_rmdi <- c(8, 9, 9.2, 10, 11, 18, 19, 20, 21, 28, 29, 30, 31, 38, 39, 40, 41)
  heights_rm99_rmdi <- c(100, 1000, 200, 10000, 6000, 100, 1000, 10000, 500, 100, 1000, 10000, 500, 100, 1000, 10000, 500)
  bsrs <- c(0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0)
  fsrs <- c(0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0)
  dsrs <- c(0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0)
  m2srs <- c(0, 0, 0, 0.02, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  bsr_fltr <- 0.5
  fsr_fltr <- 0.5
  dsr_fltr <- 0.5
  m2sr_fltr <- 0.5
  n_dropin_fltr <- 2
  h_dropin_fltr <- 10000
  n_dropin <- 0
  bool_gtc_di <- rep(FALSE, length(peaks))
  bool_gtc_di[1] <- TRUE

  # Run
  bool_gtc_di_new <-  judge_gtc_sr(gtc, peaks, peaks_rm99_rmdi, heights_rm99_rmdi,
                                   bsrs, fsrs, dsrs, m2srs,
                                   bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr,
                                   n_dropin_fltr, h_dropin_fltr, n_dropin, bool_gtc_di)

  # Test
  expect_equal(length(bool_gtc_di_new), length(bool_gtc_di))
  expect_equal(bool_gtc_di_new[1], TRUE)
  expect_equal(bool_gtc_di_new[2], FALSE)
  expect_equal(bool_gtc_di_new[3], FALSE)
  expect_equal(bool_gtc_di_new[4], FALSE)
  expect_equal(bool_gtc_di_new[5], FALSE)
  expect_equal(bool_gtc_di_new[6], TRUE)
  expect_equal(any(bool_gtc_di_new[7:length(bool_gtc_di_new)]), FALSE)
})

test_that("judge_gtc_sr, pattern 16 (dsr_fltr dropin)", {
  # Condition
  gtc <- c(10, 20, 30, 40)
  peaks <- c(8, 9, 9.2, 10, 11, 18, 19, 20, 21, 28, 29, 30, 31, 38, 39, 40, 41, 99)
  peaks_rm99_rmdi <- c(8, 9, 9.2, 10, 11, 18, 19, 20, 21, 28, 29, 30, 31, 38, 39, 40, 41)
  heights_rm99_rmdi <- c(6000, 1000, 200, 10000, 500, 100, 1000, 10000, 500, 100, 1000, 10000, 500, 100, 1000, 10000, 500)
  bsrs <- c(0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0)
  fsrs <- c(0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0)
  dsrs <- c(0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0)
  m2srs <- c(0, 0, 0, 0.02, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  bsr_fltr <- 0.5
  fsr_fltr <- 0.5
  dsr_fltr <- 0.5
  m2sr_fltr <- 0.5
  n_dropin_fltr <- 2
  h_dropin_fltr <- 10000
  n_dropin <- 0
  bool_gtc_di <- rep(FALSE, length(peaks))
  bool_gtc_di[1] <- TRUE

  # Run
  bool_gtc_di_new <-  judge_gtc_sr(gtc, peaks, peaks_rm99_rmdi, heights_rm99_rmdi,
                                   bsrs, fsrs, dsrs, m2srs,
                                   bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr,
                                   n_dropin_fltr, h_dropin_fltr, n_dropin, bool_gtc_di)

  # Test
  expect_equal(length(bool_gtc_di_new), length(bool_gtc_di))
  expect_equal(bool_gtc_di_new[1], TRUE)
  expect_equal(bool_gtc_di_new[2], TRUE)
  expect_equal(any(bool_gtc_di_new[3:length(bool_gtc_di_new)]), FALSE)
})

test_that("judge_gtc_sr, pattern 17 (m2sr_fltr dropin)", {
  # Condition
  gtc <- c(10, 20, 30, 40)
  peaks <- c(8, 9, 9.2, 10, 11, 18, 19, 20, 21, 28, 29, 30, 31, 38, 39, 40, 41, 99)
  peaks_rm99_rmdi <- c(8, 9, 9.2, 10, 11, 18, 19, 20, 21, 28, 29, 30, 31, 38, 39, 40, 41)
  heights_rm99_rmdi <- c(100, 1000, 6000, 10000, 500, 100, 1000, 10000, 500, 100, 1000, 10000, 500, 100, 1000, 10000, 500)
  bsrs <- c(0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0)
  fsrs <- c(0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0)
  dsrs <- c(0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0)
  m2srs <- c(0, 0, 0, 0.02, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  bsr_fltr <- 0.5
  fsr_fltr <- 0.5
  dsr_fltr <- 0.5
  m2sr_fltr <- 0.5
  n_dropin_fltr <- 2
  h_dropin_fltr <- 10000
  n_dropin <- 0
  bool_gtc_di <- rep(FALSE, length(peaks))
  bool_gtc_di[1] <- TRUE

  # Run
  bool_gtc_di_new <-  judge_gtc_sr(gtc, peaks, peaks_rm99_rmdi, heights_rm99_rmdi,
                                   bsrs, fsrs, dsrs, m2srs,
                                   bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr,
                                   n_dropin_fltr, h_dropin_fltr, n_dropin, bool_gtc_di)

  # Test
  expect_equal(length(bool_gtc_di_new), length(bool_gtc_di))
  expect_equal(bool_gtc_di_new[1], TRUE)
  expect_equal(bool_gtc_di_new[2], FALSE)
  expect_equal(bool_gtc_di_new[3], FALSE)
  expect_equal(bool_gtc_di_new[4], TRUE)
  expect_equal(any(bool_gtc_di_new[5:length(bool_gtc_di_new)]), FALSE)
})

test_that("judge_gtc_sr, pattern 18 (bsrs n_dropin = 1)", {
  # Condition
  gtc <- c(10, 20, 30, 40)
  peaks <- c(8, 9, 9.2, 10, 11, 18, 19, 20, 21, 28, 29, 30, 31, 38, 39, 40, 41, 99)
  peaks_rm99_rmdi <- c(8, 9, 9.2, 10, 11, 18, 19, 20, 21, 28, 29, 30, 31, 38, 39, 40, 41)
  heights_rm99_rmdi <- c(100, 1000, 200, 10000, 500, 100, 1000, 10000, 500, 100, 1000, 10000, 500, 100, 1000, 10000, 500)
  bsrs <- c(0, 0, 0, 0, 0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0)
  fsrs <- c(0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0)
  dsrs <- c(0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0)
  m2srs <- c(0, 0, 0, 0.02, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  bsr_fltr <- 0.5
  fsr_fltr <- 0.5
  dsr_fltr <- 0.5
  m2sr_fltr <- 0.5
  n_dropin_fltr <- 2
  h_dropin_fltr <- 10000
  n_dropin <- 1
  bool_gtc_di <- rep(FALSE, length(peaks))
  bool_gtc_di[1] <- TRUE

  # Run
  bool_gtc_di_new <-  judge_gtc_sr(gtc, peaks, peaks_rm99_rmdi, heights_rm99_rmdi,
                                   bsrs, fsrs, dsrs, m2srs,
                                   bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr,
                                   n_dropin_fltr, h_dropin_fltr, n_dropin, bool_gtc_di)

  # Test
  expect_equal(length(bool_gtc_di_new), length(bool_gtc_di))
  expect_equal(bool_gtc_di_new[1], TRUE)
  expect_equal(bool_gtc_di_new[2], FALSE)
  expect_equal(bool_gtc_di_new[3], TRUE)
  expect_equal(bool_gtc_di_new[4], FALSE)
  expect_equal(any(bool_gtc_di_new[5:length(bool_gtc_di_new)]), FALSE)
})

test_that("judge_gtc_sr, pattern 19 (bsrs n_dropin = 2)", {
  # Condition
  gtc <- c(10, 20, 30, 40)
  peaks <- c(8, 9, 9.2, 10, 11, 18, 19, 20, 21, 28, 29, 30, 31, 38, 39, 40, 41, 99)
  peaks_rm99_rmdi <- c(8, 9, 9.2, 10, 11, 18, 19, 20, 21, 28, 29, 30, 31, 38, 39, 40, 41)
  heights_rm99_rmdi <- c(100, 1000, 200, 10000, 500, 100, 1000, 10000, 500, 100, 1000, 10000, 500, 100, 1000, 10000, 500)
  bsrs <- c(0, 0, 0, 0, 0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0)
  fsrs <- c(0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0)
  dsrs <- c(0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0)
  m2srs <- c(0, 0, 0, 0.02, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  bsr_fltr <- 0.5
  fsr_fltr <- 0.5
  dsr_fltr <- 0.5
  m2sr_fltr <- 0.5
  n_dropin_fltr <- 2
  h_dropin_fltr <- 10000
  n_dropin <- 2
  bool_gtc_di <- rep(FALSE, length(peaks))
  bool_gtc_di[1] <- TRUE

  # Run
  bool_gtc_di_new <-  judge_gtc_sr(gtc, peaks, peaks_rm99_rmdi, heights_rm99_rmdi,
                                   bsrs, fsrs, dsrs, m2srs,
                                   bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr,
                                   n_dropin_fltr, h_dropin_fltr, n_dropin, bool_gtc_di)

  # Test
  expect_equal(length(bool_gtc_di_new), length(bool_gtc_di))
  expect_equal(bool_gtc_di_new[1], FALSE)
  expect_equal(bool_gtc_di_new[2], FALSE)
  expect_equal(bool_gtc_di_new[3], FALSE)
  expect_equal(bool_gtc_di_new[4], FALSE)
  expect_equal(any(bool_gtc_di_new[5:length(bool_gtc_di_new)]), FALSE)
})

test_that("judge_gtc_sr, pattern 20 (bsrs, two new drop-ins, n_dropin = 0)", {
  # Condition
  gtc <- c(10, 20, 30, 40)
  peaks <- c(8, 9, 9.2, 10, 11, 18, 19, 20, 21, 28, 29, 30, 31, 38, 39, 40, 41, 99)
  peaks_rm99_rmdi <- c(8, 9, 9.2, 10, 11, 18, 19, 20, 21, 28, 29, 30, 31, 38, 39, 40, 41)
  heights_rm99_rmdi <- c(100, 1000, 200, 10000, 500, 100, 1000, 10000, 500, 100, 1000, 10000, 500, 100, 1000, 10000, 500)
  bsrs <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0)
  fsrs <- c(0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0)
  dsrs <- c(0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0)
  m2srs <- c(0, 0, 0, 0.02, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  bsr_fltr <- 0.5
  fsr_fltr <- 0.5
  dsr_fltr <- 0.5
  m2sr_fltr <- 0.5
  n_dropin_fltr <- 2
  h_dropin_fltr <- 10000
  n_dropin <- 0
  bool_gtc_di <- rep(FALSE, length(peaks))
  bool_gtc_di[1] <- TRUE

  # Run
  bool_gtc_di_new <-  judge_gtc_sr(gtc, peaks, peaks_rm99_rmdi, heights_rm99_rmdi,
                                   bsrs, fsrs, dsrs, m2srs,
                                   bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr,
                                   n_dropin_fltr, h_dropin_fltr, n_dropin, bool_gtc_di)

  # Test
  expect_equal(length(bool_gtc_di_new), length(bool_gtc_di))
  expect_equal(bool_gtc_di_new[1], TRUE)
  expect_equal(bool_gtc_di_new[2], FALSE)
  expect_equal(bool_gtc_di_new[3], TRUE)
  expect_equal(bool_gtc_di_new[4], FALSE)
  expect_equal(bool_gtc_di_new[5], FALSE)
  expect_equal(bool_gtc_di_new[6], FALSE)
  expect_equal(bool_gtc_di_new[7], FALSE)
  expect_equal(bool_gtc_di_new[8], TRUE)
  expect_equal(any(bool_gtc_di_new[9:length(bool_gtc_di_new)]), FALSE)
})

test_that("judge_gtc_sr, pattern 21 (bsrs, two new drop-ins, n_dropin = 1)", {
  # Condition
  gtc <- c(10, 20, 30, 40)
  peaks <- c(8, 9, 9.2, 10, 11, 18, 19, 20, 21, 28, 29, 30, 31, 38, 39, 40, 41, 99)
  peaks_rm99_rmdi <- c(8, 9, 9.2, 10, 11, 18, 19, 20, 21, 28, 29, 30, 31, 38, 39, 40, 41)
  heights_rm99_rmdi <- c(100, 1000, 200, 10000, 500, 100, 1000, 10000, 500, 100, 1000, 10000, 500, 100, 1000, 10000, 500)
  bsrs <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0)
  fsrs <- c(0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0, 0, 0.05, 0, 0)
  dsrs <- c(0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01, 0, 0)
  m2srs <- c(0, 0, 0, 0.02, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  bsr_fltr <- 0.5
  fsr_fltr <- 0.5
  dsr_fltr <- 0.5
  m2sr_fltr <- 0.5
  n_dropin_fltr <- 2
  h_dropin_fltr <- 10000
  n_dropin <- 1
  bool_gtc_di <- rep(FALSE, length(peaks))
  bool_gtc_di[1] <- TRUE

  # Run
  bool_gtc_di_new <-  judge_gtc_sr(gtc, peaks, peaks_rm99_rmdi, heights_rm99_rmdi,
                                   bsrs, fsrs, dsrs, m2srs,
                                   bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr,
                                   n_dropin_fltr, h_dropin_fltr, n_dropin, bool_gtc_di)

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
  expect_equal(any(bool_gtc_di_new[9:length(bool_gtc_di_new)]), FALSE)
})
