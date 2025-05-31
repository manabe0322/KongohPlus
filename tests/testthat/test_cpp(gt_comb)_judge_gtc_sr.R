test_that("judge_gtc_sr, pattern 1", {
  ## Two-person mixture with many peaks at stutter positions
  ## These peaks can be explained by stutters
  ## Therefore, the genotype combination can be adopted

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

test_that("judge_gtc_sr, pattern 2", {
  ## Two-person mixture with many peaks at stutter positions
  ## The back stutter ratio for allele 10 is zero -> peak 9 cannot be explained by the back stutter of allele 10
  ## Peak 9 cannot be also explained by a drop-in peak
  ## Therefore, the genotype combination cannot be adopted

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

test_that("judge_gtc_sr, pattern 3", {
  ## Two-person mixture with many peaks at stutter positions
  ## The forward stutter ratio for allele 10 is zero -> peak 11 cannot be explained by the forward stutter of allele 10
  ## Peak 11 cannot be also explained by a drop-in peak
  ## Therefore, the genotype combination cannot be adopted

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

test_that("judge_gtc_sr, pattern 4", {
  ## Two-person mixture with many peaks at stutter positions
  ## The double-back stutter ratio for allele 10 is zero -> peak 8 cannot be explained by the double-back stutter of allele 10
  ## Peak 8 cannot be also explained by a drop-in peak
  ## Therefore, the genotype combination cannot be adopted

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

test_that("judge_gtc_sr, pattern 5", {
  ## Two-person mixture with many peaks at stutter positions
  ## The minus 2-nt stutter ratio for allele 10 allele is zero -> peak 9.2 cannot be explained by the minus 2-nt stutter of allele 10
  ## Peak 9.2 cannot be also explained by a drop-in peak
  ## Therefore, the genotype combination cannot be adopted

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

test_that("judge_gtc_sr, pattern 6", {
  ## Two-person mixture with many peaks at stutter positions
  ## The back stutter ratio for allele 10 allele is greater than the stutter filter -> peak 9 cannot be explained by the back stutter of allele 10
  ## Peak 9 cannot be also explained by a drop-in peak
  ## Therefore, the genotype combination cannot be adopted

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

test_that("judge_gtc_sr, pattern 7", {
  ## Two-person mixture with many peaks at stutter positions
  ## The forward stutter ratio for allele 10 allele is greater than the stutter filter -> peak 11 cannot be explained by the forward stutter of allele 10
  ## Peak 11 cannot be also explained by a drop-in peak
  ## Therefore, the genotype combination cannot be adopted

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

test_that("judge_gtc_sr, pattern 8", {
  ## Two-person mixture with many peaks at stutter positions
  ## The double-back stutter ratio for allele 10 allele is greater than the stutter filter -> peak 8 cannot be explained by the double-back stutter of allele 10
  ## Peak 8 cannot be also explained by a drop-in peak
  ## Therefore, the genotype combination cannot be adopted

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

test_that("judge_gtc_sr, pattern 9", {
  ## Two-person mixture with many peaks at stutter positions
  ## The minus 2-nt stutter ratio for allele 10 allele is greater than the stutter filter -> peak 9.2 cannot be explained by the minus 2-nt stutter of allele 10
  ## Peak 9.2 cannot be also explained by a drop-in peak
  ## Therefore, the genotype combination cannot be adopted

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

test_that("judge_gtc_sr, pattern 10", {
  ## Two-person mixture with many peaks at stutter positions
  ## The back stutter ratio for allele 10 is zero -> peak 9 cannot be explained by the back stutter of allele 10
  ## However, peak 9 can be explained by a drop-in peak
  ## Therefore, the genotype combination can be adopted

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

test_that("judge_gtc_sr, pattern 11", {
  ## Two-person mixture with many peaks at stutter positions
  ## The forward stutter ratio for allele 10 is zero -> peak 11 cannot be explained by the forward stutter of allele 10
  ## However, peak 11 can be explained by a drop-in peak
  ## Therefore, the genotype combination can be adopted

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

test_that("judge_gtc_sr, pattern 12", {
  ## Two-person mixture with many peaks at stutter positions
  ## The double-back stutter ratio for allele 10 is zero -> peak 8 cannot be explained by the double-back stutter of allele 10
  ## However, peak 8 can be explained by a drop-in peak
  ## Therefore, the genotype combination can be adopted

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

test_that("judge_gtc_sr, pattern 13", {
  ## Two-person mixture with many peaks at stutter positions
  ## The minus 2-nt stutter ratio for allele 10 is zero -> peak 9.2 cannot be explained by the minus 2-nt stutter of allele 10
  ## However, peak 9.2 can be explained by a drop-in peak
  ## Therefore, the genotype combination can be adopted

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

test_that("judge_gtc_sr, pattern 14", {
  ## Two-person mixture with many peaks at stutter positions
  ## The back stutter ratio for allele 10 is greater than the stutter filter -> peak 9 cannot be explained by the back stutter of allele 10
  ## However, peak 9 can be explained by a drop-in peak
  ## Therefore, the genotype combination can be adopted

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

test_that("judge_gtc_sr, pattern 15", {
  ## Two-person mixture with many peaks at stutter positions
  ## The forward stutter ratio for allele 10 is greater than the stutter filter -> peak 11 cannot be explained by the forward stutter of allele 10
  ## However, peak 11 can be explained by a drop-in peak
  ## Therefore, the genotype combination can be adopted

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

test_that("judge_gtc_sr, pattern 16", {
  ## Two-person mixture with many peaks at stutter positions
  ## The double-back stutter ratio for allele 10 is greater than the stutter filter -> peak 8 cannot be explained by the double-back stutter of allele 10
  ## However, peak 8 can be explained by a drop-in peak
  ## Therefore, the genotype combination can be adopted

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

test_that("judge_gtc_sr, pattern 17", {
  ## Two-person mixture with many peaks at stutter positions
  ## The minus 2-nt stutter ratio for allele 10 is greater than the stutter filter -> peak 9.2 cannot be explained by the minus 2-nt stutter of allele 10
  ## However, peak 9.2 can be explained by a drop-in peak
  ## Therefore, the genotype combination can be adopted

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

test_that("judge_gtc_sr, pattern 18", {
  ## Two-person mixture with many peaks at stutter positions
  ## The back stutter ratio for allele 10 is zero -> peak 9 cannot be explained by the back stutter of allele 10
  ## Peak 9 can be explained by a drop-in peak
  ## The number of drop-in peaks does not exceed the upper limit
  ## Therefore, the genotype combination can be adopted

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

test_that("judge_gtc_sr, pattern 19", {
  ## Two-person mixture with many peaks at stutter positions
  ## The back stutter ratio for allele 10 is zero -> peak 9 cannot be explained by the back stutter of allele 10
  ## Peak 9 can be explained by a drop-in peak
  ## However, the number of drop-in peaks does exceed the upper limit
  ## Therefore, the genotype combination cannot be adopted

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

test_that("judge_gtc_sr, pattern 20 (bsr = 0, two new drop-ins, n_dropin = 0)", {
  ## Two-person mixture with many peaks at stutter positions
  ## The back stutter ratios for allele 10 and allele 20 are zero -> peak 9 and peak 19 cannot be explained by the back stutters
  ## Peak 9 and peak 19 can be explained by drop-in peaks
  ## The number of drop-in peaks does not exceed the upper limit
  ## Therefore, the genotype combination can be adopted

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

test_that("judge_gtc_sr, pattern 21 (bsr = 0, two new drop-ins, n_dropin = 1)", {
  ## Two-person mixture with many peaks at stutter positions
  ## The back stutter ratios for allele 10 and allele 20 are zero -> peak 9 and peak 19 cannot be explained by the back stutters
  ## Peak 9 and peak 19 can be explained by drop-in peaks
  ## However, the number of drop-in peaks does exceed the upper limit
  ## Therefore, the genotype combination cannot be adopted

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
