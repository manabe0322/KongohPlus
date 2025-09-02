test_that("judge_gt_comb pattern 1 (no drop-in)", {
  # Condition
  gtc <- c(26.2, 28.2, 21, 31.2)
  peaks <- c(20, 20.2, 21, 26.2, 27.2, 28.2, 30.2, 31.2, 32.2, 99)
  heights <- c(884, 182, 8250, 500, 134, 164, 841, 7329, 253, 0)
  peaks_rm99 <- c(20, 20.2, 21, 26.2, 27.2, 28.2, 30.2, 31.2, 32.2)
  heights_rm99 <- c(884, 182, 8250, 500, 134, 164, 841, 7329, 253)
  bsrs <- c(0.0953, 0.053, 0.1011, 0.0896, 0.0961, 0.1017, 0.1142, 0.1204, 0.124, 0)
  fsrs <- c(0.0063, 0.0063, 0.0063, 0.0063, 0.0063, 0.0063, 0.0063, 0.0063, 0.0063, 0)
  dsrs <- c(0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0)
  m2srs <- c(0.0191, 0.0191, 0.0191, 0.0191, 0.0191, 0.0191, 0.0191, 0.0191, 0.0191, 0)
  bsr_fltr <- 0.85
  fsr_fltr <- 0.35
  dsr_fltr <- 0.2
  m2sr_fltr <- 0.1
  n_dropin_fltr <- 2
  h_dropin_fltr <- 150

  # Run
  bool_gtc_di <- judge_gt_comb(gtc, peaks, heights, peaks_rm99, heights_rm99,
                               bsrs, fsrs, dsrs, m2srs,
                               bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr, n_dropin_fltr, h_dropin_fltr)

  # Test
  expect_equal(length(bool_gtc_di), length(peaks) + 1)
  expect_equal(bool_gtc_di[1], TRUE)
  expect_equal(any(bool_gtc_di[2:length(bool_gtc_di)]), FALSE)
})

test_that("judge_gt_comb pattern 2 (one drop-in other than stutter positions)", {
  # Condition
  gtc <- c(26.2, 28.2, 21, 31.2)
  peaks <- c(20, 20.2, 21, 26.2, 27.2, 28.2, 30.2, 31.2, 32.2, 35, 99)
  heights <- c(884, 182, 8250, 500, 134, 164, 841, 7329, 253, 100, 0)
  peaks_rm99 <- c(20, 20.2, 21, 26.2, 27.2, 28.2, 30.2, 31.2, 32.2, 35)
  heights_rm99 <- c(884, 182, 8250, 500, 134, 164, 841, 7329, 253, 100)
  bsrs <- c(0.0953, 0.053, 0.1011, 0.0896, 0.0961, 0.1017, 0.1142, 0.1204, 0.124, 0.01, 0)
  fsrs <- c(0.0063, 0.0063, 0.0063, 0.0063, 0.0063, 0.0063, 0.0063, 0.0063, 0.0063, 0.01, 0)
  dsrs <- c(0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.01, 0)
  m2srs <- c(0.0191, 0.0191, 0.0191, 0.0191, 0.0191, 0.0191, 0.0191, 0.0191, 0.0191, 0.01, 0)
  bsr_fltr <- 0.85
  fsr_fltr <- 0.35
  dsr_fltr <- 0.2
  m2sr_fltr <- 0.1
  n_dropin_fltr <- 2
  h_dropin_fltr <- 150

  # Run
  bool_gtc_di <- judge_gt_comb(gtc, peaks, heights, peaks_rm99, heights_rm99,
                               bsrs, fsrs, dsrs, m2srs,
                               bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr, n_dropin_fltr, h_dropin_fltr)

  # Test
  expect_equal(length(bool_gtc_di), length(peaks) + 1)
  expect_equal(bool_gtc_di[1], TRUE)
  expect_equal(any(bool_gtc_di[2:10]), FALSE)
  expect_equal(bool_gtc_di[11], TRUE)
  expect_equal(bool_gtc_di[12], FALSE)
})

test_that("judge_gt_comb pattern 3 (two drop-ins other than stutter positions)", {
  # Condition
  gtc <- c(26.2, 28.2, 21, 31.2)
  peaks <- c(20, 20.2, 21, 26.2, 27.2, 28.2, 30.2, 31.2, 32.2, 35, 40, 99)
  heights <- c(884, 182, 8250, 500, 134, 164, 841, 7329, 253, 100, 100, 0)
  peaks_rm99 <- c(20, 20.2, 21, 26.2, 27.2, 28.2, 30.2, 31.2, 32.2, 35, 40)
  heights_rm99 <- c(884, 182, 8250, 500, 134, 164, 841, 7329, 253, 100, 100)
  bsrs <- c(0.0953, 0.053, 0.1011, 0.0896, 0.0961, 0.1017, 0.1142, 0.1204, 0.124, 0.01, 0.01, 0)
  fsrs <- c(0.0063, 0.0063, 0.0063, 0.0063, 0.0063, 0.0063, 0.0063, 0.0063, 0.0063, 0.01, 0.01, 0)
  dsrs <- c(0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.01, 0.01, 0)
  m2srs <- c(0.0191, 0.0191, 0.0191, 0.0191, 0.0191, 0.0191, 0.0191, 0.0191, 0.0191, 0.01, 0.01, 0)
  bsr_fltr <- 0.85
  fsr_fltr <- 0.35
  dsr_fltr <- 0.2
  m2sr_fltr <- 0.1
  n_dropin_fltr <- 2
  h_dropin_fltr <- 150

  # Run
  bool_gtc_di <- judge_gt_comb(gtc, peaks, heights, peaks_rm99, heights_rm99,
                               bsrs, fsrs, dsrs, m2srs,
                               bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr, n_dropin_fltr, h_dropin_fltr)

  # Test
  expect_equal(length(bool_gtc_di), length(peaks) + 1)
  expect_equal(bool_gtc_di[1], TRUE)
  expect_equal(any(bool_gtc_di[2:10]), FALSE)
  expect_equal(bool_gtc_di[11], TRUE)
  expect_equal(bool_gtc_di[12], TRUE)
  expect_equal(bool_gtc_di[13], FALSE)
})

test_that("judge_gt_comb pattern 4 (one drop-in at the stutter position of allele 28.2, over_bsr & over_fsr = TRUE of peak 27.2)", {
  # Condition
  gtc <- c(26.2, 28.2, 21, 31.2)
  peaks <- c(20, 20.2, 21, 26.2, 27.2, 28.2, 30.2, 31.2, 32.2, 99)
  heights <- c(884, 182, 8250, 500, 150, 164, 841, 7329, 253, 0)
  peaks_rm99 <- c(20, 20.2, 21, 26.2, 27.2, 28.2, 30.2, 31.2, 32.2)
  heights_rm99 <- c(884, 182, 8250, 500, 150, 164, 841, 7329, 253)
  bsrs <- c(0.0953, 0.053, 0.1011, 0.0896, 0.0961, 0.1017, 0.1142, 0.1204, 0.124, 0)
  fsrs <- c(0.0063, 0.0063, 0.0063, 0.0063, 0.0063, 0.0063, 0.0063, 0.0063, 0.0063, 0)
  dsrs <- c(0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0)
  m2srs <- c(0.0191, 0.0191, 0.0191, 0.0191, 0.0191, 0.0191, 0.0191, 0.0191, 0.0191, 0)
  bsr_fltr <- 0.85
  fsr_fltr <- 0.25
  dsr_fltr <- 0.2
  m2sr_fltr <- 0.1
  n_dropin_fltr <- 2
  h_dropin_fltr <- 150

  # Run
  bool_gtc_di <- judge_gt_comb(gtc, peaks, heights, peaks_rm99, heights_rm99,
                               bsrs, fsrs, dsrs, m2srs,
                               bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr, n_dropin_fltr, h_dropin_fltr)

  # Test
  expect_equal(length(bool_gtc_di), length(peaks) + 1)
  expect_equal(bool_gtc_di[1], TRUE)
  expect_equal(any(bool_gtc_di[2:5]), FALSE)
  expect_equal(bool_gtc_di[6], TRUE)
  expect_equal(any(bool_gtc_di[7:length(bool_gtc_di)]), FALSE)
})

test_that("judge_gt_comb pattern 5 (pattern 2 + pattern 4)", {
  # Condition
  gtc <- c(26.2, 28.2, 21, 31.2)
  peaks <- c(20, 20.2, 21, 26.2, 27.2, 28.2, 30.2, 31.2, 32.2, 35, 99)
  heights <- c(884, 182, 8250, 500, 150, 164, 841, 7329, 253, 100, 0)
  peaks_rm99 <- c(20, 20.2, 21, 26.2, 27.2, 28.2, 30.2, 31.2, 32.2, 35)
  heights_rm99 <- c(884, 182, 8250, 500, 150, 164, 841, 7329, 253, 100)
  bsrs <- c(0.0953, 0.053, 0.1011, 0.0896, 0.0961, 0.1017, 0.1142, 0.1204, 0.124, 0.01, 0)
  fsrs <- c(0.0063, 0.0063, 0.0063, 0.0063, 0.0063, 0.0063, 0.0063, 0.0063, 0.0063, 0.01, 0)
  dsrs <- c(0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.01, 0)
  m2srs <- c(0.0191, 0.0191, 0.0191, 0.0191, 0.0191, 0.0191, 0.0191, 0.0191, 0.0191, 0.01, 0)
  bsr_fltr <- 0.85
  fsr_fltr <- 0.25
  dsr_fltr <- 0.2
  m2sr_fltr <- 0.1
  n_dropin_fltr <- 2
  h_dropin_fltr <- 150

  # Run
  bool_gtc_di <- judge_gt_comb(gtc, peaks, heights, peaks_rm99, heights_rm99,
                               bsrs, fsrs, dsrs, m2srs,
                               bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr, n_dropin_fltr, h_dropin_fltr)

  # Test
  expect_equal(length(bool_gtc_di), length(peaks) + 1)
  expect_equal(bool_gtc_di[1], TRUE)
  expect_equal(any(bool_gtc_di[2:5]), FALSE)
  expect_equal(bool_gtc_di[6], TRUE)
  expect_equal(any(bool_gtc_di[7:10]), FALSE)
  expect_equal(bool_gtc_di[11], TRUE)
  expect_equal(bool_gtc_di[12], FALSE)
})

test_that("judge_gt_comb pattern 5 (pattern 3 + pattern 4)", {
  # Condition
  gtc <- c(26.2, 28.2, 21, 31.2)
  peaks <- c(20, 20.2, 21, 26.2, 27.2, 28.2, 30.2, 31.2, 32.2, 35, 40, 99)
  heights <- c(884, 182, 8250, 500, 150, 164, 841, 7329, 253, 100, 100, 0)
  peaks_rm99 <- c(20, 20.2, 21, 26.2, 27.2, 28.2, 30.2, 31.2, 32.2, 35, 40)
  heights_rm99 <- c(884, 182, 8250, 500, 150, 164, 841, 7329, 253, 100, 100)
  bsrs <- c(0.0953, 0.053, 0.1011, 0.0896, 0.0961, 0.1017, 0.1142, 0.1204, 0.124, 0.01, 0.01, 0)
  fsrs <- c(0.0063, 0.0063, 0.0063, 0.0063, 0.0063, 0.0063, 0.0063, 0.0063, 0.0063, 0.01, 0.01, 0)
  dsrs <- c(0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.01, 0.01, 0)
  m2srs <- c(0.0191, 0.0191, 0.0191, 0.0191, 0.0191, 0.0191, 0.0191, 0.0191, 0.0191, 0.01, 0.01, 0)
  bsr_fltr <- 0.85
  fsr_fltr <- 0.25
  dsr_fltr <- 0.2
  m2sr_fltr <- 0.1
  n_dropin_fltr <- 2
  h_dropin_fltr <- 150

  # Run
  bool_gtc_di <- judge_gt_comb(gtc, peaks, heights, peaks_rm99, heights_rm99,
                               bsrs, fsrs, dsrs, m2srs,
                               bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr, n_dropin_fltr, h_dropin_fltr)

  # Test
  expect_equal(length(bool_gtc_di), length(peaks) + 1)
  expect_equal(bool_gtc_di[1], FALSE)
  expect_equal(any(bool_gtc_di[2:10]), FALSE)
  expect_equal(bool_gtc_di[11], TRUE)
  expect_equal(bool_gtc_di[12], TRUE)
  expect_equal(bool_gtc_di[13], FALSE)
})
