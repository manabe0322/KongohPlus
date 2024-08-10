test_that("judge_gtc pattern 1", {
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
  adopt_gtc <- judge_gtc(gtc, peaks, heights, peaks_rm99, heights_rm99,
                         bsrs, fsrs, dsrs, m2srs,
                         bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr, n_dropin_fltr, h_dropin_fltr)

  # Test
  expect_equal(adopt_gtc, TRUE)
})
