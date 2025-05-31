test_that("make_detectable_peak pattern 1", {
  # Condition
  gtc <- c(13, 14, 15, 16)
  peaks <- c(12, 13, 14, 15, 16, 17)
  bsrs <- c(0, 0, 0, 0, 0, 0)
  fsrs <- c(0, 0, 0, 0, 0, 0)
  dsrs <- c(0, 0, 0, 0, 0, 0)
  m2srs <- c(0, 0, 0, 0, 0, 0)

  # Run
  detectable_peak <- make_detectable_peak(gtc, peaks, bsrs, fsrs, dsrs, m2srs)

  # Test
  expect_equal(length(detectable_peak), 4)
  expect_equal(detectable_peak[1], 13)
  expect_equal(detectable_peak[2], 14)
  expect_equal(detectable_peak[3], 15)
  expect_equal(detectable_peak[4], 16)
})

test_that("make_detectable_peak pattern 2", {
  # Condition
  gtc <- c(13, 14, 15, 16)
  peaks <- c(12, 13, 14, 15, 16, 17)
  bsrs <- c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
  fsrs <- c(0, 0, 0, 0, 0, 0)
  dsrs <- c(0, 0, 0, 0, 0, 0)
  m2srs <- c(0, 0, 0, 0, 0, 0)

  # Run
  detectable_peak <- make_detectable_peak(gtc, peaks, bsrs, fsrs, dsrs, m2srs)

  # Test
  expect_equal(length(detectable_peak), 5)
  expect_equal(detectable_peak[1], 12)
  expect_equal(detectable_peak[2], 13)
  expect_equal(detectable_peak[3], 14)
  expect_equal(detectable_peak[4], 15)
  expect_equal(detectable_peak[5], 16)
})

test_that("make_detectable_peak pattern 3", {
  # Condition
  gtc <- c(13, 14, 15, 16)
  peaks <- c(12, 13, 14, 15, 16, 17)
  bsrs <- c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
  fsrs <- c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
  dsrs <- c(0, 0, 0, 0, 0, 0)
  m2srs <- c(0, 0, 0, 0, 0, 0)

  # Run
  detectable_peak <- make_detectable_peak(gtc, peaks, bsrs, fsrs, dsrs, m2srs)

  # Test
  expect_equal(length(detectable_peak), 6)
  expect_equal(detectable_peak[1], 12)
  expect_equal(detectable_peak[2], 13)
  expect_equal(detectable_peak[3], 14)
  expect_equal(detectable_peak[4], 15)
  expect_equal(detectable_peak[5], 16)
  expect_equal(detectable_peak[6], 17)
})

test_that("make_detectable_peak pattern 4", {
  # Condition
  gtc <- c(13, 14, 15, 16)
  peaks <- c(12, 13, 14, 15, 16, 17)
  bsrs <- c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
  fsrs <- c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
  dsrs <- c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
  m2srs <- c(0, 0, 0, 0, 0, 0)

  # Run
  detectable_peak <- make_detectable_peak(gtc, peaks, bsrs, fsrs, dsrs, m2srs)

  # Test
  expect_equal(length(detectable_peak), 7)
  expect_equal(detectable_peak[1], 11)
  expect_equal(detectable_peak[2], 12)
  expect_equal(detectable_peak[3], 13)
  expect_equal(detectable_peak[4], 14)
  expect_equal(detectable_peak[5], 15)
  expect_equal(detectable_peak[6], 16)
  expect_equal(detectable_peak[7], 17)
})

test_that("make_detectable_peak pattern 5", {
  # Condition
  gtc <- c(13, 14, 15, 16)
  peaks <- c(12, 13, 14, 15, 16, 17)
  bsrs <- c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
  fsrs <- c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
  dsrs <- c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
  m2srs <- c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1)

  # Run
  detectable_peak <- make_detectable_peak(gtc, peaks, bsrs, fsrs, dsrs, m2srs)

  # Test
  expect_equal(length(detectable_peak), 11)
  expect_equal(detectable_peak[1], 11)
  expect_equal(detectable_peak[2], 12)
  expect_equal(detectable_peak[3], 12.2)
  expect_equal(detectable_peak[4], 13)
  expect_equal(detectable_peak[5], 13.2)
  expect_equal(detectable_peak[6], 14)
  expect_equal(detectable_peak[7], 14.2)
  expect_equal(detectable_peak[8], 15)
  expect_equal(detectable_peak[9], 15.2)
  expect_equal(detectable_peak[10], 16)
  expect_equal(detectable_peak[11], 17)
})

test_that("make_detectable_peak pattern 6 (If the assumed alleles are not detected, then the stutters are also not assumed to be detected.)", {
  # Condition
  gtc <- c(13, 14, 15, 16)
  peaks <- c(12, 17)
  bsrs <- c(0.1, 0.1)
  fsrs <- c(0.1, 0.1)
  dsrs <- c(0.1, 0.1)
  m2srs <- c(0.1, 0.1)

  # Run
  detectable_peak <- make_detectable_peak(gtc, peaks, bsrs, fsrs, dsrs, m2srs)

  # Test
  expect_equal(length(detectable_peak), 4)
  expect_equal(detectable_peak[1], 13)
  expect_equal(detectable_peak[2], 14)
  expect_equal(detectable_peak[3], 15)
  expect_equal(detectable_peak[4], 16)
})

test_that("make_detectable_peak pattern 7", {
  # Condition
  gtc <- c(13.2, 14.2, 15.2, 16)
  peaks <- c(13.2, 14.2, 15.2, 16)
  bsrs <- c(0.1, 0.1, 0.1, 0.1)
  fsrs <- c(0.1, 0.1, 0.1, 0.1)
  dsrs <- c(0.1, 0.1, 0.1, 0.1)
  m2srs <- c(0.1, 0.1, 0.1, 0.1)

  # Run
  detectable_peak <- make_detectable_peak(gtc, peaks, bsrs, fsrs, dsrs, m2srs)

  # Test
  expect_equal(length(detectable_peak), 11)
  expect_equal(detectable_peak[1], 11.2)
  expect_equal(detectable_peak[2], 12.2)
  expect_equal(detectable_peak[3], 13)
  expect_equal(detectable_peak[4], 13.2)
  expect_equal(detectable_peak[5], 14)
  expect_equal(detectable_peak[6], 14.2)
  expect_equal(detectable_peak[7], 15)
  expect_equal(detectable_peak[8], 15.2)
  expect_equal(detectable_peak[9], 16)
  expect_equal(detectable_peak[10], 16.2)
  expect_equal(detectable_peak[11], 17)
})
