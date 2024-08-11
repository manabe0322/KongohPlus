test_that("search_index_dropin, pattern 1", {
  # Condition
  peaks_rm99 <- c(16, 17, 18, 19.2)
  detectable_peak <- c(15, 16, 17, 18, 18.2, 19, 19.2, 20.2, 99)

  # Run
  index_dropin <- search_index_dropin(peaks_rm99, detectable_peak)

  # Test
  expect_equal(length(index_dropin), 0)
})

test_that("search_index_dropin, pattern 2", {
  # Condition
  peaks_rm99 <- c(16, 17, 18, 19.2)
  detectable_peak <- c(15, 16, 17, 18, 19)

  # Run
  index_dropin <- search_index_dropin(peaks_rm99, detectable_peak)

  # Test
  expect_equal(length(index_dropin), 1)
  expect_equal(index_dropin[1], 3)
})

test_that("search_index_dropin, pattern 3", {
  # Condition
  peaks_rm99 <- c(16, 17, 18, 19.2)
  detectable_peak <- c(15, 16, 17, 19)

  # Run
  index_dropin <- search_index_dropin(peaks_rm99, detectable_peak)

  # Test
  expect_equal(length(index_dropin), 2)
  expect_equal(index_dropin[1], 2)
  expect_equal(index_dropin[2], 3)
})
