test_that("remove_dropin, pattern 1", {
  # Condition
  peaks_rm99 <- c(16, 17, 18, 19)
  heights_rm99 <- c(100, 1000, 10000, 500)
  index_dropin <- integer(0)

  # Run
  tmp <- remove_dropin(peaks_rm99, heights_rm99, index_dropin)
  peaks_rm99_rmdi <- tmp[[1]]
  heights_rm99_rmdi <- tmp[[2]]

  # Test
  expect_equal(length(peaks_rm99_rmdi), 4)
  expect_equal(length(heights_rm99_rmdi), 4)
  expect_equal(peaks_rm99_rmdi[1], 16)
  expect_equal(peaks_rm99_rmdi[2], 17)
  expect_equal(peaks_rm99_rmdi[3], 18)
  expect_equal(peaks_rm99_rmdi[4], 19)
  expect_equal(heights_rm99_rmdi[1], 100)
  expect_equal(heights_rm99_rmdi[2], 1000)
  expect_equal(heights_rm99_rmdi[3], 10000)
  expect_equal(heights_rm99_rmdi[4], 500)
})

test_that("remove_dropin, pattern 2", {
  # Condition
  peaks_rm99 <- c(16, 17, 18, 19)
  heights_rm99 <- c(100, 1000, 10000, 500)
  index_dropin <- c(1, 3)

  # Run
  tmp <- remove_dropin(peaks_rm99, heights_rm99, index_dropin)
  peaks_rm99_rmdi <- tmp[[1]]
  heights_rm99_rmdi <- tmp[[2]]

  # Test
  expect_equal(length(peaks_rm99_rmdi), 2)
  expect_equal(length(heights_rm99_rmdi), 2)
  expect_equal(peaks_rm99_rmdi[1], 16)
  expect_equal(peaks_rm99_rmdi[2], 18)
  expect_equal(heights_rm99_rmdi[1], 100)
  expect_equal(heights_rm99_rmdi[2], 10000)
})
