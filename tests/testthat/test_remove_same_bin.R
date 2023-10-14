test_that("remove_same_bin", {

  ##############
  # Conditions #
  ##############

  peaks <- as.character(c(11, 12, 12, 13, 14, 14, 15))
  sizes <- c(100, 103.8, 104.2, 108, 111.8, 112.2, 116)
  heights <- c(200, 120, 1000, 140, 100, 1000, 160)

  #######
  # Run #
  #######

  tmp <- remove_same_bin(peaks, sizes, heights)
  peaks_new <- tmp[[1]]
  sizes_new <- tmp[[2]]
  heights_new <- tmp[[3]]

  ########
  # Test #
  ########

  expect_equal(length(peaks_new), 3)
  expect_equal(peaks_new[1], "11")
  expect_equal(peaks_new[2], "13")
  expect_equal(peaks_new[3], "15")
  expect_equal(sizes_new[1], 100)
  expect_equal(sizes_new[2], 108)
  expect_equal(sizes_new[3], 116)
  expect_equal(heights_new[1], 200)
  expect_equal(heights_new[2], 140)
  expect_equal(heights_new[3], 160)
})
