test_that("update_heights_list_mk", {
  # Condition
  peaks_list_mk <- list()
  peaks_list_mk[[1]] <- c(13, 14, 15, 16)
  peaks_list_mk[[2]] <- c(14, 15, 16, 17)
  heights_list_mk <- list()
  heights_list_mk[[1]] <- c(50, 500, 600, 100)
  heights_list_mk[[2]] <- c(550, 650, 150, 25)
  peaks_all_rep <- sort(unique(unlist(peaks_list_mk)))

  # Run
  height_list_mk_new <- update_heights_list_mk(peaks_all_rep, peaks_list_mk, heights_list_mk)

  # Test
  expect_equal(length(peaks_all_rep), 5)
  expect_equal(peaks_all_rep[1], 13)
  expect_equal(peaks_all_rep[2], 14)
  expect_equal(peaks_all_rep[3], 15)
  expect_equal(peaks_all_rep[4], 16)
  expect_equal(peaks_all_rep[5], 17)
  expect_equal(length(height_list_mk_new), 2)
  expect_equal(length(height_list_mk_new[[1]]), 5)
  expect_equal(height_list_mk_new[[1]][1], 50)
  expect_equal(height_list_mk_new[[1]][2], 500)
  expect_equal(height_list_mk_new[[1]][3], 600)
  expect_equal(height_list_mk_new[[1]][4], 100)
  expect_equal(height_list_mk_new[[1]][5], 0)
  expect_equal(length(height_list_mk_new[[2]]), 5)
  expect_equal(height_list_mk_new[[2]][1], 0)
  expect_equal(height_list_mk_new[[2]][2], 550)
  expect_equal(height_list_mk_new[[2]][3], 650)
  expect_equal(height_list_mk_new[[2]][4], 150)
  expect_equal(height_list_mk_new[[2]][5], 25)
})
