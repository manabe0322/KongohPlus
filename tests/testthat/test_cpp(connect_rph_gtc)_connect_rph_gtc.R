test_that("connect_rph_gtc", {
  # Condition
  gtc <- c(13, 14, 15, 16)
  bs_gtc <- c(12, 13, 14, 15)
  fs_gtc <- c(14, 15, 16, 17)
  ds_gtc <- c(11, 12, 13, 14)
  m2s_gtc <- c(12.2, 13.2, 14.2, 15.2)
  peaks <- c(11, 12, 13, 14, 15, 16, 18, 99)
  bsrs <- c(rep(0.1, 7), 0)
  fsrs <- c(rep(0.01, 7), 0)
  dsrs <- c(rep(0.005, 7), 0)
  m2srs <- c(rep(0.05, 7), 0)
  bool_gtc_di <- c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE)

  # Run
  index_rph_gtc_cpp <- connect_rph_gtc(gtc, bs_gtc, fs_gtc, ds_gtc, m2s_gtc, peaks, bsrs, fsrs, dsrs, m2srs, bool_gtc_di)

  # Test
  expect_equal(length(index_rph_gtc_cpp), 8)
  ## for peak 11
  expect_equal(length(index_rph_gtc_cpp[[1]]), 4)
  expect_equal(index_rph_gtc_cpp[[1]][1], 3 * 8 + 0)
  expect_equal(index_rph_gtc_cpp[[1]][2], 5 * 8 + 0)
  expect_equal(index_rph_gtc_cpp[[1]][3], 5 * 8 + 0)
  expect_equal(index_rph_gtc_cpp[[1]][4], 5 * 8 + 0)
  ## for peak 12
  expect_equal(length(index_rph_gtc_cpp[[2]]), 4)
  expect_equal(index_rph_gtc_cpp[[2]][1], 1 * 8 + 1)
  expect_equal(index_rph_gtc_cpp[[2]][2], 3 * 8 + 1)
  expect_equal(index_rph_gtc_cpp[[2]][3], 5 * 8 + 0)
  expect_equal(index_rph_gtc_cpp[[2]][4], 5 * 8 + 0)
  ## for peak 13
  expect_equal(length(index_rph_gtc_cpp[[3]]), 4)
  expect_equal(index_rph_gtc_cpp[[3]][1], 0 * 8 + 2)
  expect_equal(index_rph_gtc_cpp[[3]][2], 1 * 8 + 2)
  expect_equal(index_rph_gtc_cpp[[3]][3], 3 * 8 + 2)
  expect_equal(index_rph_gtc_cpp[[3]][4], 5 * 8 + 0)
  ## for peak 14
  expect_equal(length(index_rph_gtc_cpp[[4]]), 4)
  expect_equal(index_rph_gtc_cpp[[4]][1], 0 * 8 + 3)
  expect_equal(index_rph_gtc_cpp[[4]][2], 2 * 8 + 3)
  expect_equal(index_rph_gtc_cpp[[4]][3], 1 * 8 + 3)
  expect_equal(index_rph_gtc_cpp[[4]][4], 3 * 8 + 3)
  ## for peak 15
  expect_equal(length(index_rph_gtc_cpp[[5]]), 4)
  expect_equal(index_rph_gtc_cpp[[5]][1], 2 * 8 + 4)
  expect_equal(index_rph_gtc_cpp[[5]][2], 5 * 8 + 0)
  expect_equal(index_rph_gtc_cpp[[5]][3], 0 * 8 + 4)
  expect_equal(index_rph_gtc_cpp[[5]][4], 1 * 8 + 4)
  ## for peak 16
  expect_equal(length(index_rph_gtc_cpp[[6]]), 4)
  expect_equal(index_rph_gtc_cpp[[6]][1], 5 * 8 + 0)
  expect_equal(index_rph_gtc_cpp[[6]][2], 5 * 8 + 0)
  expect_equal(index_rph_gtc_cpp[[6]][3], 0 * 8 + 5)
  expect_equal(index_rph_gtc_cpp[[6]][4], 2 * 8 + 5)
  ## for peak 18
  expect_equal(length(index_rph_gtc_cpp[[7]]), 4)
  expect_equal(index_rph_gtc_cpp[[7]][1], 5 * 8 + 1)
  expect_equal(index_rph_gtc_cpp[[7]][2], 5 * 8 + 1)
  expect_equal(index_rph_gtc_cpp[[7]][3], 5 * 8 + 1)
  expect_equal(index_rph_gtc_cpp[[7]][4], 5 * 8 + 1)
  ## for peak 99
  expect_equal(length(index_rph_gtc_cpp[[8]]), 4)
  expect_equal(index_rph_gtc_cpp[[8]][1], 5 * 8 + 0)
  expect_equal(index_rph_gtc_cpp[[8]][2], 5 * 8 + 0)
  expect_equal(index_rph_gtc_cpp[[8]][3], 5 * 8 + 0)
  expect_equal(index_rph_gtc_cpp[[8]][4], 5 * 8 + 0)
})
