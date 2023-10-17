test_that("calc_sr 1 repeat difference", {

  ##############
  # Conditions #
  ##############

  peaks <- c(12, 13, 13.2, 14, 14.2, 15, 16)
  heights <- c(100, 200, 150, 3000, 150, 3000, 100)
  gt <- c(14, 15)
  sr_consider <- rep(TRUE, 4)
  names(sr_consider) <- c("BS", "FS", "DS", "M2S")

  #######
  # Run #
  #######

  tmp <- calc_sr(peaks, heights, gt, sr_consider)
  bsr <- tmp[[1]]
  fsr <- tmp[[2]]
  dsr <- tmp[[3]]
  m2sr <- tmp[[4]]
  al <- tmp[[5]]
  tac <- tmp[[6]]

  ########
  # Test #
  ########

  expect_equal(as.numeric(bsr[1]), 0)
  expect_equal(as.numeric(bsr[2]), 0)
  expect_equal(as.numeric(fsr[1]), 0)
  expect_equal(as.numeric(fsr[2]), 0)
  expect_equal(as.numeric(dsr[1]), 0)
  expect_equal(as.numeric(dsr[2]), 0)
  expect_equal(as.numeric(m2sr[1]), 0)
  expect_equal(as.numeric(m2sr[2]), 0)
  expect_equal(as.numeric(al[1]), 0)
  expect_equal(as.numeric(al[2]), 0)
  expect_equal(as.numeric(tac[1]), 0)
  expect_equal(as.numeric(tac[2]), 0)
})

test_that("calc_sr 2 repeat difference", {

  ##############
  # Conditions #
  ##############

  peaks <- c(12, 13, 13.2, 14, 15, 15.2, 16, 17)
  heights <- c(100, 200, 150, 3000, 300, 150, 3000, 100)
  gt <- c(14, 16)
  sr_consider <- rep(TRUE, 4)
  names(sr_consider) <- c("BS", "FS", "DS", "M2S")

  #######
  # Run #
  #######

  tmp <- calc_sr(peaks, heights, gt, sr_consider)
  bsr <- tmp[[1]]
  fsr <- tmp[[2]]
  dsr <- tmp[[3]]
  m2sr <- tmp[[4]]
  al <- tmp[[5]]
  tac <- tmp[[6]]

  ########
  # Test #
  ########

  expect_equal(as.numeric(bsr[1]), 200 / 3000)
  expect_equal(as.numeric(bsr[2]), 0)
  expect_equal(as.numeric(fsr[1]), 0)
  expect_equal(as.numeric(fsr[2]), 100 / 3000)
  expect_equal(as.numeric(dsr[1]), 100 / 3000)
  expect_equal(as.numeric(dsr[2]), 0)
  expect_equal(as.numeric(m2sr[1]), 150 / 3000)
  expect_equal(as.numeric(m2sr[2]), 150 / 3000)
  expect_equal(as.numeric(al[1]), 14)
  expect_equal(as.numeric(al[2]), 16)
  expect_equal(as.numeric(tac[1]), 100 + 200 + 150 + 3000)
  expect_equal(as.numeric(tac[2]), 150 + 3000 + 100)
})

test_that("calc_sr 3 repeat difference", {

  ##############
  # Conditions #
  ##############

  peaks <- c(12, 13, 13.2, 14, 15, 16, 16.2, 17, 18)
  heights <- c(100, 200, 150, 3000, 200, 300, 150, 3000, 100)
  gt <- c(14, 17)
  sr_consider <- rep(TRUE, 4)
  names(sr_consider) <- c("BS", "FS", "DS", "M2S")

  #######
  # Run #
  #######

  tmp <- calc_sr(peaks, heights, gt, sr_consider)
  bsr <- tmp[[1]]
  fsr <- tmp[[2]]
  dsr <- tmp[[3]]
  m2sr <- tmp[[4]]
  al <- tmp[[5]]
  tac <- tmp[[6]]

  ########
  # Test #
  ########

  expect_equal(as.numeric(bsr[1]), 200 / 3000)
  expect_equal(as.numeric(bsr[2]), 300 / 3000)
  expect_equal(as.numeric(fsr[1]), 0)
  expect_equal(as.numeric(fsr[2]), 100 / 3000)
  expect_equal(as.numeric(dsr[1]), 100 / 3000)
  expect_equal(as.numeric(dsr[2]), 0)
  expect_equal(as.numeric(m2sr[1]), 150 / 3000)
  expect_equal(as.numeric(m2sr[2]), 150 / 3000)
  expect_equal(as.numeric(al[1]), 14)
  expect_equal(as.numeric(al[2]), 17)
  expect_equal(as.numeric(tac[1]), 100 + 200 + 150 + 3000)
  expect_equal(as.numeric(tac[2]), 300 + 150 + 3000 + 100)
})

test_that("calcSr 0.2 repeat difference", {

  ##############
  # Conditions #
  ##############

  peaks <- c(12, 12.2, 13, 13.2, 14, 14.2, 15, 15.2)
  heights <- c(200, 100, 400, 500, 4000, 4000, 300, 200)
  gt <- c(14, 14.2)
  sr_consider <- rep(TRUE, 4)
  names(sr_consider) <- c("BS", "FS", "DS", "M2S")

  #######
  # Run #
  #######

  tmp <- calc_sr(peaks, heights, gt, sr_consider)
  bsr <- tmp[[1]]
  fsr <- tmp[[2]]
  dsr <- tmp[[3]]
  m2sr <- tmp[[4]]
  al <- tmp[[5]]
  tac <- tmp[[6]]

  ########
  # Test #
  ########

  expect_equal(as.numeric(bsr[1]), 400 / 4000)
  expect_equal(as.numeric(bsr[2]), 0)
  expect_equal(as.numeric(fsr[1]), 300 / 4000)
  expect_equal(as.numeric(fsr[2]), 200 / 4000)
  expect_equal(as.numeric(dsr[1]), 200 / 4000)
  expect_equal(as.numeric(dsr[2]), 100 / 4000)
  expect_equal(as.numeric(m2sr[1]), 0)
  expect_equal(as.numeric(m2sr[2]), 0)
  expect_equal(as.numeric(al[1]), 14)
  expect_equal(as.numeric(al[2]), 14.2)
  expect_equal(as.numeric(tac[1]), 200 + 400 + 4000 + 300)
  expect_equal(as.numeric(tac[2]), 100 + 4000 + 200)
})

test_that("calcSr 0.8 repeat difference", {

  ##############
  # Conditions #
  ##############

  peaks <- c(12.2, 13, 13.2, 14, 14.2, 15, 15.2, 16)
  heights <- c(200, 100, 400, 500, 4000, 4000, 300, 200)
  gt <- c(14.2, 15)
  sr_consider <- rep(TRUE, 4)
  names(sr_consider) <- c("BS", "FS", "DS", "M2S")

  #######
  # Run #
  #######

  tmp <- calc_sr(peaks, heights, gt, sr_consider)
  bsr <- tmp[[1]]
  fsr <- tmp[[2]]
  dsr <- tmp[[3]]
  m2sr <- tmp[[4]]
  al <- tmp[[5]]
  tac <- tmp[[6]]

  ########
  # Test #
  ########

  expect_equal(as.numeric(bsr[1]), 400 / 4000)
  expect_equal(as.numeric(bsr[2]), 0)
  expect_equal(as.numeric(fsr[1]), 300 / 4000)
  expect_equal(as.numeric(fsr[2]), 200 / 4000)
  expect_equal(as.numeric(dsr[1]), 200 / 4000)
  expect_equal(as.numeric(dsr[2]), 100 / 4000)
  expect_equal(as.numeric(m2sr[1]), 0)
  expect_equal(as.numeric(m2sr[2]), 0)
  expect_equal(as.numeric(al[1]), 14.2)
  expect_equal(as.numeric(al[2]), 15)
  expect_equal(as.numeric(tac[1]), 200 + 400 + 4000 + 300)
  expect_equal(as.numeric(tac[2]), 100 + 4000 + 200)
})

test_that("calcSr homozygote", {

  ##############
  # Conditions #
  ##############

  peaks <- c(12, 13, 14, 15)
  heights <- c(200, 400, 4000, 300)
  gt <- c(14)
  sr_consider <- rep(TRUE, 4)
  names(sr_consider) <- c("BS", "FS", "DS", "M2S")

  #######
  # Run #
  #######

  tmp <- calc_sr(peaks, heights, gt, sr_consider)
  bsr <- tmp[[1]]
  fsr <- tmp[[2]]
  dsr <- tmp[[3]]
  m2sr <- tmp[[4]]
  al <- tmp[[5]]
  tac <- tmp[[6]]

  ########
  # Test #
  ########

  expect_equal(as.numeric(bsr[1]), 400 / 4000)
  expect_equal(as.numeric(bsr[2]), 0)
  expect_equal(as.numeric(fsr[1]), 300 / 4000)
  expect_equal(as.numeric(fsr[2]), 0)
  expect_equal(as.numeric(dsr[1]), 200 / 4000)
  expect_equal(as.numeric(dsr[2]), 0)
  expect_equal(as.numeric(m2sr[1]), 0)
  expect_equal(as.numeric(m2sr[2]), 0)
  expect_equal(as.numeric(al[1]), 14)
  expect_equal(as.numeric(al[2]), 0)
  expect_equal(as.numeric(tac[1]), 200 + 400 + 4000 + 300)
  expect_equal(as.numeric(tac[2]), 0)
})

test_that("calcSr drop-out", {

  ##############
  # Conditions #
  ##############

  peaks <- c(12, 13, 14, 15)
  heights <- c(200, 400, 4000, 300)
  gt <- c(14, 16)
  sr_consider <- rep(TRUE, 4)
  names(sr_consider) <- c("BS", "FS", "DS", "M2S")

  #######
  # Run #
  #######

  tmp <- calc_sr(peaks, heights, gt, sr_consider)
  bsr <- tmp[[1]]
  fsr <- tmp[[2]]
  dsr <- tmp[[3]]
  m2sr <- tmp[[4]]
  al <- tmp[[5]]
  tac <- tmp[[6]]

  ########
  # Test #
  ########

  expect_equal(as.numeric(bsr[1]), 400 / 4000)
  expect_equal(as.numeric(bsr[2]), 0)
  expect_equal(as.numeric(fsr[1]), 0)
  expect_equal(as.numeric(fsr[2]), 0)
  expect_equal(as.numeric(dsr[1]), 200 / 4000)
  expect_equal(as.numeric(dsr[2]), 0)
  expect_equal(as.numeric(m2sr[1]), 0)
  expect_equal(as.numeric(m2sr[2]), 0)
  expect_equal(as.numeric(al[1]), 14)
  expect_equal(as.numeric(al[2]), 0)
  expect_equal(as.numeric(tac[1]), 200 + 400 + 4000)
  expect_equal(as.numeric(tac[2]), 0)
})
