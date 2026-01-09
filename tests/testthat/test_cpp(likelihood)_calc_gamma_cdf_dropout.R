test_that("calc_gamma_cdf_dropout", {
  # Condition
  eph_dedup1_dropout <- c(100, 200, 50, 200)
  scale_base_mk1 <- c(50, 30, 50)
  index_dedup1 <- c(0, 2, 1, 2)  # C++ scale
  index_ignore <- 2  # C++ scale
  at <- 30

  # Run
  log_cd_sum <- calc_gamma_cdf_dropout(eph_dedup1_dropout, scale_base_mk1, index_dedup1, index_ignore, at)

  # Test
  p1 <- pgamma(at, shape = 2, scale = 50)
  p2 <- pgamma(at, shape = 5/3, scale = 30)
  expect_equal(round(log_cd_sum, 7), round(log(p1) + log(p2), 7))
})
