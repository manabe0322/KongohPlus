test_that("calc_gamma_param_cpp", {
  # Condition
  eph_dedup1 <- c(1000, 150)
  scale_base_mk1 <- c(50, 30)
  index_dedup1 <- c(0, 1) # C++ scale

  # Run
  gamma_param <- calc_gamma_param_cpp(eph_dedup1, scale_base_mk1, index_dedup1)

  # Test
  expect_equal(round(gamma_param[1], 5), 24.26606)
  expect_equal(round(gamma_param[2], 4), 47.3913)
})
