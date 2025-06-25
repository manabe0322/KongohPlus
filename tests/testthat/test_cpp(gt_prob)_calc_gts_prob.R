test_that("calc_gts_prob, HWE", {
  # Condition
  gtc <- c(11, 11, 12, 13)
  pop_al <- c(11, 12, 13, 99)
  pop_prob <- c(0.1, 0.2, 0.3, 0.4)
  theta <- 0
  prof_ref <- numeric(0)
  prof_ref_99 <- numeric(0)
  hyp_id <- c(0, 0)

  # Run
  n_obs <- make_n_obs(pop_al, prof_ref)
  obs_al <- as.numeric(names(n_obs))
  gtc_prob <- calc_gts_prob(gtc, pop_al, pop_prob, theta, prof_ref_99, hyp_id, n_obs, obs_al)

  # Test
  expect_equal(gtc_prob, 0.01 * 0.12)
})

test_that("calc_gts_prob, HWE, drop-out, pattern1", {
  # Condition
  gtc <- c(11, 11, 13, 99)
  pop_al <- c(11, 12, 13, 99)
  pop_prob <- c(0.1, 0.2, 0.3, 0.4)
  theta <- 0
  prof_ref <- c(13, 14)
  prof_ref_99 <- c(13, 99)
  hyp_id <- c(0, 1)

  # Run
  n_obs <- make_n_obs(pop_al, prof_ref)
  obs_al <- as.numeric(names(n_obs))
  gtc_prob <- calc_gts_prob(gtc, pop_al, pop_prob, theta, prof_ref_99, hyp_id, n_obs, obs_al)

  # Test
  expect_equal(gtc_prob, 0.01)
})

test_that("calc_gts_prob, HWE, drop-out, pattern2", {
  # Condition
  gtc <- c(11, 11, 13, 99)
  pop_al <- c(11, 12, 13, 99)
  pop_prob <- c(0.1, 0.2, 0.3, 0.4)
  theta <- 0
  prof_ref <- c(13, 14)
  prof_ref_99 <- c(12, 99)
  hyp_id <- c(0, 1)

  # Run
  n_obs <- make_n_obs(pop_al, prof_ref)
  obs_al <- as.numeric(names(n_obs))
  gtc_prob <- calc_gts_prob(gtc, pop_al, pop_prob, theta, prof_ref_99, hyp_id, n_obs, obs_al)

  # Test
  expect_equal(gtc_prob, 0)
})

test_that("calc_gts_prob, theta = 0.01, pattern1", {
  # Condition
  gtc <- c(11, 11, 12, 13)
  pop_al <- c(11, 12, 13, 99)
  pop_prob <- c(0.1, 0.2, 0.3, 0.4)
  theta <- 0.01
  prof_ref <- numeric(0)
  prof_ref_99 <- numeric(0)
  hyp_id <- c(0, 0)

  # Run
  n_obs <- make_n_obs(pop_al, prof_ref)
  obs_al <- as.numeric(names(n_obs))
  gtc_prob <- calc_gts_prob(gtc, pop_al, pop_prob, theta, prof_ref_99, hyp_id, n_obs, obs_al)

  # Test
  expect_equal(gtc_prob, 0.1 * (theta + (1 - theta) * 0.1) * 2 * ((1 - theta) * 0.2) / (1 + theta) * ((1 - theta) * 0.3) / (1 + 2 * theta))
})

test_that("calc_gts_prob, theta = 0.01, pattern2", {
  # Condition
  gtc <- c(11, 11, 12, 13)
  pop_al <- c(11, 12, 13, 99)
  pop_prob <- c(0.1, 0.2, 0.3, 0.4)
  theta <- 0.01
  prof_ref <- c(11, 11, 12, 13)
  prof_ref_99 <- c(11, 11, 12, 13)
  hyp_id <- c(0, 0)

  # Run
  n_obs <- make_n_obs(pop_al, prof_ref)
  obs_al <- as.numeric(names(n_obs))
  gtc_prob <- calc_gts_prob(gtc, pop_al, pop_prob, theta, prof_ref_99, hyp_id, n_obs, obs_al)

  # Test
  expect_equal(gtc_prob, (2 * theta + (1 - theta) * 0.1) / (1 + 3 * theta) * (3 * theta + (1 - theta) * 0.1) / (1 + 4 * theta) * 2 * (theta + (1 - theta) * 0.2) / (1 + 5 * theta) * (theta + (1 - theta) * 0.3) / (1 + 6 * theta))
})

test_that("calc_gts_prob, theta = 0.01, pattern3", {
  # Condition
  gtc <- c(11, 11, 12, 13)
  pop_al <- c(11, 12, 13, 99)
  pop_prob <- c(0.1, 0.2, 0.3, 0.4)
  theta <- 0.01
  prof_ref <- c(11, 11, 12, 13)
  prof_ref_99 <- c(11, 11, 12, 13)
  hyp_id <- c(1, 0)

  # Run
  n_obs <- make_n_obs(pop_al, prof_ref)
  obs_al <- as.numeric(names(n_obs))
  gtc_prob <- calc_gts_prob(gtc, pop_al, pop_prob, theta, prof_ref_99, hyp_id, n_obs, obs_al)

  # Test
  expect_equal(gtc_prob, 2 * (theta + (1 - theta) * 0.2) / (1 + 3 * theta) * (theta + (1 - theta) * 0.3) / (1 + 4 * theta))
})

test_that("calc_gts_prob, theta = 0.01, pattern4", {
  # Condition
  gtc <- c(11, 11, 12, 13)
  pop_al <- c(11, 12, 13, 99)
  pop_prob <- c(0.1, 0.2, 0.3, 0.4)
  theta <- 0.01
  prof_ref <- c(11, 11, 12, 13)
  prof_ref_99 <- c(11, 11, 12, 13)
  hyp_id <- c(0, 1)

  # Run
  n_obs <- make_n_obs(pop_al, prof_ref)
  obs_al <- as.numeric(names(n_obs))
  gtc_prob <- calc_gts_prob(gtc, pop_al, pop_prob, theta, prof_ref_99, hyp_id, n_obs, obs_al)

  # Test
  expect_equal(gtc_prob, 0)
})

test_that("calc_gts_prob, theta = 0.01, pattern5", {
  # Condition
  gtc <- c(11, 11, 12, 13)
  pop_al <- c(11, 12, 13, 99)
  pop_prob <- c(0.1, 0.2, 0.3, 0.4)
  theta <- 0.01
  prof_ref <- c(11, 11, 12, 13)
  prof_ref_99 <- c(11, 11, 12, 13)
  hyp_id <- c(2, 0)

  # Run
  n_obs <- make_n_obs(pop_al, prof_ref)
  obs_al <- as.numeric(names(n_obs))
  gtc_prob <- calc_gts_prob(gtc, pop_al, pop_prob, theta, prof_ref_99, hyp_id, n_obs, obs_al)

  # Test
  expect_equal(gtc_prob, 0)
})

test_that("calc_gts_prob, theta = 0.01, pattern6", {
  # Condition
  gtc <- c(11, 11, 12, 13)
  pop_al <- c(11, 12, 13, 99)
  pop_prob <- c(0.1, 0.2, 0.3, 0.4)
  theta <- 0.01
  prof_ref <- c(11, 11, 12, 13)
  prof_ref_99 <- c(11, 11, 12, 13)
  hyp_id <- c(0, 2)

  # Run
  n_obs <- make_n_obs(pop_al, prof_ref)
  obs_al <- as.numeric(names(n_obs))
  gtc_prob <- calc_gts_prob(gtc, pop_al, pop_prob, theta, prof_ref_99, hyp_id, n_obs, obs_al)

  # Test
  expect_equal(gtc_prob, (2 * theta + (1 - theta) * 0.1) / (1 + 3 * theta) * (3 * theta + (1 - theta) * 0.1) / (1 + 4 * theta))
})

test_that("calc_gts_prob, theta = 0.01, drop-out, pattern1", {
  # Condition
  gtc <- c(11, 11, 13, 99)
  pop_al <- c(11, 12, 13, 99)
  pop_prob <- c(0.1, 0.2, 0.3, 0.4)
  theta <- 0.01
  prof_ref <- c(11, 11, 13, 14)
  prof_ref_99 <- c(11, 11, 13, 99)
  hyp_id <- c(0, 0)

  # Run
  n_obs <- make_n_obs(pop_al, prof_ref)
  obs_al <- as.numeric(names(n_obs))
  gtc_prob <- calc_gts_prob(gtc, pop_al, pop_prob, theta, prof_ref_99, hyp_id, n_obs, obs_al)

  # Test
  expect_equal(gtc_prob, (2 * theta + (1 - theta) * 0.1) / (1 + 3 * theta) * (3 * theta + (1 - theta) * 0.1) / (1 + 4 * theta) * 2 * (theta + (1 - theta) * 0.3) / (1 + 5 * theta) * ((1 - theta) * 0.4) / (1 + 6 * theta))
})

test_that("calc_gts_prob, theta = 0.01, drop-out, pattern2", {
  # Condition
  gtc <- c(11, 11, 99, 99)
  pop_al <- c(11, 12, 13, 99)
  pop_prob <- c(0.1, 0.2, 0.3, 0.4)
  theta <- 0.01
  prof_ref <- c(11, 11, 14, 15)
  prof_ref_99 <- c(11, 11, 99, 99)
  hyp_id <- c(0, 0)

  # Run
  n_obs <- make_n_obs(pop_al, prof_ref)
  obs_al <- as.numeric(names(n_obs))
  gtc_prob <- calc_gts_prob(gtc, pop_al, pop_prob, theta, prof_ref_99, hyp_id, n_obs, obs_al)

  # Test
  expect_equal(gtc_prob, (2 * theta + (1 - theta) * 0.1) / (1 + 3 * theta) * (3 * theta + (1 - theta) * 0.1) / (1 + 4 * theta) * ((1 - theta) * 0.4) / (1 + 5 * theta) * (theta + (1 - theta) * 0.4) / (1 + 6 * theta))
})

