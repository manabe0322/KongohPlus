test_that("calc_gt_prob, HWE, hetero", {
  # Condition
  gt <- c(12, 13)
  pop_al <- c(11, 12, 13, 99)
  pop_prob <- c(0.1, 0.2, 0.3, 0.4)
  theta <- 0
  prof_ref <- c(12, 13)

  # Run
  n_obs <- make_n_obs(pop_al, prof_ref)
  gt_prob <- calc_gt_prob(gt, pop_al, pop_prob, theta, n_obs)

  # Test
  expect_equal(gt_prob, 0.12)
})

test_that("calc_gt_prob, HWE, homo", {
  # Condition
  gt <- c(12, 12)
  pop_al <- c(11, 12, 13, 99)
  pop_prob <- c(0.1, 0.2, 0.3, 0.4)
  theta <- 0
  prof_ref <- c(12, 12)

  # Run
  n_obs <- make_n_obs(pop_al, prof_ref)
  gt_prob <- calc_gt_prob(gt, pop_al, pop_prob, theta, n_obs)

  # Test
  expect_equal(gt_prob, 0.04)
})

test_that("calc_gt_prob, theta = 0.01, hetero", {
  # Condition
  gt <- c(12, 13)
  pop_al <- c(11, 12, 13, 99)
  pop_prob <- c(0.1, 0.2, 0.3, 0.4)
  theta <- 0.01
  prof_ref <- c(12, 13)

  # Run
  n_obs <- make_n_obs(pop_al, prof_ref)
  gt_prob <- calc_gt_prob(gt, pop_al, pop_prob, theta, n_obs)

  # Test
  expect_equal(gt_prob, 2 * (theta + (1 - theta) * 0.2) * (theta + (1 - theta) * 0.3) / ((1 + theta) * (1 + 2 * theta)))
})

test_that("calc_gt_prob, theta = 0.01, homo", {
  # Condition
  gt <- c(12, 12)
  pop_al <- c(11, 12, 13, 99)
  pop_prob <- c(0.1, 0.2, 0.3, 0.4)
  theta <- 0.01
  prof_ref <- c(12, 12)

  # Run
  n_obs <- make_n_obs(pop_al, prof_ref)
  gt_prob <- calc_gt_prob(gt, pop_al, pop_prob, theta, n_obs)

  # Test
  expect_equal(gt_prob, (2 * theta + (1 - theta) * 0.2) * (3 * theta + (1 - theta) * 0.2) / ((1 + theta) * (1 + 2 * theta)))
})
