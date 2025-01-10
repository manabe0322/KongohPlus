test_that("calc_gtc_prob", {
  # Condition
  gtc <- c(11, 11, 12, 13)
  pop_al <- c(11, 12, 13, 99)
  pop_prob <- c(0.1, 0.2, 0.3, 0.4)
  theta <- 0.01
  prof_ref <- c(11, 11, 12, 13)
  hyp_id <- c(0, 0)

  # Run
  n_obs <- make_n_obs(pop_al, prof_ref)
  gt_prob_1 <- calc_gt_prob(gtc[1:2], pop_al, pop_prob, theta, n_obs)
  gt_prob_2 <- calc_gt_prob(gtc[3:4], pop_al, pop_prob, theta, n_obs)
  gtc_prob <- calc_gtc_prob(gtc, pop_al, pop_prob, theta, prof_ref, hyp_id, n_obs)

  # Test
  expect_equal(gtc_prob, gt_prob_1 * gt_prob_2)
})
