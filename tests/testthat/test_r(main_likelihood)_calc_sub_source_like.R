test_that("calc_sub_source_like (vs calc_subsub_source_like, hnc = 2)", {
  ## Number of markers: 3
  ## Two person mixture with one major and one minor (e.g., 10:1 mixture)
  ## Reference A is the actual major contributor.
  ## If Reference A is assumed to be the minor, log10 likelihoods in each marker are expected to be -Inf.
  ## In that case, the sub-sub source level LR is expected to be twice as large as the sub source level LR.

  # Condition
  hyp_id_all_hp <- make_hyp_id_all(2, "A", c("A", "B"))
  hyp_id_all_hd <- make_hyp_id_all(2, NULL, c("A", "B"))
  log10like_all_hyp_hp <- matrix(c(-10, -20, -30, -Inf, -Inf, -Inf), 3, 2) # n_mk = 3, n_hyp = 2
  log10like_all_hyp_hd <- matrix(c(-15, -30, -45), 3, 1) # n_mk = 3, n_hyp = 1

  # Run
  sub_source_like_hp <- calc_sub_source_like(hyp_id_all_hp, log10like_all_hyp_hp) # n_unk = 1 -> factorial(n_unk) = 1
  sub_source_like_hd <- calc_sub_source_like(hyp_id_all_hd, log10like_all_hyp_hd) # n_unk = 2 -> factorial(n_unk) = 2
  subsub_source_like_hp <- calc_subsub_source_like(log10like_all_hyp_hp)
  subsub_source_like_hd <- calc_subsub_source_like(log10like_all_hyp_hd)
  lr_sub <- 10^(sub_source_like_hp - sub_source_like_hd)
  lr_subsub <- 10^(subsub_source_like_hp - subsub_source_like_hd)

  # Test
  expect_equal(lr_subsub[1] / lr_sub[1], 2)
  expect_equal(lr_subsub[2] / lr_sub[2], 2)
  expect_equal(lr_subsub[3] / lr_sub[3], 2)
})

test_that("calc_sub_source_like (vs calc_subsub_source_like, hnc = 3)", {
  ## Number of markers: 3
  ## Three person mixture with one major and two minors (e.g., 10:1:1 mixture)
  ## Reference A is the actual major contributor.
  ## If Reference A is assumed to be a minor, log10 likelihoods in each marker are expected to be -Inf.
  ## In that case, the sub-sub source level LR is expected to be three times greater than the sub source level LR.

  # Condition
  hyp_id_all_hp <- make_hyp_id_all(3, "A", c("A", "B"))
  hyp_id_all_hd <- make_hyp_id_all(3, NULL, c("A", "B"))
  log10like_all_hyp_hp <- matrix(c(-10, -20, -30, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf), 3, 3) # n_mk = 3, n_hyp = 3, person A assumed to be a major contributor
  log10like_all_hyp_hd <- matrix(c(-15, -30, -45), 3, 1) # n_mk = 3, n_hyp = 1

  # Run
  sub_source_like_hp <- calc_sub_source_like(hyp_id_all_hp, log10like_all_hyp_hp) # n_unk = 2 -> factorial(n_unk) = 2
  sub_source_like_hd <- calc_sub_source_like(hyp_id_all_hd, log10like_all_hyp_hd) # n_unk = 3 -> factorial(n_unk) = 6
  subsub_source_like_hp <- calc_subsub_source_like(log10like_all_hyp_hp)
  subsub_source_like_hd <- calc_subsub_source_like(log10like_all_hyp_hd)
  lr_sub <- 10^(sub_source_like_hp - sub_source_like_hd)
  lr_subsub <- 10^(subsub_source_like_hp - subsub_source_like_hd)

  # Test
  expect_equal(lr_subsub[1] / lr_sub[1], 3)
  expect_equal(lr_subsub[2] / lr_sub[2], 3)
  expect_equal(lr_subsub[3] / lr_sub[3], 3)
})
