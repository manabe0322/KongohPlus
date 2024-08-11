test_that("pick_rep with repeat motif", {
  # Condition
  seq <- "TCTA [TCTG]2 [TCTA]8"

  # Run
  reps <- pick_rep(seq)

  # Test
  expect_equal(length(reps), 2)
  expect_equal(reps[1], 2)
  expect_equal(reps[2], 8)
})

test_that("pick_rep without repeat motif", {
  # Condition
  seq <- "TCTA TCTG TCTA"

  # Run
  reps <- pick_rep(seq)

  # Test
  expect_equal(reps, integer(0))
})
