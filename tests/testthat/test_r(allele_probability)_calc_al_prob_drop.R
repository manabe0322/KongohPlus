test_that("calc_al_prob_drop with drop-out allele", {
  # Condition
  pop_al <- 1:4
  pop_prob <- c(0.05, 0.25, 0.3, 0.4)
  al_csp <- c(2, 3)

  # Run
  pop_prob_drop <- calc_al_prob_drop(al_csp, pop_al, pop_prob)

  # Test
  expect_equal(pop_prob_drop, 0.45)
})

test_that("calc_al_prob_drop without drop-out allele", {
  # Condition
  pop_al <- 1:4
  pop_prob <- c(0.05, 0.25, 0.3, 0.4)
  al_csp <- 1:4

  # Run
  pop_prob_drop <- calc_al_prob_drop(al_csp, pop_al, pop_prob)

  # Test
  expect_equal(pop_prob_drop, 0)
})
