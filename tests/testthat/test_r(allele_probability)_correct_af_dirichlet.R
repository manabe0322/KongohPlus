test_that("correct_af_dirichlet, no unobserved allele", {
  # Condition
  pop_al <- 12:21
  pop_freq <- c(6, 3, 79, 1192, 910, 608, 194, 8, 1, 1)
  unobs_al <- numeric(0)

  # Run
  tmp <- correct_af_dirichlet(pop_al, pop_freq, unobs_al, add_unobs = TRUE)
  pop_al <- tmp$pop_al
  pop_prob <- tmp$pop_prob

  # Test
  expect_equal(length(pop_al), 11)
  expect_equal(pop_al[1], 0)
  expect_equal(pop_al[2], 12)
  expect_equal(pop_al[3], 13)
  expect_equal(pop_al[4], 14)
  expect_equal(pop_al[5], 15)
  expect_equal(pop_al[6], 16)
  expect_equal(pop_al[7], 17)
  expect_equal(pop_al[8], 18)
  expect_equal(pop_al[9], 19)
  expect_equal(pop_al[10], 20)
  expect_equal(pop_al[11], 21)
  expect_equal(pop_prob[1], 1 / 3013)
  expect_equal(pop_prob[2], 7 / 3013)
  expect_equal(pop_prob[3], 4 / 3013)
  expect_equal(pop_prob[4], 80 / 3013)
  expect_equal(pop_prob[5], 1193 / 3013)
  expect_equal(pop_prob[6], 911 / 3013)
  expect_equal(pop_prob[7], 609 / 3013)
  expect_equal(pop_prob[8], 195 / 3013)
  expect_equal(pop_prob[9], 9 / 3013)
  expect_equal(pop_prob[10], 2 / 3013)
  expect_equal(pop_prob[11], 2 / 3013)
})

test_that("correct_af_dirichlet, one unobserved allele", {
  # Condition
  pop_al <- 12:21
  pop_freq <- c(6, 3, 79, 1192, 910, 608, 194, 8, 1, 1)
  unobs_al <- 15.2

  # Run
  tmp <- correct_af_dirichlet(pop_al, pop_freq, unobs_al, add_unobs = TRUE)
  pop_al <- tmp$pop_al
  pop_prob <- tmp$pop_prob

  # Test
  expect_equal(length(pop_al), 11)
  expect_equal(pop_al[1], 12)
  expect_equal(pop_al[2], 13)
  expect_equal(pop_al[3], 14)
  expect_equal(pop_al[4], 15)
  expect_equal(pop_al[5], 15.2)
  expect_equal(pop_al[6], 16)
  expect_equal(pop_al[7], 17)
  expect_equal(pop_al[8], 18)
  expect_equal(pop_al[9], 19)
  expect_equal(pop_al[10], 20)
  expect_equal(pop_al[11], 21)
  expect_equal(pop_prob[1], 7 / 3013)
  expect_equal(pop_prob[2], 4 / 3013)
  expect_equal(pop_prob[3], 80 / 3013)
  expect_equal(pop_prob[4], 1193 / 3013)
  expect_equal(pop_prob[5], 1 / 3013)
  expect_equal(pop_prob[6], 911 / 3013)
  expect_equal(pop_prob[7], 609 / 3013)
  expect_equal(pop_prob[8], 195 / 3013)
  expect_equal(pop_prob[9], 9 / 3013)
  expect_equal(pop_prob[10], 2 / 3013)
  expect_equal(pop_prob[11], 2 / 3013)
})

test_that("correct_af_dirichlet, two unobserved alleles", {
  # Condition
  pop_al <- 12:21
  pop_freq <- c(6, 3, 79, 1192, 910, 608, 194, 8, 1, 1)
  unobs_al <- c(15.2, 22)

  # Run
  tmp <- correct_af_dirichlet(pop_al, pop_freq, unobs_al, add_unobs = TRUE)
  pop_al <- tmp$pop_al
  pop_prob <- tmp$pop_prob

  # Test
  expect_equal(length(pop_al), 12)
  expect_equal(pop_al[1], 12)
  expect_equal(pop_al[2], 13)
  expect_equal(pop_al[3], 14)
  expect_equal(pop_al[4], 15)
  expect_equal(pop_al[5], 15.2)
  expect_equal(pop_al[6], 16)
  expect_equal(pop_al[7], 17)
  expect_equal(pop_al[8], 18)
  expect_equal(pop_al[9], 19)
  expect_equal(pop_al[10], 20)
  expect_equal(pop_al[11], 21)
  expect_equal(pop_al[12], 22)
  expect_equal(pop_prob[1], 7 / 3014)
  expect_equal(pop_prob[2], 4 / 3014)
  expect_equal(pop_prob[3], 80 / 3014)
  expect_equal(pop_prob[4], 1193 / 3014)
  expect_equal(pop_prob[5], 1 / 3014)
  expect_equal(pop_prob[6], 911 / 3014)
  expect_equal(pop_prob[7], 609 / 3014)
  expect_equal(pop_prob[8], 195 / 3014)
  expect_equal(pop_prob[9], 9 / 3014)
  expect_equal(pop_prob[10], 2 / 3014)
  expect_equal(pop_prob[11], 2 / 3014)
  expect_equal(pop_prob[12], 1 / 3014)
})

test_that("correct_af_dirichlet, no unobserved allele, add_unobs = FALSE", {
  # Condition
  pop_al <- 12:21
  pop_freq <- c(6, 3, 79, 1192, 910, 608, 194, 8, 1, 1)
  unobs_al <- numeric(0)

  # Run
  tmp <- correct_af_dirichlet(pop_al, pop_freq, unobs_al, add_unobs = FALSE)
  pop_al <- tmp$pop_al
  pop_prob <- tmp$pop_prob

  # Test
  expect_equal(length(pop_al), 10)
  expect_equal(pop_al[1], 12)
  expect_equal(pop_al[2], 13)
  expect_equal(pop_al[3], 14)
  expect_equal(pop_al[4], 15)
  expect_equal(pop_al[5], 16)
  expect_equal(pop_al[6], 17)
  expect_equal(pop_al[7], 18)
  expect_equal(pop_al[8], 19)
  expect_equal(pop_al[9], 20)
  expect_equal(pop_al[10], 21)
  expect_equal(pop_prob[1], 7 / 3012)
  expect_equal(pop_prob[2], 4 / 3012)
  expect_equal(pop_prob[3], 80 / 3012)
  expect_equal(pop_prob[4], 1193 / 3012)
  expect_equal(pop_prob[5], 911 / 3012)
  expect_equal(pop_prob[6], 609 / 3012)
  expect_equal(pop_prob[7], 195 / 3012)
  expect_equal(pop_prob[8], 9 / 3012)
  expect_equal(pop_prob[9], 2 / 3012)
  expect_equal(pop_prob[10], 2 / 3012)
})
