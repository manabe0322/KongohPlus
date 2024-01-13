test_that("set_mp, pattern 1", {
  # Condition
  mp_c1 <- c(0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5)

  # Run
  mp_2 <- set_mp(2, mp_c1)
  mp_3 <- set_mp(3, mp_c1)
  mp_4 <- set_mp(4, mp_c1)

  # Test
  expect_equal(nrow(mp_2), 9)
  expect_equal(nrow(mp_3), 34)
  expect_equal(nrow(mp_4), 97)
})

test_that("set_mp, pattern 2", {
  # Condition
  mp_c1 <- c(0.0025, 0.005, 0.0075, 0.01, 0.025, 0.05, 0.075, 0.1, 0.2, 0.3, 0.4, 0.5)

  # Run
  mp_2 <- set_mp(2, mp_c1)
  mp_3 <- set_mp(3, mp_c1)
  mp_4 <- set_mp(4, mp_c1)

  # Test
  expect_equal(nrow(mp_2), 12)
  expect_equal(nrow(mp_3), 64)
  expect_equal(nrow(mp_4), 254)
})
