test_that("calc_m2s_one pattern 1", {
  # Condition
  al <- 15

  # Run
  m2s <- calc_m2s_one(al)

  # Test
  expect_equal(length(m2s), 1)
  expect_equal(m2s, 14.2)
})

test_that("calc_m2s_one pattern 2", {
  # Condition
  al <- 15.1

  # Run
  m2s <- calc_m2s_one(al)

  # Test
  expect_equal(length(m2s), 1)
  expect_equal(m2s, 14.3)
})

test_that("calc_m2s_one pattern 3", {
  # Condition
  al <- 15.2

  # Run
  m2s <- calc_m2s_one(al)

  # Test
  expect_equal(length(m2s), 1)
  expect_equal(m2s, 15)
})

test_that("calc_m2s_one pattern 4", {
  # Condition
  al <- 15.3

  # Run
  m2s <- calc_m2s_one(al)

  # Test
  expect_equal(length(m2s), 1)
  expect_equal(m2s, 15.1)
})
