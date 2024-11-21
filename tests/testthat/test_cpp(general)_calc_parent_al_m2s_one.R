test_that("calc_parent_al_m2s_one pattern 1", {
  # Condition
  m2s <- 15

  # Run
  al <- calc_parent_al_m2s_one(m2s)

  # Test
  expect_equal(length(al), 1)
  expect_equal(al, 15.2)
})

test_that("calc_parent_al_m2s_one pattern 2", {
  # Condition
  m2s <- 15.1

  # Run
  al <- calc_parent_al_m2s_one(m2s)

  # Test
  expect_equal(length(al), 1)
  expect_equal(al, 15.3)
})

test_that("calc_parent_al_m2s_one pattern 3", {
  # Condition
  m2s <- 15.2

  # Run
  al <- calc_parent_al_m2s_one(m2s)

  # Test
  expect_equal(length(al), 1)
  expect_equal(al, 16)
})

test_that("calc_parent_al_m2s_one pattern 4", {
  # Condition
  m2s <- 15.3

  # Run
  al <- calc_parent_al_m2s_one(m2s)

  # Test
  expect_equal(length(al), 1)
  expect_equal(al, 16.1)
})
