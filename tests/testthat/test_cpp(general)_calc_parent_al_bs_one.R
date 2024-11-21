test_that("calc_parent_al_bs_one pattern 1", {
  # Condition
  bs <- 15

  # Run
  al <- calc_parent_al_bs_one(bs)

  # Test
  expect_equal(length(al), 1)
  expect_equal(al, 16)
})

test_that("calc_parent_al_bs_one pattern 2", {
  # Condition
  bs <- 15.1

  # Run
  al <- calc_parent_al_bs_one(bs)

  # Test
  expect_equal(length(al), 1)
  expect_equal(al, 16.1)
})

test_that("calc_parent_al_bs_one pattern 3", {
  # Condition
  bs <- 15.2

  # Run
  al <- calc_parent_al_bs_one(bs)

  # Test
  expect_equal(length(al), 1)
  expect_equal(al, 16.2)
})

test_that("calc_parent_al_bs_one pattern 4", {
  # Condition
  bs <- 15.3

  # Run
  al <- calc_parent_al_bs_one(bs)

  # Test
  expect_equal(length(al), 1)
  expect_equal(al, 16.3)
})
