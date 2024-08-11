test_that("calc_bs_one pattern 1", {
  # Condition
  al <- 15

  # Run
  bs <- calc_bs_one(al)

  # Test
  expect_equal(length(bs), 1)
  expect_equal(bs, 14)
})

test_that("calc_bs_one pattern 2", {
  # Condition
  al <- 15.1

  # Run
  bs <- calc_bs_one(al)

  # Test
  expect_equal(length(bs), 1)
  expect_equal(bs, 14.1)
})

test_that("calc_bs_one pattern 3", {
  # Condition
  al <- 15.2

  # Run
  bs <- calc_bs_one(al)

  # Test
  expect_equal(length(bs), 1)
  expect_equal(bs, 14.2)
})

test_that("calc_bs_one pattern 4", {
  # Condition
  al <- 15.3

  # Run
  bs <- calc_bs_one(al)

  # Test
  expect_equal(length(bs), 1)
  expect_equal(bs, 14.3)
})
