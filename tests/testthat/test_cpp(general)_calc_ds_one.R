test_that("calc_ds_one pattern 1", {
  # Condition
  al <- 15

  # Run
  ds <- calc_ds_one(al)

  # Test
  expect_equal(length(ds), 1)
  expect_equal(ds, 13)
})

test_that("calc_ds_one pattern 2", {
  # Condition
  al <- 15.1

  # Run
  ds <- calc_ds_one(al)

  # Test
  expect_equal(length(ds), 1)
  expect_equal(ds, 13.1)
})

test_that("calc_ds_one pattern 3", {
  # Condition
  al <- 15.2

  # Run
  ds <- calc_ds_one(al)

  # Test
  expect_equal(length(ds), 1)
  expect_equal(ds, 13.2)
})

test_that("calc_ds_one pattern 4", {
  # Condition
  al <- 15.3

  # Run
  ds <- calc_ds_one(al)

  # Test
  expect_equal(length(ds), 1)
  expect_equal(ds, 13.3)
})
