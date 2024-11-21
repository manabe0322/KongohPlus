test_that("calc_parent_al_ds_one pattern 1", {
  # Condition
  ds <- 15

  # Run
  al <- calc_parent_al_ds_one(ds)

  # Test
  expect_equal(length(al), 1)
  expect_equal(al, 17)
})

test_that("calc_parent_al_ds_one pattern 2", {
  # Condition
  ds <- 15.1

  # Run
  al <- calc_parent_al_ds_one(ds)

  # Test
  expect_equal(length(al), 1)
  expect_equal(al, 17.1)
})

test_that("calc_parent_al_ds_one pattern 3", {
  # Condition
  ds <- 15.2

  # Run
  al <- calc_parent_al_ds_one(ds)

  # Test
  expect_equal(length(al), 1)
  expect_equal(al, 17.2)
})

test_that("calc_parent_al_ds_one pattern 4", {
  # Condition
  ds <- 15.3

  # Run
  al <- calc_parent_al_ds_one(ds)

  # Test
  expect_equal(length(al), 1)
  expect_equal(al, 17.3)
})
