test_that("calc_parent_al_fs_one pattern 1", {
  # Condition
  fs <- 15

  # Run
  al <- calc_parent_al_fs_one(fs)

  # Test
  expect_equal(length(al), 1)
  expect_equal(al, 14)
})

test_that("calc_parent_al_fs_one pattern 2", {
  # Condition
  fs <- 15.1

  # Run
  al <- calc_parent_al_fs_one(fs)

  # Test
  expect_equal(length(al), 1)
  expect_equal(al, 14.1)
})

test_that("calc_parent_al_fs_one pattern 3", {
  # Condition
  fs <- 15.2

  # Run
  al <- calc_parent_al_fs_one(fs)

  # Test
  expect_equal(length(al), 1)
  expect_equal(al, 14.2)
})

test_that("calc_parent_al_fs_one pattern 4", {
  # Condition
  fs <- 15.3

  # Run
  al <- calc_parent_al_fs_one(fs)

  # Test
  expect_equal(length(al), 1)
  expect_equal(al, 14.3)
})
