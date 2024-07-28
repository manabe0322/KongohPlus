test_that("calc_fs_one pattern 1", {
  # Condition
  al <- 15

  # Run
  fs <- calc_fs_one(al)

  # Test
  expect_equal(length(fs), 1)
  expect_equal(fs, 16)
})

test_that("calc_fs_one pattern 2", {
  # Condition
  al <- 15.1

  # Run
  fs <- calc_fs_one(al)

  # Test
  expect_equal(length(fs), 1)
  expect_equal(fs, 16.1)
})

test_that("calc_fs_one pattern 3", {
  # Condition
  al <- 15.2

  # Run
  fs <- calc_fs_one(al)

  # Test
  expect_equal(length(fs), 1)
  expect_equal(fs, 16.2)
})

test_that("calc_fs_one pattern 4", {
  # Condition
  al <- 15.3

  # Run
  fs <- calc_fs_one(al)

  # Test
  expect_equal(length(fs), 1)
  expect_equal(fs, 16.3)
})
