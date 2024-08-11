test_that("setdiff_double pattern 1", {
  # Condition
  v_al <- c(11, 12)
  r_al <- c(11, 12)

  # Run
  ans1 <- setdiff_double(v_al, r_al)
  ans2 <- setdiff_double(r_al, v_al)

  # Test
  expect_equal(ans1, numeric(0))
  expect_equal(ans2, numeric(0))
})

test_that("setdiff_double pattern 2", {
  # Condition
  v_al <- c(11, 12)
  r_al <- c(11)

  # Run
  ans1 <- setdiff_double(v_al, r_al)
  ans2 <- setdiff_double(r_al, v_al)

  # Test
  expect_equal(length(ans1), 1)
  expect_equal(ans1, 12)
  expect_equal(ans2, numeric(0))
})

test_that("setdiff_double pattern 3", {
  # Condition
  v_al <- c(11, 12)
  r_al <- c(10, 11)

  # Run
  ans1 <- setdiff_double(v_al, r_al)
  ans2 <- setdiff_double(r_al, v_al)

  # Test
  expect_equal(length(ans1), 1)
  expect_equal(ans1, 12)
  expect_equal(length(ans2), 1)
  expect_equal(ans2, 10)
})

test_that("setdiff_double pattern 3", {
  # Condition
  v_al <- c(11, 12, 13, 14)
  r_al <- c(10, 11)

  # Run
  ans1 <- setdiff_double(v_al, r_al)
  ans2 <- setdiff_double(r_al, v_al)

  # Test
  expect_equal(length(ans1), 3)
  expect_equal(ans1[1], 12)
  expect_equal(ans1[2], 13)
  expect_equal(ans1[3], 14)
  expect_equal(length(ans2), 1)
  expect_equal(ans2, 10)
})
