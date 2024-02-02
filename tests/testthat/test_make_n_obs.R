test_that("make_n_obs, pattern 1", {
  # Conditions
  pop_al <- c(16, 17, 18, 19, 99)
  prof_ref <- c(16, 17, 18, 19)

  # Run
  n_obs <- make_n_obs(pop_al, prof_ref)

  # Test
  expect_equal(length(n_obs), 5)
  expect_equal(n_obs[1], 1)
  expect_equal(n_obs[2], 1)
  expect_equal(n_obs[3], 1)
  expect_equal(n_obs[4], 1)
  expect_equal(n_obs[5], 0)
})

test_that("make_n_obs, pattern 2", {
  # Conditions
  pop_al <- c(16, 17, 18, 19, 99)
  prof_ref <- c(16, 17, 16, 17)

  # Run
  n_obs <- make_n_obs(pop_al, prof_ref)

  # Test
  expect_equal(length(n_obs), 5)
  expect_equal(n_obs[1], 2)
  expect_equal(n_obs[2], 2)
  expect_equal(n_obs[3], 0)
  expect_equal(n_obs[4], 0)
  expect_equal(n_obs[5], 0)
})

test_that("make_n_obs, pattern 3", {
  # Conditions
  pop_al <- c(16, 17, 18, 19, 99)
  prof_ref <- c(16, 17, 15, 20)

  # Run
  n_obs <- make_n_obs(pop_al, prof_ref)

  # Test
  expect_equal(length(n_obs), 5)
  expect_equal(n_obs[1], 1)
  expect_equal(n_obs[2], 1)
  expect_equal(n_obs[3], 0)
  expect_equal(n_obs[4], 0)
  expect_equal(n_obs[5], 2)
})
