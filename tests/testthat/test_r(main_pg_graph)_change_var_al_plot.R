test_that("change_var_al_plot, tri", {
  # Condition
  pk <- c(10, 10.1, 10.2)
  rep_unit_len <- 3

  # Run
  pk_new <- change_var_al_plot(pk, rep_unit_len)

  # Test
  expect_equal(length(pk_new), 3)
  expect_equal(pk_new[1], 10)
  expect_equal(pk_new[2], 10.333)
  expect_equal(pk_new[3], 10.667)
})

test_that("change_var_al_plot, tetra", {
  # Condition
  pk <- c(10, 10.1, 10.2, 10.3)
  rep_unit_len <- 4

  # Run
  pk_new <- change_var_al_plot(pk, rep_unit_len)

  # Test
  expect_equal(length(pk_new), 4)
  expect_equal(pk_new[1], 10)
  expect_equal(pk_new[2], 10.25)
  expect_equal(pk_new[3], 10.5)
  expect_equal(pk_new[4], 10.75)
})

test_that("change_var_al_plot, penta", {
  # Condition
  pk <- c(10, 10.1, 10.2, 10.3, 10.4)
  rep_unit_len <- 5

  # Run
  pk_new <- change_var_al_plot(pk, rep_unit_len)

  # Test
  expect_equal(length(pk_new), 5)
  expect_equal(pk_new[1], 10)
  expect_equal(pk_new[2], 10.2)
  expect_equal(pk_new[3], 10.4)
  expect_equal(pk_new[4], 10.6)
  expect_equal(pk_new[5], 10.8)
})
