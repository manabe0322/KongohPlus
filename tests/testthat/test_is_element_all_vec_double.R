test_that("is_element_all_vec_double, ans = TRUE", {
  # Condition
  vec_a <- c(3, 1, 5)
  vec_b <- 1:5

  # Run
  ans <- is_element_all_vec_double(vec_a, vec_b)

  # Test
  expect_equal(ans, TRUE)
})

test_that("is_element_all_vec_double, ans = FALSE", {
  # Condition
  vec_a <- c(1, 6, 2)
  vec_b <- 1:5

  # Run
  ans <- is_element_all_vec_double(vec_a, vec_b)

  # Test
  expect_equal(ans, FALSE)
})
