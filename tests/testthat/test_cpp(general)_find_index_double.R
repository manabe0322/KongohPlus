test_that("find_index_double", {
  # Condition
  elements <- 1:6
  vec <- c(3, 1, 5, 2, 4)

  # Run
  index1 <- find_index_double(elements[1], vec)
  index2 <- find_index_double(elements[2], vec)
  index3 <- find_index_double(elements[3], vec)
  index4 <- find_index_double(elements[4], vec)
  index5 <- find_index_double(elements[5], vec)
  index6 <- find_index_double(elements[6], vec)

  # Test
  expect_equal(index1, 1)
  expect_equal(index2, 3)
  expect_equal(index3, 0)
  expect_equal(index4, 4)
  expect_equal(index5, 2)
  expect_equal(index6, -1)
})
