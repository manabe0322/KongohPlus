test_that("make_hyp_id_all pattern 1", {
  # Condition
  hnc <- 2
  hyp <- NULL
  refs <- c("A", "B")

  # Run
  hyp_id_all <- make_hyp_id_all(hnc, hyp, refs)

  # Test
  expect_equal(nrow(hyp_id_all), 1)
  expect_equal(ncol(hyp_id_all), hnc)
  expect_equal(hyp_id_all[1, 1], 0)
  expect_equal(hyp_id_all[1, 2], 0)
})

test_that("make_hyp_id_all pattern 2", {
  # Condition
  hnc <- 2
  hyp <- "A"
  refs <- c("A", "B")

  # Run
  hyp_id_all <- make_hyp_id_all(hnc, hyp, refs)

  # Test
  expect_equal(nrow(hyp_id_all), 2)
  expect_equal(ncol(hyp_id_all), hnc)
  expect_equal(hyp_id_all[1, 1], 1)
  expect_equal(hyp_id_all[1, 2], 0)
  expect_equal(hyp_id_all[2, 1], 0)
  expect_equal(hyp_id_all[2, 2], 1)
})

test_that("make_hyp_id_all pattern 3", {
  # Condition
  hnc <- 2
  hyp <- "B"
  refs <- c("A", "B")

  # Run
  hyp_id_all <- make_hyp_id_all(hnc, hyp, refs)

  # Test
  expect_equal(nrow(hyp_id_all), 2)
  expect_equal(ncol(hyp_id_all), hnc)
  expect_equal(hyp_id_all[1, 1], 2)
  expect_equal(hyp_id_all[1, 2], 0)
  expect_equal(hyp_id_all[2, 1], 0)
  expect_equal(hyp_id_all[2, 2], 2)
})

test_that("make_hyp_id_all pattern 4", {
  # Condition
  hnc <- 2
  hyp <- c("A", "B")
  refs <- c("A", "B")

  # Run
  hyp_id_all <- make_hyp_id_all(hnc, hyp, refs)

  # Test
  expect_equal(nrow(hyp_id_all), 2)
  expect_equal(ncol(hyp_id_all), hnc)
  expect_equal(hyp_id_all[1, 1], 1)
  expect_equal(hyp_id_all[1, 2], 2)
  expect_equal(hyp_id_all[2, 1], 2)
  expect_equal(hyp_id_all[2, 2], 1)
})

test_that("make_hyp_id_all pattern 5", {
  # Condition
  hnc <- 3
  hyp <- NULL
  refs <- c("A", "B")

  # Run
  hyp_id_all <- make_hyp_id_all(hnc, hyp, refs)

  # Test
  expect_equal(nrow(hyp_id_all), 1)
  expect_equal(ncol(hyp_id_all), hnc)
  expect_equal(hyp_id_all[1, 1], 0)
  expect_equal(hyp_id_all[1, 2], 0)
  expect_equal(hyp_id_all[1, 3], 0)
})

test_that("make_hyp_id_all pattern 6", {
  # Condition
  hnc <- 3
  hyp <- "A"
  refs <- c("A", "B")

  # Run
  hyp_id_all <- make_hyp_id_all(hnc, hyp, refs)

  # Test
  expect_equal(nrow(hyp_id_all), 3)
  expect_equal(ncol(hyp_id_all), hnc)
  expect_equal(hyp_id_all[1, 1], 1)
  expect_equal(hyp_id_all[1, 2], 0)
  expect_equal(hyp_id_all[1, 3], 0)
  expect_equal(hyp_id_all[2, 1], 0)
  expect_equal(hyp_id_all[2, 2], 1)
  expect_equal(hyp_id_all[2, 3], 0)
  expect_equal(hyp_id_all[3, 1], 0)
  expect_equal(hyp_id_all[3, 2], 0)
  expect_equal(hyp_id_all[3, 3], 1)
})

test_that("make_hyp_id_all pattern 7", {
  # Condition
  hnc <- 3
  hyp <- "B"
  refs <- c("A", "B")

  # Run
  hyp_id_all <- make_hyp_id_all(hnc, hyp, refs)

  # Test
  expect_equal(nrow(hyp_id_all), 3)
  expect_equal(ncol(hyp_id_all), hnc)
  expect_equal(hyp_id_all[1, 1], 2)
  expect_equal(hyp_id_all[1, 2], 0)
  expect_equal(hyp_id_all[1, 3], 0)
  expect_equal(hyp_id_all[2, 1], 0)
  expect_equal(hyp_id_all[2, 2], 2)
  expect_equal(hyp_id_all[2, 3], 0)
  expect_equal(hyp_id_all[3, 1], 0)
  expect_equal(hyp_id_all[3, 2], 0)
  expect_equal(hyp_id_all[3, 3], 2)
})

test_that("make_hyp_id_all pattern 8", {
  # Condition
  hnc <- 3
  hyp <- c("A", "B")
  refs <- c("A", "B")

  # Run
  hyp_id_all <- make_hyp_id_all(hnc, hyp, refs)

  # Test
  expect_equal(nrow(hyp_id_all), 6)
  expect_equal(ncol(hyp_id_all), hnc)
  expect_equal(hyp_id_all[1, 1], 1)
  expect_equal(hyp_id_all[1, 2], 2)
  expect_equal(hyp_id_all[1, 3], 0)
  expect_equal(hyp_id_all[2, 1], 1)
  expect_equal(hyp_id_all[2, 2], 0)
  expect_equal(hyp_id_all[2, 3], 2)
  expect_equal(hyp_id_all[3, 1], 2)
  expect_equal(hyp_id_all[3, 2], 1)
  expect_equal(hyp_id_all[3, 3], 0)
  expect_equal(hyp_id_all[4, 1], 0)
  expect_equal(hyp_id_all[4, 2], 1)
  expect_equal(hyp_id_all[4, 3], 2)
  expect_equal(hyp_id_all[5, 1], 2)
  expect_equal(hyp_id_all[5, 2], 0)
  expect_equal(hyp_id_all[5, 3], 1)
  expect_equal(hyp_id_all[6, 1], 0)
  expect_equal(hyp_id_all[6, 2], 2)
  expect_equal(hyp_id_all[6, 3], 1)
})
