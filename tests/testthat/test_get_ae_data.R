test_that("get_ae_data", {

  ##############
  # Conditions #
  ##############

  # Peak heights of sample 1
  height_list_1 <- list()
  height_list_1[[1]] <- numeric(0)
  height_list_1[[2]] <- 200
  height_list_1[[3]] <- c(201, 199)
  height_list_1[[4]] <- c(302, 298)
  height_list_1[[5]] <- c(403, 397)
  names(height_list_1) <- paste0("mk", 1:5)

  # Peak heights of sample 2
  height_list_2 <- list()
  height_list_2[[1]] <- c(501, 499)
  height_list_2[[2]] <- 1000
  height_list_2[[3]] <- c(2010, 1990)
  height_list_2[[4]] <- c(3020, 2980)
  height_list_2[[5]] <- c(4030, 3970)
  names(height_list_2) <- paste0("mk", 1:5)

  height_list <- list()
  height_list[[1]] <- height_list_1
  height_list[[2]] <- height_list_2
  names(height_list) <- c("sample1", "sample2")

  #######
  # Run #
  #######

  ae_list <- get_ae_data(height_list)
  ae_data <- ae_list[[1]]
  height_sum <- ae_list[[2]]

  ########
  # Test #
  ########

  expect_equal(nrow(ae_data), 1)
  expect_equal(ae_data[1, 1], 0.25)
  expect_equal(ae_data[1, 2], 0.25)
  expect_equal(ae_data[1, 3], 1)
  expect_equal(ae_data[1, 4], 1.5)
  expect_equal(ae_data[1, 5], 2)
  expect_equal(nrow(height_sum), 1)
  expect_equal(height_sum[1, 1], 1000)
  expect_equal(height_sum[1, 2], 1000)
  expect_equal(height_sum[1, 3], 4000)
  expect_equal(height_sum[1, 4], 6000)
  expect_equal(height_sum[1, 5], 8000)
})
