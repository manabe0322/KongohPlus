test_that("remove_pu pattern1", {

  ##############
  # Conditions #
  ##############

  # Information on peaks
  peaks <- as.character(c(11, 12, 13, 14))
  sizes <- c(100, 104, 108, 112)
  heights <- c(200, 1000, 100, 1005)

  # Genotype
  gt <- as.character(c(12, 14))

  # Dye
  dye <- "B"

  # Information on alleles
  al_size <- c(104, 112, 200)
  al_dye <- c("B", "B", "G")

  # The maximum distance between a pull-up peak and the parent allele
  max_dist_pu <- 1

  #######
  # Run #
  #######

  tmp <- remove_pu(peaks, sizes, heights, gt, dye, al_size, al_dye, max_dist_pu)
  peaks_new <- tmp[[1]]
  sizes_new <- tmp[[2]]
  heights_new <- tmp[[3]]

  ########
  # Test #
  ########

  expect_equal(length(peaks_new), 4)
  expect_equal(peaks_new[1], "11")
  expect_equal(peaks_new[2], "12")
  expect_equal(peaks_new[3], "13")
  expect_equal(peaks_new[4], "14")
  expect_equal(sizes_new[1], 100)
  expect_equal(sizes_new[2], 104)
  expect_equal(sizes_new[3], 108)
  expect_equal(sizes_new[4], 112)
  expect_equal(heights_new[1], 200)
  expect_equal(heights_new[2], 1000)
  expect_equal(heights_new[3], 100)
  expect_equal(heights_new[4], 1005)
})

test_that("remove_pu pattern2", {

  ##############
  # Conditions #
  ##############

  # Information on peaks
  peaks <- as.character(c(11, 12, 13, 14))
  sizes <- c(100, 104, 108, 112)
  heights <- c(200, 1000, 100, 1005)

  # Genotype
  gt <- as.character(c(12, 14))

  # Dye
  dye <- "B"

  # Information on alleles
  al_size <- c(104, 112, 101)
  al_dye <- c("B", "B", "G")

  # The maximum distance between a pull-up peak and the parent allele
  max_dist_pu <- 1

  #######
  # Run #
  #######

  tmp <- remove_pu(peaks, sizes, heights, gt, dye, al_size, al_dye, max_dist_pu)
  peaks_new <- tmp[[1]]
  sizes_new <- tmp[[2]]
  heights_new <- tmp[[3]]

  ########
  # Test #
  ########

  expect_equal(length(peaks_new), 3)
  expect_equal(peaks_new[1], "12")
  expect_equal(peaks_new[2], "13")
  expect_equal(peaks_new[3], "14")
  expect_equal(sizes_new[1], 104)
  expect_equal(sizes_new[2], 108)
  expect_equal(sizes_new[3], 112)
  expect_equal(heights_new[1], 1000)
  expect_equal(heights_new[2], 100)
  expect_equal(heights_new[3], 1005)
})

test_that("remove_pu pattern3", {

  ##############
  # Conditions #
  ##############

  # Information on peaks
  peaks <- as.character(c(11, 12, 13, 14))
  sizes <- c(100, 104, 108, 112)
  heights <- c(200, 1000, 100, 1005)

  # Genotype
  gt <- as.character(c(12, 14))

  # Dye
  dye <- "B"

  # Information on alleles
  al_size <- c(104, 112, 103)
  al_dye <- c("B", "B", "G")

  # The maximum distance between a pull-up peak and the parent allele
  max_dist_pu <- 1

  #######
  # Run #
  #######

  tmp <- remove_pu(peaks, sizes, heights, gt, dye, al_size, al_dye, max_dist_pu)
  peaks_new <- tmp[[1]]
  sizes_new <- tmp[[2]]
  heights_new <- tmp[[3]]

  ########
  # Test #
  ########

  expect_equal(length(peaks_new), 4)
  expect_equal(peaks_new[1], "11")
  expect_equal(peaks_new[2], "12")
  expect_equal(peaks_new[3], "13")
  expect_equal(peaks_new[4], "14")
  expect_equal(sizes_new[1], 100)
  expect_equal(sizes_new[2], 104)
  expect_equal(sizes_new[3], 108)
  expect_equal(sizes_new[4], 112)
  expect_equal(heights_new[1], 200)
  expect_equal(heights_new[2], 1000)
  expect_equal(heights_new[3], 100)
  expect_equal(heights_new[4], 1005)
})

test_that("remove_pu pattern4", {

  ##############
  # Conditions #
  ##############

  # Information on peaks
  peaks <- as.character(c(11, 12, 13, 14))
  sizes <- c(100, 104, 108, 112)
  heights <- c(200, 1000, 100, 1005)

  # Genotype
  gt <- as.character(c(12, 14))

  # Dye
  dye <- "B"

  # Information on alleles
  al_size <- c(104, 112, 101, 109)
  al_dye <- c("B", "B", "G", "G")

  # The maximum distance between a pull-up peak and the parent allele
  max_dist_pu <- 1

  #######
  # Run #
  #######

  tmp <- remove_pu(peaks, sizes, heights, gt, dye, al_size, al_dye, max_dist_pu)
  peaks_new <- tmp[[1]]
  sizes_new <- tmp[[2]]
  heights_new <- tmp[[3]]

  ########
  # Test #
  ########

  expect_equal(length(peaks_new), 2)
  expect_equal(peaks_new[1], "12")
  expect_equal(peaks_new[2], "14")
  expect_equal(sizes_new[1], 104)
  expect_equal(sizes_new[2], 112)
  expect_equal(heights_new[1], 1000)
  expect_equal(heights_new[2], 1005)
})
