test_that("pick_al_size, AT = 30", {

  ##############
  # Conditions #
  ##############

  # Experimental data
  peak_all_mk <- size_all_mk <- height_all_mk <- matrix("", 5, 20)
  peak_all_mk[1, ] <- c("OL", "OL", "OL", "OL", "12", "OL", "14", "15", "16", "17", "OL", "OL", "", "", "", "", "", NA, "", NA)
  peak_all_mk[2, ] <- c("X", "Y", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", NA, "", NA)
  peak_all_mk[3, ] <- c("9", "OL", "10", "11", "12", "OL", "16", "", "", "", "", "", "", "", "", "", "", NA, "", NA)
  peak_all_mk[4, ] <- c("10", "11", "12", "13", "14", "15", "16", "", "", "", "", "", "", "", "", "", "", NA, "", NA)
  peak_all_mk[5, ] <- c("11", "12", "13", "14", "15", "", "", "", "", "", "", "", "", "", "", "", "", NA, "", NA)
  size_all_mk[1, ] <- c("97.34", "98.63", "104.56", "106.03", "109.28", "114.16", "117.36", "121.37", "125.40", "129.40", "138.33", "146.66", NA, NA, NA, NA, NA, NA, NA, NA)
  size_all_mk[2, ] <- c("98.79", "104.89", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
  size_all_mk[3, ] <- c("80.79", "82.44", "84.97", "89.10", "93.25", "106.19", "109.20", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
  size_all_mk[4, ] <- c("94.45", "97.42", "100.41", "103.50", "106.52", "109.45", "112.46", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
  size_all_mk[5, ] <- c("97.99", "101.96", "106.11", "110.10", "114.16", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
  height_all_mk[1, ] <- c("18", "18", "20", "22", "18", "26", "42", "709", "10164", "50", "49", "27", NA, NA, NA, NA, NA, NA, NA, NA)
  height_all_mk[2, ] <- c("6814", "6947", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
  height_all_mk[3, ] <- c("250", "39", "5303", "4504", "27", "21", "21", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
  height_all_mk[4, ] <- c("113", "4383", "61", "20", "310", "4286", "172", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
  height_all_mk[5, ] <- c("28", "292", "4464", "396", "4899", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
  dt_exp_right <- as.data.frame(cbind(peak_all_mk, size_all_mk, height_all_mk))
  setDT(dt_exp_right)
  names(dt_exp_right) <- c(paste0("Allele", 1:20), paste0("Size", 1:20), paste0("Height", 1:20))
  dt_exp_left <- data.table(SampleName = rep("Example1", 5),
                            Marker = c("D3S1358", "AMEL", "D2S441", "D22S1045", "D10S1248"))
  dt_exp_one <- cbind(dt_exp_left, dt_exp_right)
  dt_exp_one <- change_data_type_exp(dt_exp_one)

  # Genotype data
  dt_gt_one <- data.table(SampleName = rep("Example1", 5),
                          Marker = c("D3S1358", "AMEL", "D2S441", "D22S1045", "D10S1248"),
                          Allele1 = c("16", "X", "10", "11", "13"),
                          Allele2 = c("16", "Y", "11", "15", "15"))

  # Analytical thresholds
  at_all_mk <- rep(30, 5)
  names(at_all_mk) <- c("D3S1358", "AMEL", "D2S441", "D22S1045", "D10S1248")

  # Dyes
  dye_all_mk <- c("B", "G", "Y", "R", "P")
  names(dye_all_mk) <- c("D3S1358", "AMEL", "D2S441", "D22S1045", "D10S1248")

  #######
  # Run #
  #######

  al_size_info <- pick_al_size(dt_exp_one, dt_gt_one, at_all_mk, dye_all_mk)

  ########
  # Test #
  ########

  expect_equal(ncol(al_size_info), 9)
  expect_equal(colnames(al_size_info)[1], "D3S1358")
  expect_equal(colnames(al_size_info)[2], "AMEL")
  expect_equal(colnames(al_size_info)[3], "AMEL")
  expect_equal(colnames(al_size_info)[4], "D2S441")
  expect_equal(colnames(al_size_info)[5], "D2S441")
  expect_equal(colnames(al_size_info)[6], "D22S1045")
  expect_equal(colnames(al_size_info)[7], "D22S1045")
  expect_equal(colnames(al_size_info)[8], "D10S1248")
  expect_equal(colnames(al_size_info)[9], "D10S1248")
  expect_equal(al_size_info[1, 1], "16")
  expect_equal(al_size_info[1, 2], "X")
  expect_equal(al_size_info[1, 3], "Y")
  expect_equal(al_size_info[1, 4], "10")
  expect_equal(al_size_info[1, 5], "11")
  expect_equal(al_size_info[1, 6], "11")
  expect_equal(al_size_info[1, 7], "15")
  expect_equal(al_size_info[1, 8], "13")
  expect_equal(al_size_info[1, 9], "15")
  expect_equal(as.numeric(al_size_info[2, 1]), 125.40)
  expect_equal(as.numeric(al_size_info[2, 2]), 98.79)
  expect_equal(as.numeric(al_size_info[2, 3]), 104.89)
  expect_equal(as.numeric(al_size_info[2, 4]), 84.97)
  expect_equal(as.numeric(al_size_info[2, 5]), 89.10)
  expect_equal(as.numeric(al_size_info[2, 6]), 97.42)
  expect_equal(as.numeric(al_size_info[2, 7]), 109.45)
  expect_equal(as.numeric(al_size_info[2, 8]), 106.11)
  expect_equal(as.numeric(al_size_info[2, 9]), 114.16)
  expect_equal(al_size_info[3, 1], "B")
  expect_equal(al_size_info[3, 2], "G")
  expect_equal(al_size_info[3, 3], "G")
  expect_equal(al_size_info[3, 4], "Y")
  expect_equal(al_size_info[3, 5], "Y")
  expect_equal(al_size_info[3, 6], "R")
  expect_equal(al_size_info[3, 7], "R")
  expect_equal(al_size_info[3, 8], "P")
  expect_equal(al_size_info[3, 9], "P")
})

test_that("pick_al_size, AT = 4500", {

  ##############
  # Conditions #
  ##############

  # Experimental data
  peak_all_mk <- size_all_mk <- height_all_mk <- matrix("", 5, 20)
  peak_all_mk[1, ] <- c("OL", "OL", "OL", "OL", "12", "OL", "14", "15", "16", "17", "OL", "OL", "", "", "", "", "", NA, "", NA)
  peak_all_mk[2, ] <- c("X", "Y", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", NA, "", NA)
  peak_all_mk[3, ] <- c("9", "OL", "10", "11", "12", "OL", "16", "", "", "", "", "", "", "", "", "", "", NA, "", NA)
  peak_all_mk[4, ] <- c("10", "11", "12", "13", "14", "15", "16", "", "", "", "", "", "", "", "", "", "", NA, "", NA)
  peak_all_mk[5, ] <- c("11", "12", "13", "14", "15", "", "", "", "", "", "", "", "", "", "", "", "", NA, "", NA)
  size_all_mk[1, ] <- c("97.34", "98.63", "104.56", "106.03", "109.28", "114.16", "117.36", "121.37", "125.40", "129.40", "138.33", "146.66", NA, NA, NA, NA, NA, NA, NA, NA)
  size_all_mk[2, ] <- c("98.79", "104.89", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
  size_all_mk[3, ] <- c("80.79", "82.44", "84.97", "89.10", "93.25", "106.19", "109.20", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
  size_all_mk[4, ] <- c("94.45", "97.42", "100.41", "103.50", "106.52", "109.45", "112.46", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
  size_all_mk[5, ] <- c("97.99", "101.96", "106.11", "110.10", "114.16", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
  height_all_mk[1, ] <- c("18", "18", "20", "22", "18", "26", "42", "709", "10164", "50", "49", "27", NA, NA, NA, NA, NA, NA, NA, NA)
  height_all_mk[2, ] <- c("6814", "6947", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
  height_all_mk[3, ] <- c("250", "39", "5303", "4504", "27", "21", "21", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
  height_all_mk[4, ] <- c("113", "4383", "61", "20", "310", "4286", "172", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
  height_all_mk[5, ] <- c("28", "292", "4464", "396", "4899", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
  dt_exp_right <- as.data.frame(cbind(peak_all_mk, size_all_mk, height_all_mk))
  setDT(dt_exp_right)
  names(dt_exp_right) <- c(paste0("Allele", 1:20), paste0("Size", 1:20), paste0("Height", 1:20))
  dt_exp_left <- data.table(SampleName = rep("Example1", 5),
                            Marker = c("D3S1358", "AMEL", "D2S441", "D22S1045", "D10S1248"))
  dt_exp_one <- cbind(dt_exp_left, dt_exp_right)
  dt_exp_one <- change_data_type_exp(dt_exp_one)

  # Genotype data
  dt_gt_one <- data.table(SampleName = rep("Example1", 5),
                          Marker = c("D3S1358", "AMEL", "D2S441", "D22S1045", "D10S1248"),
                          Allele1 = c("16", "X", "10", "11", "13"),
                          Allele2 = c("16", "Y", "11", "15", "15"))

  # Analytical thresholds
  at_all_mk <- rep(4500, 5)
  names(at_all_mk) <- c("D3S1358", "AMEL", "D2S441", "D22S1045", "D10S1248")

  # Dyes
  dye_all_mk <- c("B", "G", "Y", "R", "P")
  names(dye_all_mk) <- c("D3S1358", "AMEL", "D2S441", "D22S1045", "D10S1248")

  #######
  # Run #
  #######

  al_size_info <- pick_al_size(dt_exp_one, dt_gt_one, at_all_mk, dye_all_mk)

  ########
  # Test #
  ########

  expect_equal(ncol(al_size_info), 6)
  expect_equal(colnames(al_size_info)[1], "D3S1358")
  expect_equal(colnames(al_size_info)[2], "AMEL")
  expect_equal(colnames(al_size_info)[3], "AMEL")
  expect_equal(colnames(al_size_info)[4], "D2S441")
  expect_equal(colnames(al_size_info)[5], "D2S441")
  expect_equal(colnames(al_size_info)[6], "D10S1248")
  expect_equal(al_size_info[1, 1], "16")
  expect_equal(al_size_info[1, 2], "X")
  expect_equal(al_size_info[1, 3], "Y")
  expect_equal(al_size_info[1, 4], "10")
  expect_equal(al_size_info[1, 5], "11")
  expect_equal(al_size_info[1, 6], "15")
  expect_equal(as.numeric(al_size_info[2, 1]), 125.40)
  expect_equal(as.numeric(al_size_info[2, 2]), 98.79)
  expect_equal(as.numeric(al_size_info[2, 3]), 104.89)
  expect_equal(as.numeric(al_size_info[2, 4]), 84.97)
  expect_equal(as.numeric(al_size_info[2, 5]), 89.10)
  expect_equal(as.numeric(al_size_info[2, 6]), 114.16)
  expect_equal(al_size_info[3, 1], "B")
  expect_equal(al_size_info[3, 2], "G")
  expect_equal(al_size_info[3, 3], "G")
  expect_equal(al_size_info[3, 4], "Y")
  expect_equal(al_size_info[3, 5], "Y")
  expect_equal(al_size_info[3, 6], "P")
})
