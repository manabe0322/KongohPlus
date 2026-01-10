test_that("calc_eph (for other than deg)", {
  # Condition
  mk <- "D3S1358"

  gt_select_list <- list()
  gt_select_list[[1]] <- c(14, 15)
  gt_select_list[[2]] <- c(15, 16)
  gt_select_list[[3]] <- c(16, 99)
  gt_select_list[[4]] <- c(99, 99)

  temps <- c(10000, 1000, 100, 10)
  deg <- c(1, 1, 1, 1)
  rep_ae_one <- 2
  ae <- 1.1

  dt_average_sr <- data.table("Allele" = c("9", "10", "11", "12", "13", "14", "15", "15.2", "16", "17", "18", "19", "20"),
                              "Corrected allele (BSR)" = c(9, 10, 11, 12, 13, 14, 15, 15.2, 16, 17, 18, 19, 20),
                              "Average BSR" = c(0.02550202, 0.03267084, 0.03983965, 0.04700847, 0.05417728, 0.06134610, 0.06851491, 0.06994868, 0.07568373, 0.08285254, 0.09002136, 0.09719017, 0.10435899),
                              "Corrected allele (FSR)" = c(9, 10, 11, 12, 13, 14, 15, 15.2, 16, 17, 18, 19, 20),
                              "Average FSR" = rep(0.007929056, 13),
                              "Corrected allele (DSR)" = c(9, 10, 11, 12, 13, 14, 15, 15.2, 16, 17, 18, 19, 20),
                              "Average DSR" = rep(0.006014601, 13),
                              "Corrected allele (M2SR)" = rep(NA, 13),
                              "Average M2SR" = rep(0, 13))

  dt_size <- data.table("Marker" = rep("D3S1358", 13),
                        "Repeat Unit Length" = rep(4, 13),
                        "Chromosome" = rep("Autosomal", 13),
                        "Allele" = c("9", "10", "11", "12", "13", "14", "15", "15.2", "16", "17", "18", "19", "20"),
                        "Corrected Allele" = rep(NA, 13),
                        "Size" = c(96.5200, 100.6850, 104.8550, 108.8750, 113.1250, 117.1150, 121.0000, 123.1353, 125.1200, 129.2550, 133.2900, 137.2750, 141.5900))

  base_size <- 100

  # Run
  eph <- calc_eph(gt_select_list, temps, deg, rep_ae_one, ae, dt_average_sr, dt_size, base_size)

  # Test
  expect_equal(length(eph), 8)
  expect_equal(names(eph)[1], "12")
  expect_equal(names(eph)[2], "13")
  expect_equal(names(eph)[3], "14")
  expect_equal(names(eph)[4], "15")
  expect_equal(names(eph)[5], "16")
  expect_equal(names(eph)[6], "17")
  expect_equal(names(eph)[7], "993")
  expect_equal(names(eph)[8], "994")

  tmp_14 <- 1 / (1 + 0.06134610 + 0.007929056 + 0.006014601)
  tmp_15 <- 1 / (1 + 0.06851491 + 0.007929056 + 0.006014601)
  tmp_16 <- 1 / (1 + 0.07568373 + 0.007929056 + 0.006014601)
  expect_equal(round(as.numeric(eph[1]), 5), round(22000 * 0.006014601 * tmp_14, 5))
  expect_equal(round(as.numeric(eph[2]), 5), round(22000 * 0.06134610 * tmp_14 + (22000 + 2200) * 0.006014601 * tmp_15, 5))
  expect_equal(round(as.numeric(eph[3]), 5), round(22000 * tmp_14 + (22000 + 2200) * 0.06851491 * tmp_15 + (2200 + 220) * 0.006014601 * tmp_16, 5))
  expect_equal(round(as.numeric(eph[4]), 5), round(22000 * 0.007929056 * tmp_14 + (22000 + 2200) * tmp_15 + (2200 + 220) * 0.07568373 * tmp_16, 5))
  expect_equal(round(as.numeric(eph[5]), 5), round((22000 + 2200) * 0.007929056 * tmp_15 + (2200 + 220) * tmp_16, 5))
  expect_equal(round(as.numeric(eph[6]), 5), round((2200 + 220) * 0.007929056 * tmp_16, 5))
  expect_equal(round(as.numeric(eph[7]), 5), 220)
  expect_equal(round(as.numeric(eph[8]), 5), 22)
})

test_that("calc_eph (for deg)", {
  # For testing this function when considering deg = 0.1, I assumed that the fragment size increases by 100 bp for each increment in allele number, using allele 14 as the reference point (base_size = 600 bp).

  # Condition
  mk <- "D3S1358"

  gt_select_list <- list()
  gt_select_list[[1]] <- c(14, 15)
  gt_select_list[[2]] <- c(15, 16)
  gt_select_list[[3]] <- c(16, 99)
  gt_select_list[[4]] <- c(99, 99)

  temps <- c(10000, 1000, 100, 10)
  deg <- c(0.1, 0.1, 0.1, 0.1)
  rep_ae_one <- 2
  ae <- 1.1

  dt_average_sr <- data.table("Allele" = c("9", "10", "11", "12", "13", "14", "15", "15.2", "16", "17", "18", "19", "20"),
                              "Corrected allele (BSR)" = c(9, 10, 11, 12, 13, 14, 15, 15.2, 16, 17, 18, 19, 20),
                              "Average BSR" = c(0.02550202, 0.03267084, 0.03983965, 0.04700847, 0.05417728, 0.06134610, 0.06851491, 0.06994868, 0.07568373, 0.08285254, 0.09002136, 0.09719017, 0.10435899),
                              "Corrected allele (FSR)" = c(9, 10, 11, 12, 13, 14, 15, 15.2, 16, 17, 18, 19, 20),
                              "Average FSR" = rep(0.007929056, 13),
                              "Corrected allele (DSR)" = c(9, 10, 11, 12, 13, 14, 15, 15.2, 16, 17, 18, 19, 20),
                              "Average DSR" = rep(0.006014601, 13),
                              "Corrected allele (M2SR)" = rep(NA, 13),
                              "Average M2SR" = rep(0, 13))

  dt_size <- data.table("Marker" = rep("D3S1358", 13),
                        "Repeat Unit Length" = rep(4, 13),
                        "Chromosome" = rep("Autosomal", 13),
                        "Allele" = c("9", "10", "11", "12", "13", "14", "15", "15.2", "16", "17", "18", "19", "20"),
                        "Corrected Allele" = rep(NA, 13),
                        "Size" = c(100, 200, 300, 400, 500, 600, 700, 750, 800, 900, 1000, 1100, 1200))

  base_size <- 600

  # Run
  eph <- calc_eph(gt_select_list, temps, deg, rep_ae_one, ae, dt_average_sr, dt_size, base_size)

  # Test
  expect_equal(length(eph), 8)
  expect_equal(names(eph)[1], "12")
  expect_equal(names(eph)[2], "13")
  expect_equal(names(eph)[3], "14")
  expect_equal(names(eph)[4], "15")
  expect_equal(names(eph)[5], "16")
  expect_equal(names(eph)[6], "17")
  expect_equal(names(eph)[7], "993")
  expect_equal(names(eph)[8], "994")

  tmp_14 <- 1 / (1 + 0.06134610 + 0.007929056 + 0.006014601)
  tmp_15 <- 1 / (1 + 0.06851491 + 0.007929056 + 0.006014601)
  tmp_16 <- 1 / (1 + 0.07568373 + 0.007929056 + 0.006014601)
  expect_equal(round(as.numeric(eph[1]), 5), round(22000 * 0.006014601 * tmp_14, 5))
  expect_equal(round(as.numeric(eph[2]), 5), round(22000 * 0.06134610 * tmp_14 + (2200 + 220) * 0.006014601 * tmp_15, 5))
  expect_equal(round(as.numeric(eph[3]), 5), round(22000 * tmp_14 + (2200 + 220) * 0.06851491 * tmp_15 + (22 + 2.2) * 0.006014601 * tmp_16, 5))
  expect_equal(round(as.numeric(eph[4]), 5), round(22000 * 0.007929056 * tmp_14 + (2200 + 220) * tmp_15 + (22 + 2.2) * 0.07568373 * tmp_16, 5))
  expect_equal(round(as.numeric(eph[5]), 5), round((2200 + 220) * 0.007929056 * tmp_15 + (22 + 2.2) * tmp_16, 5))
  expect_equal(round(as.numeric(eph[6]), 5), round((22 + 2.2) * 0.007929056 * tmp_16, 5))
  expect_equal(round(as.numeric(eph[7]), 5), 22)
  expect_equal(round(as.numeric(eph[8]), 5), 2.2)
})
