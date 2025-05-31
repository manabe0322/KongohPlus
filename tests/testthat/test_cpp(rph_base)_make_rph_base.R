test_that("make_rph_base", {
  # Condition
  peaks <- c(11, 12, 13, 13.2, 99)
  bsrs <- c(0.09, 0.1, 0.11, 0.12, 0)
  fsrs <- c(0.019, 0.02, 0.021, 0.022, 0)
  dsrs <- c(0.009, 0.01, 0.011, 0.012, 0)
  m2srs <- c(0.014, 0.015, 0.016, 0.017, 0)
  ae <- 1.1

  # Run
  rph_base <- make_rph_base(peaks, bsrs, fsrs, dsrs, m2srs, ae)

  # Test
  expect_equal(length(rph_base), 27)

  ## Allele
  for(i in 1:5){
    expect_equal(rph_base[i], ae / (1 + bsrs[i] + fsrs[i] + dsrs[i] + m2srs[i]))
  }

  ## Back stutter
  expect_equal(rph_base[6], ae * bsrs[2] / (1 + bsrs[2] + fsrs[2] + dsrs[2] + m2srs[2])) # back stutter of allele 12
  expect_equal(rph_base[7], ae * bsrs[3] / (1 + bsrs[3] + fsrs[3] + dsrs[3] + m2srs[3])) # back stutter of allele 13
  expect_equal(rph_base[8], 0)
  expect_equal(rph_base[9], 0)
  expect_equal(rph_base[10], 0)

  ## Forward stutter
  expect_equal(rph_base[11], 0)
  expect_equal(rph_base[12], ae * fsrs[1] / (1 + bsrs[1] + fsrs[1] + dsrs[1] + m2srs[1])) # forward stutter of allele 11
  expect_equal(rph_base[13], ae * fsrs[2] / (1 + bsrs[2] + fsrs[2] + dsrs[2] + m2srs[2])) # forward stutter of allele 12
  expect_equal(rph_base[14], 0)
  expect_equal(rph_base[15], 0)

  ## Double-back stutter
  expect_equal(rph_base[16], ae * dsrs[3] / (1 + bsrs[3] + fsrs[3] + dsrs[3] + m2srs[3])) # double-back stutter of allele 13
  expect_equal(rph_base[17], 0)
  expect_equal(rph_base[18], 0)
  expect_equal(rph_base[19], 0)
  expect_equal(rph_base[20], 0)

  # Minus 2-nt stutter
  expect_equal(rph_base[21], 0)
  expect_equal(rph_base[22], 0)
  expect_equal(rph_base[23], ae * m2srs[4] / (1 + bsrs[4] + fsrs[4] + dsrs[4] + m2srs[4])) # minus 2-nt stutter of allele 13.2
  expect_equal(rph_base[24], 0)
  expect_equal(rph_base[25], 0)

  # Ignore
  expect_equal(rph_base[26], 0)

  # Drop-in
  expect_equal(rph_base[27], -1)
})
