test_that("check_mk_est_par, Condition 1", {
  # Condition
  dt <- data.table(SampleName = c(rep("Sample_A", 3), rep("Sample_B", 3), rep("Sample_C", 3)),
                   Marker = rep(c("Marker_A", "Marker_B", "Marker_C"), 3))

  # Run
  sn_mk_problem <- check_mk_est_par(dt)

  # Test
  expect_equal(length(sn_mk_problem), 0)
})

test_that("check_mk_est_par, Condition 2", {
  # Condition
  dt <- data.table(SampleName = c(rep("Sample_A", 3), rep("Sample_B", 4), rep("Sample_C", 3)),
                   Marker = c("Marker_A", "Marker_B", "Marker_C", "Marker_A", "Marker_B", "Marker_C", "Marker_D", "Marker_A", "Marker_B", "Marker_C"))

  # Run
  sn_mk_problem <- check_mk_est_par(dt)

  # Test
  expect_equal(length(sn_mk_problem), 1)
  expect_equal(sn_mk_problem[1], "Sample_B")
})

test_that("check_mk_est_par, Condition 3", {
  # Condition
  dt <- data.table(SampleName = c(rep("Sample_A", 4), rep("Sample_B", 3), rep("Sample_C", 3)),
                   Marker = c("Marker_A", "Marker_B", "Marker_C", "Marker_D", "Marker_A", "Marker_B", "Marker_C", "Marker_A", "Marker_B", "Marker_C"))

  # Run
  sn_mk_problem <- check_mk_est_par(dt)

  # Test
  expect_equal(length(sn_mk_problem), 1)
  expect_equal(sn_mk_problem[1], "Sample_A")
})

test_that("check_mk_est_par, Condition 4", {
  # Condition
  dt <- data.table(SampleName = c(rep("Sample_A", 4), rep("Sample_B", 3), rep("Sample_C", 3), rep("Sample_D", 4)),
                   Marker = c("Marker_A", "Marker_B", "Marker_C", "Marker_D", "Marker_A", "Marker_B", "Marker_C", "Marker_A", "Marker_B", "Marker_C", "Marker_A", "Marker_B", "Marker_C", "Marker_D"))

  # Run
  sn_mk_problem <- check_mk_est_par(dt)

  # Test
  expect_equal(length(sn_mk_problem), 2)
  expect_equal(sn_mk_problem[1], "Sample_A")
  expect_equal(sn_mk_problem[2], "Sample_D")
})
