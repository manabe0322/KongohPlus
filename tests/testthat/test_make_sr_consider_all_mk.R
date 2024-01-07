test_that("make_sr_consider_all_mk", {

  # Conditions
  markers_auto <- c("A", "B", "C")
  dt_parameter <- data.table(Marker = rep(c("B", "C", "A"), 4),
                             Parameter = c(rep("BSR", 3), rep("FSR", 3), rep("DSR", 3), rep("M2SR", 3)),
                             Model = c("LUS", "LUS", "LUS",
                                       "Not consider", "LUS", "LUS",
                                       "Not consider", "Not consider", "LUS",
                                       "Not consider", "Not consider", "Not consider"))

  # Run
  sr_consider_all_mk <- make_sr_consider_all_mk(markers_auto, dt_parameter)

  # Test
  expect_equal(nrow(sr_consider_all_mk), 3)
  expect_equal(ncol(sr_consider_all_mk), 4)
  expect_equal(sr_consider_all_mk[1, 1], TRUE)
  expect_equal(sr_consider_all_mk[2, 1], TRUE)
  expect_equal(sr_consider_all_mk[3, 1], TRUE)
  expect_equal(sr_consider_all_mk[1, 2], TRUE)
  expect_equal(sr_consider_all_mk[2, 2], FALSE)
  expect_equal(sr_consider_all_mk[3, 2], TRUE)
  expect_equal(sr_consider_all_mk[1, 3], TRUE)
  expect_equal(sr_consider_all_mk[2, 3], FALSE)
  expect_equal(sr_consider_all_mk[3, 3], FALSE)
  expect_equal(sr_consider_all_mk[1, 4], FALSE)
  expect_equal(sr_consider_all_mk[2, 4], FALSE)
  expect_equal(sr_consider_all_mk[3, 4], FALSE)
})
