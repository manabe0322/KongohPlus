perform_deconvo <- function(dt_csp, dt_af, dt_analysis_method, dt_parameter, dt_alcor, dt_kit, dt_at, hnc){
  dt_kit_auto <- extract_kit_auto(dt_kit)
  markers_auto <- unique(dt_kit_auto[, Marker])

  tmp <- make_peak_data(dt_csp, dt_at, dt_kit, markers_auto)
  peaks_list <- tmp$peaks_list
  sizes_list <- tmp$sizes_list
  heights_list <- tmp$heights_list

  sr_consider_all_mk <- make_sr_consider_all_mk(markers_auto, dt_parameter)

  m2sr_all_mk <- extract_m2sr_all_mk(markers_auto, dt_parameter)

  eph_base_all <- calc_eph_base_all(hnc, markers_auto, peaks_list, sr_consider_all_mk, dt_alcor, m2sr_all_mk)

  print(sapply(eph_base_all[[1]], nrow))
}
