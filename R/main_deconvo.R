perform_deconvo <- function(dt_csp, dt_af, dt_analysis_method, dt_parameter, dt_alcor, dt_kit, dt_at, hnc){
  markers_auto <- extract_auto_mk(dt_kit)

  tmp <- make_peak_data(dt_csp, dt_at, markers_auto)
  peaks_list <- tmp$peaks_list
  sizes_list <- tmp$sizes_list
  heights_list <- tmp$heights_list

  sr_consider_all_mk <- make_sr_consider_all_mk(markers_auto, dt_parameter)

  gt_comb_all <- make_gt_comb_all(hnc, peaks_list, sr_consider_all_mk)
}
