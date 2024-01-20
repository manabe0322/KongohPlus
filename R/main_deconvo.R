perform_deconvo <- function(dt_csp, dt_af, dt_analysis_method, dt_parameter, dt_alcor, dt_kit, dt_at, hnc){
  dt_kit_auto <- extract_kit_auto(dt_kit)
  markers_auto <- unique(dt_kit_auto[, Marker])
  n_mk <- length(markers_auto)

  tmp <- make_peak_data(dt_csp, dt_at, dt_kit, markers_auto)
  peaks_list <- tmp$peaks_list
  sizes_list <- tmp$sizes_list
  heights_list <- tmp$heights_list
  n_sn <- length(peaks_list)

  sr_consider_all_mk <- make_sr_consider_all_mk(markers_auto, dt_parameter)

  ae_all_mk <- extract_myu_all_mk("AE", markers_auto, dt_parameter)
  m2sr_all_mk <- extract_myu_all_mk("M2SR", markers_auto, dt_parameter)

  setkey(dt_alcor, Marker)

  eph_base_list <- list()
  for(i in 1:n_sn){
    peaks_list_sn1 <- peaks_list[[i]]
    heights_list_sn1 <- heights_list[[i]]

    eph_base_list_sn1 <- list()
    for(j in 1:n_mk){
      mk <- markers_auto[j]

      peaks <- peaks_list_sn1[[mk]]

      # All candidate genotypes in each contributor
      cand_gts <- combinations(n = length(peaks), r = 2, v = peaks, repeats.allowed = TRUE)
      cand_gts <- asplit(cand_gts, 1)

      sr_consider <- sr_consider_all_mk[mk, ]

      dt_alcor_mk <- dt_alcor[mk]
      alcor_alleles <- dt_alcor_mk[, Allele]
      alcor_bsrs <- dt_alcor_mk[, Average_BSR]
      alcor_fsrs <- dt_alcor_mk[, Average_FSR]
      alcor_dsrs <- dt_alcor_mk[, Average_DSR]
      index_alcor_peaks <- match(peaks, alcor_alleles)
      bsrs <- alcor_bsrs[index_alcor_peaks]
      bsrs[is.na(bsrs)] <- 0
      fsrs <- alcor_fsrs[index_alcor_peaks]
      fsrs[is.na(fsrs)] <- 0
      dsrs <- alcor_dsrs[index_alcor_peaks]
      dsrs[is.na(dsrs)] <- 0
      m2sr <- m2sr_all_mk[mk]

      ae <- ae_all_mk[mk]

      eph_base_list_sn1[[j]] <- make_eph_base(peaks, cand_gts, sr_consider, bsrs, fsrs, dsrs, m2sr, ae)
    }
    eph_base_list[[i]] <- eph_base_list_sn1
  }

  save(list = c("eph_base_list"), file = "test.Rdata")
}
