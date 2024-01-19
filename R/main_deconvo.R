perform_deconvo <- function(dt_csp, dt_af, dt_analysis_method, dt_parameter, dt_alcor, dt_kit, dt_at, hnc,
                            mp_cand_c1 = c(0.0025, 0.005, 0.0075, 0.01, 0.025, 0.05, 0.075, 0.1, 0.2, 0.3, 0.4, 0.5)){
  dt_kit_auto <- extract_kit_auto(dt_kit)
  markers_auto <- unique(dt_kit_auto[, Marker])
  n_mk <- length(markers_auto)

  tmp <- make_peak_data(dt_csp, dt_at, dt_kit, markers_auto)
  peaks_list <- tmp$peaks_list
  sizes_list <- tmp$sizes_list
  heights_list <- tmp$heights_list
  n_sn <- length(peaks_list)

  tmp <- set_pop_prob(markers_auto, dt_af, peaks_list)
  pop_al_list <- tmp$pop_al_list
  pop_prob_list <- tmp$pop_prob_list

  sr_consider_all_mk <- make_sr_consider_all_mk(markers_auto, dt_parameter)

  m2sr_all_mk <- extract_m2sr_all_mk(markers_auto, dt_parameter)
  beta <- extract_gamma_beta(dt_parameter)

  mp_cand <- set_mp(hnc, mp_cand_c1)
  mp_cand <- asplit(mp_cand, 1)

  setkey(dt_alcor, Marker)
  setkey(dt_at, Marker)

  gtc_list <- gtc_prob_list <- product_list <- list()
  for(i in 1:n_sn){
    peaks_list_sn1 <- peaks_list[[i]]
    heights_list_sn1 <- heights_list[[i]]

    gtc_list_sn1 <- gtc_prob_list_sn1 <- product_list_sn1 <- list()
    for(j in 1:n_mk){
      mk <- markers_auto[j]
      dt_alcor_mk <- dt_alcor[mk]
      dt_at_mk <- dt_at[mk]
      at <- as.integer(dt_at_mk[, AT])

      peaks <- peaks_list_sn1[[mk]]
      oph <- heights_list_sn1[[mk]]

      pop_al <- pop_al_list[[mk]]
      pop_prob <- pop_prob_list[[mk]]

      sr_consider <- sr_consider_all_mk[mk, ]

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

      # All candidate genotypes
      cand_gts <- combinations(n = length(peaks), r = 2, v = peaks, repeats.allowed = TRUE)
      cand_gts <- asplit(cand_gts, 1)

      # Indices of genotypes for candidate genotype combinations of each contributor
      index_gt_comb <- permutations(n = length(cand_gts), r = hnc, repeats.allowed = TRUE)
      index_gt_comb <- index_gt_comb - 1 # for C++ index
      index_gt_comb <- asplit(index_gt_comb, 1)

      tmp <- perform_discrete_analysis(peaks, cand_gts, index_gt_comb, pop_al, pop_prob, sr_consider, bsrs, fsrs, dsrs, m2sr, mp_cand, oph, beta, at)
      gtc_list_sn1[[j]] <- tmp[[1]]
      gtc_prob_list_sn1[[j]] <- tmp[[2]]
      product_list_sn1[[j]] <- tmp[[3]]
    }
    gtc_list[[i]] <- gtc_list_sn1
    gtc_prob_list[[i]] <- gtc_prob_list_sn1
    product_list[[i]] <- product_list_sn1
  }

  save(list = c("gtc_list", "gtc_prob_list", "product_list", "pop_al_list", "pop_prob_list"), file = "test.Rdata")
}
