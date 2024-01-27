perform_deconvo <- function(path_pack, dt_csp, dt_af, dt_parameter, dt_alcor, dt_kit, dt_at, bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr, phr_fltr, st, hnc){
  dt_kit_auto <- extract_kit_auto(dt_kit)
  markers_auto <- unique(dt_kit_auto[, Marker])
  n_mk <- length(markers_auto)

  tmp <- make_peak_data(dt_csp, dt_at, dt_kit, markers_auto)
  peaks_list <- tmp$peaks_list
  sizes_list <- tmp$sizes_list
  heights_list <- tmp$heights_list
  n_csp <- length(peaks_list)

  sr_consider_all_mk <- make_sr_consider_all_mk(markers_auto, dt_parameter)

  ae_all_mk <- extract_myu_all_mk("AE", markers_auto, dt_parameter)
  m2sr_all_mk <- extract_myu_all_mk("M2SR", markers_auto, dt_parameter)

  setkey(dt_alcor, Marker)

  eph_base_all <- list()

  size_eph_all <- numeric(0)
  index_csp <- index_mk <- index_gt <- index_peak <- integer(0)
  n_csp_mk <- n_csp * n_mk
  n_gt_all <- n_peak_all <- n_gt_comb_all <- rep(0, n_csp_mk)
  index_gt_comb_all <- NULL
  cand_gts_all <- list()
  for(i in 1:n_csp){
    peaks_list_sn1 <- peaks_list[[i]]
    sizes_list_sn1 <- sizes_list[[i]]
    heights_list_sn1 <- heights_list[[i]]

    for(j in 1:n_mk){
      mk <- markers_auto[j]

      peaks <- peaks_list_sn1[[mk]]
      sizes <- sizes_list_sn1[[mk]]
      heights <- heights_list_sn1[[mk]]

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

      eph_base <- make_eph_base(peaks, bsrs, fsrs, dsrs, m2sr, ae)
      eph_base_all[[n_mk * (i - 1) + j]] <- eph_base


      # All candidate genotypes in each contributor
#      cand_gts <- combinations(n = length(peaks), r = 2, v = peaks, repeats.allowed = TRUE)
#      cand_gts_all[[n_mk * (i - 1) + j]] <- cand_gts
#      cand_gts <- asplit(cand_gts, 1)

#      sr_consider <- sr_consider_all_mk[mk, ]

#      n_eph_base <- length(eph_base)
#      n_gt <- length(cand_gts)
#      n_peak <- length(peaks)
#      index_csp <- c(index_csp, rep(i, n_eph_base))
#      index_mk <- c(index_mk, rep(j, n_eph_base))
#      index_gt <- c(index_gt, as.vector(sapply(1:n_gt, rep, n_peak)))
#      index_peak <- c(index_peak, rep(1:n_peak, n_gt))

#      n_gt_all[n_mk * (i - 1) + j] <- n_gt
#      n_peak_all[n_mk * (i - 1) + j] <- n_peak

#      size_eph_all <- c(size_eph_all, rep(sizes, n_gt))


      # Indices of genotypes for candidate genotype combinations of each contributor
#      index_gt_comb <- permutations(n = length(cand_gts), r = hnc, repeats.allowed = TRUE)
#      index_gt_comb <- index_gt_comb - 1 # for C++ index
#      index_gt_comb <- asplit(index_gt_comb, 1)

#      index_gt_comb_adopt <- judge_gt_comb_mk(cand_gts, index_gt_comb, peaks, heights, sr_consider, bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr, phr_fltr, st)
#      n_gt_comb_all[n_mk * (i - 1) + j] <- length(index_gt_comb_adopt)
#      index_gt_comb_adopt <- unlist(index_gt_comb_adopt) + 1 # for R and stan index
#      index_gt_comb_all <- rbind(index_gt_comb_all, matrix(index_gt_comb_adopt, ncol = hnc, byrow = TRUE))
    }
  }

#  n_eph_base_all <- n_gt_all * n_peak_all

#  gamma_beta <- extract_beta(dt_parameter)

#  at_all <- make_at_all(n_csp, markers_auto, dt_at)

  ## check n_gt_comb_all == 0

#  data_list_main <- list(Hnc = hnc, N_Csp_Mk = n_csp_mk, N_Peak_All = n_peak_all, N_Eph_Base_All = n_eph_base_all,
#                         Eph_Base_All = eph_base_all, Size_Eph_All = size_eph_all,
#                         N_Gt_Comb_All = n_gt_comb_all, Index_Gt_Comb_All = index_gt_comb_all,
#                         Oph_All = unlist(heights_list), AT_All = at_all, Beta = gamma_beta)

#  model_main <- readRDS(paste0(path_pack, "/model/model_main.rds"))
#  fit_main <- model_main$sample(
#    data_list_main,
#    iter_warmup = 1000,
#    iter_sampling = 1000,
#    chains = 4,
#    parallel_chains = detectCores()
#  )

  save(list = c("eph_base_all", "m2sr_all_mk"), file = "test.Rdata")
}
