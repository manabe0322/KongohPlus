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

  n_csp_mk <- n_csp * n_mk
  eph_base_all <- list()
  index_dedup_eph_all <- NULL
  index_dedup_gt_all <- list()
  n_dedup_eph_all <- rep(0, n_csp_mk)
  size_eph_dedup_all <- list()
  oph_dedup_all <- list()
  n_peak_all <- rep(0, n_csp_mk)

  for(i in 1:n_csp){
    # Extract peak data of a crime stain profile
    peaks_list_sn1 <- peaks_list[[i]]
    sizes_list_sn1 <- sizes_list[[i]]
    heights_list_sn1 <- heights_list[[i]]

    for(j in 1:n_mk){
      mk <- markers_auto[j]

      # Extract peak data of one marker
      peaks <- peaks_list_sn1[[mk]]
      sizes <- sizes_list_sn1[[mk]]
      heights <- heights_list_sn1[[mk]]
      n_peak_all[n_mk * (i - 1) + j] <- n_peak <- length(peaks)

      # Extract parameters of one marker
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

      # The base of the expected peak heights
      eph_base <- make_eph_base(peaks, bsrs, fsrs, dsrs, m2sr, ae)
      eph_base_all[[n_mk * (i - 1) + j]] <- eph_base

      # All candidate genotypes of one contributor
      cand_gts <- combinations(n = n_peak, r = 2, v = peaks, repeats.allowed = TRUE)
      cand_gts <- asplit(cand_gts, 1)

      # Indices of genotypes of each contributor in all candidate genotype combinations
      index_gt_comb <- permutations(n = length(cand_gts), r = hnc, repeats.allowed = TRUE)
      index_gt_comb <- index_gt_comb - 1 # for C++ index
      index_gt_comb <- asplit(index_gt_comb, 1)

      # Indices of genotypes of each contributor in each adopted genotype combination
      index_gt_comb_adopt <- judge_gt_comb_mk(cand_gts, index_gt_comb, peaks, heights, sr_consider, bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr, phr_fltr, st)

      # Connect eph_base and genotype combinations
      index_eph_gtc <- connect_eph_gtc_mk(cand_gts, index_gt_comb_adopt, peaks)
      index_eph_gtc <- matrix(unlist(index_eph_gtc) + 1, ncol = 2 * hnc, byrow = TRUE) # for R and stan index
      index_eph_gtc[which(index_eph_gtc == 0, arr.ind = TRUE)] <- 5 * n_peak + 1
      index_dedup_eph <- index_eph_gtc[!duplicated(index_eph_gtc), , drop = FALSE]
      n_dedup_eph_all[n_mk * (i - 1) + j] <- nrow(index_dedup_eph)
      index_dedup_eph_all <- rbind(index_dedup_eph_all, index_dedup_eph)

      n_gtc_adopt <- length(index_gt_comb_adopt)
      size_eph <- rep(sizes, n_gtc_adopt)
      size_eph_dedup_all[[n_mk * (i - 1) + j]] <- size_eph[!duplicated(index_eph_gtc)]

      oph <- rep(heights, n_gtc_adopt)
      oph_dedup_all[[n_mk * (i - 1) + j]] <- oph[!duplicated(index_eph_gtc)]

      n_index_eph_gtc <- nrow(index_eph_gtc)
      n_al <- 2 * hnc
      index_dedup_gt <- rep(0, n_index_eph_gtc)
      bool_mat <- matrix(FALSE, nrow(index_dedup_eph), n_al)
      for(k in 1:n_index_eph_gtc){
        ind_eph_gtc1 <- index_eph_gtc[k, ]
        for(l in 1:n_al){
          bool_mat[, l] <- index_dedup_eph[, l] == ind_eph_gtc1[l]
        }
        index_dedup_gt[k] <- which(apply(bool_mat, 1, all))
      }
      index_dedup_gt_all[[n_mk * (i - 1) + j]] <- index_dedup_gt
    }
  }

  n_eph_base_all <- sapply(eph_base_all, length)

  n_index_dedup_gt_all <- sapply(index_dedup_gt_all, length)

  gamma_beta <- extract_beta(dt_parameter)

  at_all <- make_at_all(n_csp, markers_auto, dt_at)

  save(list = c("n_eph_base_all", "n_index_dedup_gt_all", "n_dedup_eph_all", "eph_base_all", "index_dedup_eph_all", "index_dedup_gt_all", "size_eph_dedup_all", "n_peak_all"), file = "test.Rdata")

  ## check n_gt_comb_all == 0

  eph_base_all <- unlist(eph_base_all)
  index_dedup_gt_all <- unlist(index_dedup_gt_all)
  size_eph_dedup_all <- unlist(size_eph_dedup_all)
  oph_dedup_all <- unlist(oph_dedup_all)

  data_list_main <- list(Hnc = hnc, N_Csp_Mk = n_csp_mk,
                         N_Eph_Base_All = n_eph_base_all, Eph_Base_All = eph_base_all,
                         N_Dedup_Eph_All = n_dedup_eph_all, Index_Dedup_Eph_All = index_dedup_eph_all, Size_Eph_Dedup_All = size_eph_dedup_all, Oph_Dedup_All = oph_dedup_all,
                         N_Index_Dedup_Gt_All = n_index_dedup_gt_all, Index_Dedup_Gt_All = index_dedup_gt_all,
                         N_Peak_All = n_peak_all,
                         AT_All = at_all, Beta = gamma_beta)

#  model_main2 <- readRDS(paste0(path_pack, "/model/model_main2.rds"))
#  fit_main <- model_main2$sample(
#    data_list_main,
#    iter_warmup = 1000,
#    iter_sampling = 1000,
#    chains = 4,
#    parallel_chains = detectCores()
#  )
}
