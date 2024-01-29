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

  n_cspmk <- n_csp * n_mk
  eph_base <- list()
  index_dedup_1 <- NULL
  index_dedup_2 <- list()
  n_dedup_1_cspmk <- rep(0, n_cspmk)
  size_dedup_1 <- list()
  oph_dedup_1 <- list()
  n_peak_cspmk <- rep(0, n_cspmk)

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
      n_peak_cspmk[n_mk * (i - 1) + j] <- n_peak <- length(peaks)

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
      eph_base[[n_mk * (i - 1) + j]] <- make_eph_base(peaks, bsrs, fsrs, dsrs, m2sr, ae)

      # All candidate genotypes of one contributor
      cand_gt <- combinations(n = n_peak, r = 2, v = peaks, repeats.allowed = TRUE)
      cand_gt <- asplit(cand_gt, 1)

      # Indices of genotypes of each contributor in all candidate genotype combinations
      index_gtc <- permutations(n = length(cand_gt), r = hnc, repeats.allowed = TRUE)
      index_gtc <- index_gtc - 1 # for C++ index
      index_gtc <- asplit(index_gtc, 1)

      # Indices of genotypes of each contributor in each adopted genotype combination
      index_gtc_adopt <- judge_gt_comb_mk(cand_gt, index_gtc, peaks, heights, sr_consider, bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr, phr_fltr, st)
      n_gtc_adopt <- length(index_gtc_adopt)

      # Connect eph_base and genotype combinations
      index_eph_gtc <- connect_eph_gtc_mk(cand_gt, index_gtc_adopt, peaks)
      index_eph_gtc <- matrix(unlist(index_eph_gtc) + 1, ncol = 2 * hnc, byrow = TRUE) # for R and stan index
      index_eph_gtc[which(index_eph_gtc == 0, arr.ind = TRUE)] <- 5 * n_peak + 1 # for density = 0

      # Judge duplicated eph combination
      bool_duplicated <- !duplicated(index_eph_gtc)

      # Delete duplicated indices -> de-duplicated objects (pattern 1)
      index_dedup_1_mk <- index_eph_gtc[bool_duplicated, , drop = FALSE]
      n_dedup_1_cspmk[n_mk * (i - 1) + j] <- nrow(index_dedup_1_mk)
      index_dedup_1 <- rbind(index_dedup_1, index_dedup_1_mk)

      # Sizes for de-duplicated objects (pattern 1)
      size_tmp <- rep(sizes, n_gtc_adopt)
      size_dedup_1[[n_mk * (i - 1) + j]] <- size_tmp[bool_duplicated]

      # Observed peak heights for de-duplicated objects (pattern 1)
      oph_tmp <- rep(heights, n_gtc_adopt)
      oph_dedup_1[[n_mk * (i - 1) + j]] <- oph_tmp[bool_duplicated]

      n_index_eph_gtc <- nrow(index_eph_gtc)
      n_al <- 2 * hnc
      index_dedup_2_mk <- rep(0, n_index_eph_gtc)
      bool_mat <- matrix(FALSE, nrow(index_dedup_1_mk), n_al)
      for(k in 1:n_index_eph_gtc){
        ind_eph_gtc1 <- index_eph_gtc[k, ]
        for(l in 1:n_al){
          bool_mat[, l] <- index_dedup_1_mk[, l] == ind_eph_gtc1[l]
        }
        index_dedup_2_mk[k] <- which(apply(bool_mat, 1, all))
      }
      index_dedup_2[[n_mk * (i - 1) + j]] <- index_dedup_2_mk
    }
  }

  n_eph_base <- sum(sapply(eph_base, length))

  n_dedup_2_cspmk <- sapply(index_dedup_2, length)
  n_gtc_cspmk <- n_dedup_2_cspmk / n_peak_cspmk

  gamma_beta <- extract_beta(dt_parameter)

  at_cspmk <- make_at_cspmk(n_csp, markers_auto, dt_at)

  save(list = c("n_eph_base", "n_dedup_2_cspmk", "n_dedup_1_cspmk", "eph_base", "index_dedup_1", "index_dedup_2", "size_dedup_1", "n_peak_cspmk", "n_gtc_cspmk"), file = "test.Rdata")

  ## check n_gt_comb_all == 0

  eph_base <- unlist(eph_base)
  index_dedup_2 <- unlist(index_dedup_2)
  size_dedup_1 <- unlist(size_dedup_1)
  oph_dedup_1 <- unlist(oph_dedup_1)

  data_list_main <- list(Hnc = hnc, N_Cspmk = n_cspmk,
                         N_Eph_Base = n_eph_base, Eph_Base = eph_base,
                         N_Dedup_1 = sum(n_dedup_1_cspmk), N_Dedup_1_Cspmk = n_dedup_1_cspmk, Index_Dedup_1 = index_dedup_1, Size_Dedup_1 = size_dedup_1, Oph_Dedup_1 = oph_dedup_1,
                         N_Dedup_2 = sum(n_dedup_2_cspmk), Index_Dedup_2 = index_dedup_2,
                         N_Peak_Cspmk = n_peak_cspmk, N_Gtc_Cspmk = n_gtc_cspmk,
                         At_Cspmk = at_cspmk, Beta = gamma_beta)

  model_main <- readRDS(paste0(path_pack, "/model/model_main.rds"))
  fit_main <- model_main$sample(
    data_list_main,
    iter_warmup = 1000,
    iter_sampling = 1000,
    chains = 4,
    parallel_chains = detectCores()
  )

  fit_main$summary()
}
