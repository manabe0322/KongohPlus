perform_deconvo <- function(path_pack, dt_csp, dt_af, dt_kit, dt_at, dt_param, dt_alcor, dt_gtc_fltr, dt_gtc_prob, hnc){
  # Extract autosomal markers
  dt_kit_auto <- extract_kit_auto(dt_kit)
  markers_auto <- unique(dt_kit_auto[, Marker])
  n_mk <- length(markers_auto)

  # Extract AT
  n_csp <- length(unique(dt_csp[, SampleName]))
  at_cspmk <- make_at_cspmk(n_csp, markers_auto, dt_at)

  # Extract peak data from crime stain profiles
  tmp <- make_peak_data(dt_csp, dt_kit, markers_auto, at_cspmk)
  peaks_list <- tmp$peaks_list
  sizes_list <- tmp$sizes_list
  heights_list <- tmp$heights_list

  # Information on consideration of each stutter
  sr_consider_all_mk <- make_sr_consider_all_mk(markers_auto, dt_param)

  # Extract the 'myu' parameters for AE and M2SR
  ae_all_mk <- extract_myu_all_mk("AE", markers_auto, dt_param)
  m2sr_all_mk <- extract_myu_all_mk("M2SR", markers_auto, dt_param)

  # Extract filters to remove unrealistic genotype combinations
  gtc_fltr_type <- dt_gtc_fltr[, Type]
  gtc_fltr_value <- dt_gtc_fltr[, Value]
  bsr_fltr <- gtc_fltr_value[gtc_fltr_type == "BSR"]
  fsr_fltr <- gtc_fltr_value[gtc_fltr_type == "FSR"]
  dsr_fltr <- gtc_fltr_value[gtc_fltr_type == "DSR"]
  m2sr_fltr <- gtc_fltr_value[gtc_fltr_type == "M2SR"]
  phr_fltr <- gtc_fltr_value[gtc_fltr_type == "PHR"]
  st <- gtc_fltr_value[gtc_fltr_type == "ST"]

  # Extract methods to estimate genotype probabilities
  gtc_prob_param <- dt_gtc_prob[, Parameter]
  gtc_prob_value <- dt_gtc_prob[, Value]
  theta <- as.numeric(gtc_prob_value[gtc_prob_param == "Theta"])

  # Set allele probabilities
  tmp <- set_al_prob(markers_auto, dt_af, peaks_list)
  pop_al_list <- tmp$pop_al_list
  pop_prob_list <- tmp$pop_prob_list

  # Extract gamma parameter
  gamma_beta <- extract_beta(dt_param)

  setkey(dt_alcor, Marker)

  n_cspmk <- n_csp * n_mk
  n_eph_base_all <- rep(0, n_cspmk)
  eph_base <- vector(mode = "list", length = n_cspmk)
  n_dedup_1_cspmk <- rep(0, n_cspmk)
  index_dedup_1 <- NULL
  size_dedup_1 <- vector(mode = "list", length = n_cspmk)
  oph_dedup_1 <- vector(mode = "list", length = n_cspmk)
  n_dedup_2_cspmk <- rep(0, n_cspmk)
  index_dedup_2 <- vector(mode = "list", length = n_cspmk)
  n_peak_cspmk <- rep(0, n_cspmk)
  n_gtc_cspmk <- rep(0, n_cspmk)
  gtc_logprob_all <- vector(mode = "list", length = n_cspmk)
  gtc_all <- vector(mode = "list", length = n_cspmk)

  sumup_index_dedup_1 <- 0
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

      # Extract allele probabilities
      pop_al <- pop_al_list[[mk]]
      pop_prob <- pop_prob_list[[mk]]

      # The base of the expected peak heights
      eph_base_mk <- make_eph_base(peaks, bsrs, fsrs, dsrs, m2sr, ae)
      eph_base[[n_mk * (i - 1) + j]] <- eph_base_mk
      n_eph_base_all[n_mk * (i - 1) + j] <- length(eph_base_mk)

      # All candidate genotypes of one contributor
      cand_gt <- combinations(n = n_peak, r = 2, v = peaks, repeats.allowed = TRUE)
      cand_gt <- asplit(cand_gt, 1)

      # Indices of genotypes of each contributor in all candidate genotype combinations
      index_gtc <- permutations(n = length(cand_gt), r = hnc, repeats.allowed = TRUE)
      index_gtc_cpp <- index_gtc - 1 # C++ index
      index_gtc_cpp <- asplit(index_gtc_cpp, 1)

      # Indices of genotypes of each contributor in each adopted genotype combination
      index_gtc_adopt_cpp <- judge_gt_comb_mk(cand_gt, index_gtc_cpp, peaks, heights, sr_consider, bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr, phr_fltr, st)
      n_gtc_cspmk[n_mk * (i - 1) + j] <- n_gtc <- length(index_gtc_adopt_cpp)

      if(n_gtc != 0){
        # Extract genotype combinations
        gtc_mk <- extract_gtc_adopt(cand_gt, index_gtc_adopt_cpp)
        gtc_all[[n_mk * (i - 1) + j]] <- matrix(unlist(gtc_mk), ncol = 2 * hnc, byrow = TRUE)

        # Calculate genotype probabilities
        n_obs <- make_n_obs(pop_al, prof_ref = numeric(0))
        gtc_prob_mk <- calc_gtc_prob_mk(gtc_mk, pop_al, pop_prob, theta, prof_ref = numeric(0), hyp_id = c(0, 0), n_obs)
        gtc_logprob_all[[n_mk * (i - 1) + j]] <- log(gtc_prob_mk)

        # Connect eph_base and genotype combinations
        index_eph_gtc_cpp <- connect_eph_gtc_mk(cand_gt, index_gtc_adopt_cpp, peaks)
        index_eph_gtc <- matrix(unlist(index_eph_gtc_cpp) + 1, ncol = 2 * hnc, byrow = TRUE) # R and stan index
        index_eph_gtc[which(index_eph_gtc == 0, arr.ind = TRUE)] <- 5 * n_peak + 1 # for density = 0

        # Judge duplicated eph combination
        bool_duplicated <- !duplicated(index_eph_gtc)

        # Delete duplicated indices -> index for deduplication (pattern 1: Before summation of the eph_base)
        index_dedup_1_mk <- index_eph_gtc[bool_duplicated, , drop = FALSE]
        n_dedup_1_cspmk[n_mk * (i - 1) + j] <- nrow(index_dedup_1_mk)

        # Sizes for de-duplicated objects (pattern 1)
        size_tmp <- rep(sizes, n_gtc)
        size_dedup_1[[n_mk * (i - 1) + j]] <- size_tmp[bool_duplicated]

        # Observed peak heights for de-duplicated objects (pattern 1)
        oph_tmp <- rep(heights, n_gtc)
        oph_dedup_1[[n_mk * (i - 1) + j]] <- oph_tmp[bool_duplicated]

        # Index for deduplication (pattern 2: After summation of the eph_base)
        index_dedup_2_mk <- make_index_dedup_2_mk(index_eph_gtc, 2 * hnc, index_dedup_1_mk)
        index_dedup_2[[n_mk * (i - 1) + j]] <- index_dedup_2_mk
        n_dedup_2_cspmk[n_mk * (i - 1) + j] <- length(index_dedup_2_mk)

        # Sequential numbering of index_dedup_1_mk
        index_dedup_1_mk <- index_dedup_1_mk + sumup_index_dedup_1 # index for combining all loci
        index_dedup_1 <- rbind(index_dedup_1, index_dedup_1_mk)
        sumup_index_dedup_1 <- sumup_index_dedup_1 + length(eph_base_mk)
      }
    }
  }

  ## check n_gt_comb_all == 0

  data_list_main <- list(Hnc = hnc, N_Cspmk = n_cspmk,
                         N_Eph_Base_All = sum(n_eph_base_all), Eph_Base = unlist(eph_base),
                         N_Dedup_1_All = sum(n_dedup_1_cspmk), N_Dedup_1_Cspmk = n_dedup_1_cspmk, Index_Dedup_1 = index_dedup_1, Size_Dedup_1 = unlist(size_dedup_1), Oph_Dedup_1 = unlist(oph_dedup_1),
                         N_Dedup_2_All = sum(n_dedup_2_cspmk), N_Dedup_2_Cspmk = n_dedup_2_cspmk, Index_Dedup_2 = unlist(index_dedup_2),
                         N_Peak_Cspmk = n_peak_cspmk, N_Gtc_Cspmk = n_gtc_cspmk, Gtc_Logprob_All = unlist(gtc_logprob_all),
                         At_Cspmk = at_cspmk, Beta = gamma_beta)

  # Perform MCMC
  model_main <- readRDS(paste0(path_pack, "/model/model_main.rds"))
  fit_main <- model_main$sample(
    data_list_main,
    iter_warmup = 1000,
    iter_sampling = 1000,
    chains = 4,
    parallel_chains = detectCores()
  )

  # Trace plot
  draws <- fit_main$draws()
  draws %>% mcmc_trace(pars = c("lp__")) -> traceplot

  # Extract parameters
  fit_main_summary <- fit_main$summary()
  tph <- fit_main_summary$median[which(fit_main_summary$variable == "tph")]
  rhat_lp <- fit_main_summary$rhat[which(fit_main_summary$variable == "lp__")]

  # Determine the mixture proportions and the degradation parameters
  sampling_data <- as_draws_matrix(draws)
  mpdeg_final <- determine_mpdeg(sampling_data, n_sampling = 1000, n_chain = 4, upper_rhat = 1.1)
  mp <- mpdeg_final$mp_final
  deg <- mpdeg_final$deg_final

  # Calculate likelihood
  index_dedup_1_cpp <- asplit(index_dedup_1 - 1, 2) # C++ index
  eph_dedup_1 <- calc_eph_dedup_1_cpp(tph, mp, deg, hnc, sum(n_dedup_1_cspmk), unlist(eph_base), index_dedup_1_cpp, unlist(size_dedup_1))
  index_dedup_2_cpp <- list()
  for(i in 1:length(index_dedup_2)){
    index_dedup_2_cpp[[i]] <- index_dedup_2[[i]] - 1
  }

  lpd_gtc_all <- calc_lpd_gtc_all(n_cspmk, n_dedup_1_cspmk, unlist(oph_dedup_1), eph_dedup_1, at_cspmk, gamma_beta, n_peak_cspmk, n_gtc_cspmk, index_dedup_2_cpp)

  dt_pg_list <- make_dt_pg_list(markers_auto, gtc_all, lpd_gtc_all, gtc_logprob_all)

  likelihood_all <- calc_likelihood_all(lpd_gtc_all, gtc_logprob_all)

  dt_likelihood <- make_dt_likelihood(markers_auto, likelihood_all)

  dt_mpdeg <- data.table(Contributor = paste0("U", 1:hnc),
                         MP = mp,
                         Deg = deg)

  return(list("dt_pg_list" = dt_pg_list,
              "dt_likelihood" = dt_likelihood,
              "dt_mpdeg" = dt_mpdeg,
              "tph" = tph,
              "rhat_lp" = rhat_lp,
              "traceplot" = traceplot))
}
