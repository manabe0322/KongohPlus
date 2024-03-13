perform_deconvo <- function(path_pack, dt_csp, dt_af, dt_kit, dt_at, dt_param, dt_alcor, dt_gtc_fltr, dt_gtc_prob, hnc,
                            warmup, sampling, chains, di_alpha = 1.22, di_beta = 0.0503, logprob_di = log(0.00335), di_min = 49){
  # Start computation
  start_time <- proc.time()

  # Extract autosomal markers
  dt_kit_auto <- extract_kit_auto(dt_kit)
  markers_auto <- unique(dt_kit_auto[, Marker])
  n_mk <- length(markers_auto)
  n_rep <- length(unique(dt_csp[, SampleName]))

  # Extract AT
  at_all_mk <- extract_at_all_mk(markers_auto, dt_at)

  # Extract peak data from crime stain profiles
  tmp <- make_peak_data(dt_csp, dt_kit_auto, markers_auto, at_all_mk)
  peaks_list <- tmp$peaks_list
  sizes_list <- tmp$sizes_list
  heights_list <- tmp$heights_list

  # Information on consideration of each stutter
  sr_consider_all_mk <- make_sr_consider_all_mk(markers_auto, dt_param)

  # Extract parameters for AE and M2SR
  ae_alpha_all_mk <- extract_param_all_mk("AE", "Alpha", markers_auto, dt_param)
  m2sr_all_mk <- extract_param_all_mk("M2SR", "Myu", markers_auto, dt_param)

  # Extract filters to remove unrealistic genotype combinations
  gtc_fltr_type <- dt_gtc_fltr[, Type]
  gtc_fltr_value <- dt_gtc_fltr[, Value]
  bsr_fltr <- gtc_fltr_value[gtc_fltr_type == "BSR"]
  fsr_fltr <- gtc_fltr_value[gtc_fltr_type == "FSR"]
  dsr_fltr <- gtc_fltr_value[gtc_fltr_type == "DSR"]
  m2sr_fltr <- gtc_fltr_value[gtc_fltr_type == "M2SR"]
  n_dropin_fltr <- gtc_fltr_value[gtc_fltr_type == "N_DropIn"]
  h_dropin_fltr <- gtc_fltr_value[gtc_fltr_type == "H_DropIn"]

  # Extract methods to estimate genotype probabilities
  gtc_prob_param <- dt_gtc_prob[, Parameter]
  gtc_prob_value <- dt_gtc_prob[, Value]
  theta <- as.numeric(gtc_prob_value[gtc_prob_param == "Theta"])

  # Set allele probabilities
  tmp <- set_al_prob(markers_auto, dt_af, peaks_list)
  pop_al_list <- tmp$pop_al_list
  pop_prob_list <- tmp$pop_prob_list

  # Extract gamma parameter
  cv <- extract_cv(dt_param)

  setkey(dt_alcor, Marker)

  # Define objects
  rph_base_all_mk <- vector(mode = "list", length = n_mk)
  index_rph_gtc_all_mk <- vector(mode = "list", length = n_mk)
  index_dedup1_general_all_mk <- vector(mode = "list", length = n_mk)
  index_dedup1_dropout_all_mk <- vector(mode = "list", length = n_mk)
  index_dedup1_dropin_all_mk <- vector(mode = "list", length = n_mk)
  index_dedup1_ignore_all_mk <- vector(mode = "list", length = n_mk)
  size_dedup1_general_all_mk <- vector(mode = "list", length = n_mk)
  size_dedup1_dropout_all_mk <- vector(mode = "list", length = n_mk)
  oph_dedup1_general_all_mk <- vector(mode = "list", length = n_mk)
  lpd_di_all_mk <- vector(mode = "list", length = n_mk)
  n_peak_all_mk <- rep(NA, n_mk)
  n_gtc_all_mk <- rep(NA, n_mk)
  gtc_all_mk <- vector(mode = "list", length = n_mk)

  for(i in 1:n_mk){
    # Extract peak data of one marker
    mk <- markers_auto[i]
    peaks <- peaks_list[[mk]]
    sizes <- sizes_list[[mk]]
    heights_all_rep <- heights_list[[mk]]
    heights <- heights_all_rep[[1]]
    n_peak_all_mk[i] <- n_peak <- length(peaks)

    # Extract parameters of one marker
    at <- at_all_mk[mk]
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

    # Extract allele probabilities
    pop_al <- pop_al_list[[mk]]
    pop_prob <- pop_prob_list[[mk]]

    # The base of the relative peak heights
    rph_base_all_mk[[i]] <- rph_base_mk <- make_rph_base(peaks, bsrs, fsrs, dsrs, m2sr)

    # All candidate genotypes of one contributor
    cand_gt <- combinations(n = n_peak, r = 2, v = peaks, repeats.allowed = TRUE)
    cand_gt <- asplit(cand_gt, 1)

    # Indices of genotypes of each contributor in all candidate genotype combinations
    index_gtc <- permutations(n = length(cand_gt), r = hnc, repeats.allowed = TRUE)
    index_gtc_cpp <- index_gtc - 1 # C++ index
    index_gtc_cpp <- asplit(index_gtc_cpp, 1)

    # Indices of genotypes of each contributor in each adopted genotype combination
    index_gtc_adopt_cpp <- judge_gtc_mk(cand_gt, index_gtc_cpp, peaks, heights_all_rep, sr_consider, bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr, n_dropin_fltr, h_dropin_fltr)
    n_gtc_all_mk[i] <- n_gtc <- length(index_gtc_adopt_cpp)

    if(n_gtc != 0){
      # Extract genotype combinations
      gtc_all_mk[[i]] <- gtc_mk <- extract_gtc_adopt(cand_gt, index_gtc_adopt_cpp)

      # Calculate genotype probabilities
      n_obs <- make_n_obs(pop_al, prof_ref = numeric(0))
      gtc_prob_mk <- calc_gtc_prob_mk(gtc_mk, pop_al, pop_prob, theta, prof_ref = numeric(0), hyp_id = rep(0, hnc), n_obs)
      gtc_logprob_mk <- log(gtc_prob_mk)

      # Connect rph_base and genotype combinations
      index_rph_gtc_cpp <- connect_rph_gtc_mk(cand_gt, index_gtc_adopt_cpp, peaks, sr_consider)
      index_rph_gtc_all_mk[[i]] <- index_rph_gtc <- matrix(unlist(index_rph_gtc_cpp) + 1, nrow = 2 * hnc) # R and stan index

      # Delete duplicated indices -> index for deduplication (pattern 1: Before summation of the rph_base)
      bool_duplicated_adopt <- !duplicated(t(index_rph_gtc))
      index_dedup1_tmp <- index_rph_gtc[, bool_duplicated_adopt, drop = FALSE]

      # Sizes for de-duplicated objects (pattern 1)
      size_dedup1_tmp <- rep(sizes, n_gtc)
      size_dedup1_tmp <- size_dedup1_tmp[bool_duplicated_adopt]

      # Observe peak heights for de-duplicated objects (pattern 1)
      oph_dedup1_tmp <- rep(heights, n_gtc)
      oph_dedup1_tmp <- oph_dedup1_tmp[bool_duplicated_adopt]

      # Indices of each peak type
      rph_dedup1_tmp <- matrix(rph_base_mk[index_dedup1_tmp], nrow = 2 * hnc)
      rph_dedup1_tmp <- apply(rph_dedup1_tmp, 2, sum)
      index_ignore <- which(rph_dedup1_tmp == 0)
      index_dropin <- which(rph_dedup1_tmp < 0)
      index_oph0 <- which(oph_dedup1_tmp == 0)
      index_dropout <- setdiff(index_oph0, index_ignore)
      index_general <- setdiff(1:length(rph_dedup1_tmp), c(index_ignore, index_dropin, index_dropout))

      index_dedup1_general_all_mk[[i]] <- index_dedup1_tmp[, index_general, drop = FALSE]
      index_dedup1_dropout_all_mk[[i]] <- index_dedup1_tmp[, index_dropout, drop = FALSE]
      index_dedup1_dropin_all_mk[[i]] <- index_dedup1_tmp[, index_dropin, drop = FALSE]
      index_dedup1_ignore_all_mk[[i]] <- index_dedup1_tmp[, index_ignore, drop = FALSE]

      size_dedup1_general_all_mk[[i]] <- size_dedup1_tmp[index_general]
      size_dedup1_dropout_all_mk[[i]] <- size_dedup1_tmp[index_dropout]
      oph_dedup1_general_all_mk[[i]] <- oph_dedup1_tmp[index_general]

      if(length(index_dropin) != 0){
        lpd_di_all_mk[[i]] <- log(dgamma(oph_dedup1_tmp[index_dropin] - di_min, di_alpha, di_beta)) + logprob_di
      }
    }
  }

  n_rph_base_all_mk <- sapply(rph_base_all_mk, length)
  n_general_all_mk <- sapply(index_dedup1_general_all_mk, ncol)
  n_dropout_all_mk <- sapply(index_dedup1_dropout_all_mk, ncol)
  n_dropin_all_mk <- sapply(index_dedup1_dropin_all_mk, ncol)
  n_ignore_all_mk <- sapply(index_dedup1_ignore_all_mk, ncol)

  n_general_sum <- sum(n_general_all_mk)
  n_dropout_sum <- sum(n_dropout_all_mk)
  n_dropin_sum <- sum(n_dropin_all_mk)
  n_ignore_sum <- sum(n_ignore_all_mk)

  index_dedup1_all_mk <- matrix(NA, nrow = 2 * hnc, ncol = sum(c(n_general_sum, n_dropout_sum + n_dropin_sum + n_ignore_sum)))
  index_general_start <- 1
  index_dropout_start <- n_general_sum + 1
  index_dropin_start <- n_general_sum + n_dropout_sum + 1
  index_ignore_start <- n_general_sum + n_dropout_sum + n_dropin_sum + 1

  sumup_index_dedup1 <- 0
  for(i in 1:n_mk){
    n_general <- n_general_all_mk[i]
    if(n_general != 0){
      index_dedup1_all_mk[, index_general_start:(index_general_start + n_general - 1)] <- index_dedup1_general_all_mk[[i]] + sumup_index_dedup1
    }

    n_dropout <- n_dropout_all_mk[i]
    if(n_dropout != 0){
      index_dedup1_all_mk[, index_dropout_start:(index_dropout_start + n_dropout - 1)] <- index_dedup1_dropout_all_mk[[i]] + sumup_index_dedup1
    }

    n_dropin <- n_dropin_all_mk[i]
    if(n_dropin != 0){
      index_dedup1_all_mk[, index_dropin_start:(index_dropin_start + n_dropin - 1)] <- index_dedup1_dropin_all_mk[[i]] + sumup_index_dedup1
    }

    n_ignore <- n_ignore_all_mk[i]
    if(n_ignore != 0){
      index_dedup1_all_mk[, index_ignore_start:(index_ignore_start + n_ignore - 1)] <- index_dedup1_ignore_all_mk[[i]] + sumup_index_dedup1
    }

    index_general_start <- index_general_start + n_general
    index_dropout_start <- index_dropout_start + n_dropout
    index_dropin_start <- index_dropin_start + n_dropin
    index_ignore_start <- index_ignore_start + n_ignore
    sumup_index_dedup1 <- sumup_index_dedup1 + n_rph_base_all_mk[i]
  }

  n_dedup2_all_mk <- n_peak_all_mk * n_gtc_all_mk
  index_dedup2_all_mk <- rep(NA, sum(n_dedup2_all_mk))
  index_dedup2_start <- 1
  sumup_index_dedup2 <- 0
  for(i in 1:n_mk){
    index_rph_gtc <- index_rph_gtc_all_mk[[i]] + sumup_index_dedup2
    n_dedup2 <- n_dedup2_all_mk[i]
    index_dedup2_all_mk[index_dedup2_start:(index_dedup2_start + n_dedup2 - 1)] <- make_index_dedup2_mk(index_rph_gtc, index_dedup1_all_mk)

    index_dedup2_start <- index_dedup2_start + n_dedup2
    sumup_index_dedup2 <- sumup_index_dedup2 + n_rph_base_all_mk[i]
  }
}
