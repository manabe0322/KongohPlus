perform_deconvo <- function(dt_csp, dt_af, dt_kit_auto, dt_info_mk, dt_par, dt_sr, dt_al_prob, dt_gt_fltr, dt_dropin,
                            hnc, warmup, sampling, chains){
  # Start computation
  start_time <- proc.time()

  # Number of replicates
  n_rep <- length(unique(dt_csp[, SampleName]))

  # Extract AT
  markers_use <- dt_info_mk$Marker[!dt_info_mk$Ignore]
  ats_use <- dt_info_mk$Marker[!dt_info_mk$Ignore]

  # Extract peak data from crime stain profiles
  tmp <- make_peak_data(dt_csp, dt_kit_auto, markers_use, at_all_mk)
  peaks_list <- tmp$peaks_list
  sizes_list <- tmp$sizes_list
  heights_list <- tmp$heights_list
  min_size <- min(unlist(sizes_list))

  # Information on consideration of each stutter
  #sr_consider_all_mk <- make_sr_consider_all_mk(markers_auto, dt_par)

  # Extract parameters for AE
  ae_alpha_all_mk <- dt_par$Alpha[dt_par$Parameter == "AE"]
  #m2sr_all_mk <- extract_param_all_mk("M2SR", "Myu", markers_auto, dt_par)

  # Extract methods to estimate allele probabilities
  al_prob_subcategory <- dt_al_prob[, Subcategory]
  al_prob_value <- dt_al_prob[, Value]
  method_al_prob <- al_prob_value[al_prob_subcategory == "Method"]
  maf <- as.numeric(al_prob_value[al_prob_subcategory == "Maf"])
  theta <- as.numeric(al_prob_value[al_prob_param == "Theta"])

  # Set allele probabilities
  tmp <- set_al_prob(markers_auto, dt_af, peaks_list, method_al_prob, maf)
  pop_al_list <- tmp$pop_al_list
  pop_prob_list <- tmp$pop_prob_list

  # Extract filters to remove unrealistic genotype combinations
  gt_fltr_subcategory <- dt_gt_fltr[, Subcategory]
  gt_fltr_value <- as.numeric(dt_gt_fltr[, Value])
  bsr_fltr <- gt_fltr_value[gt_fltr_subcategory == "BSR"]
  fsr_fltr <- gt_fltr_value[gt_fltr_subcategory == "FSR"]
  dsr_fltr <- gt_fltr_value[gt_fltr_subcategory == "DSR"]
  m2sr_fltr <- gt_fltr_value[gt_fltr_subcategory == "M2SR"]
  n_dropin_fltr <- gt_fltr_value[gt_fltr_subcategory == "N_Dropin"]
  h_dropin_fltr <- gt_fltr_value[gt_fltr_subcategory == "H_Dropin"]

  # Extract drop-in parameters
  dropin_subcategory <- dt_dropin[, Subcategory]
  dropin_value <- as.numeric(dt_dropin[, Value])
  logprob_di <- log(dropin_value[dropin_subcategory == "Probability"])
  shape_di <- dropin_value[dropin_subcategory == "Shape"]
  rate_di <- dropin_value[dropin_subcategory == "Rate"]
  mt_di <- dropin_value[dropin_subcategory == "MT"]

  # Extract gamma parameter
  cv <- as.numeric(dt_par$Cv[dt_par$Parameter == "Gamma"])

  setkey(dt_alcor, Marker)

  # Define objects
  rph_base_all_mk <- vector(mode = "list", length = n_mk)
  gtc_logprob_all_mk <- vector(mode = "list", length = n_mk)
  index_rph_gtc_all_mk <- vector(mode = "list", length = n_mk)
  index_dedup1_general_all_mk <- vector(mode = "list", length = n_mk)
  index_dedup1_dropout_all_mk <- vector(mode = "list", length = n_mk)
  index_dedup1_dropin_all_mk <- vector(mode = "list", length = n_mk)
  index_dedup1_ignore_all_mk <- vector(mode = "list", length = n_mk)
  size_dedup1_general_all_mk <- vector(mode = "list", length = n_mk)
  size_dedup1_dropout_all_mk <- vector(mode = "list", length = n_mk)
  size_dedup1_dropin_all_mk <- vector(mode = "list", length = n_mk)
  size_dedup1_ignore_all_mk <- vector(mode = "list", length = n_mk)
  oph_dedup1_general_all_mk <- vector(mode = "list", length = n_mk)
  oph_dedup1_dropout_all_mk <- vector(mode = "list", length = n_mk)
  oph_dedup1_dropin_all_mk <- vector(mode = "list", length = n_mk)
  oph_dedup1_ignore_all_mk <- vector(mode = "list", length = n_mk)
  loglike_di_all_mk <- vector(mode = "list", length = n_mk)
  n_peak_all_mk <- rep(NA, n_mk)
  n_gtc_all_mk <- rep(NA, n_mk)
  gtc_all_mk <- vector(mode = "list", length = n_mk)
  max_tph <- 0

  names(gtc_all_mk) <- markers_auto

  for(i in 1:n_mk){
    # Extract peak data of one marker
    mk <- markers_auto[i]
    peaks <- peaks_list[[mk]]
    sizes <- sizes_list[[mk]]
    heights_all_rep <- heights_list[[mk]]
    heights <- heights_all_rep[[1]]
    oph_sum <- sum(heights)
    if(oph_sum > max_tph){
      max_tph <- oph_sum
    }
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
      gtc_logprob_all_mk[[i]] <- log(gtc_prob_mk)

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
      size_dedup1_dropin_all_mk[[i]] <- size_dedup1_tmp[index_dropin]
      size_dedup1_ignore_all_mk[[i]] <- size_dedup1_tmp[index_ignore]

      oph_dedup1_general_all_mk[[i]] <- oph_dedup1_tmp[index_general]
      oph_dedup1_dropout_all_mk[[i]] <- oph_dedup1_tmp[index_dropout]
      oph_dedup1_dropin_all_mk[[i]] <- oph_dedup1_tmp[index_dropin]
      oph_dedup1_ignore_all_mk[[i]] <- oph_dedup1_tmp[index_ignore]

      if(length(index_dropin) != 0){
        loglike_di_all_mk[[i]] <- log(dgamma(oph_dedup1_tmp[index_dropin] - mt_di + 1, shape_di, rate_di)) + logprob_di
      }
    }
  }

  ## check n_gt_comb_all == 0

  n_rph_base_all_mk <- sapply(rph_base_all_mk, length)
  n_general_all_mk <- sapply(index_dedup1_general_all_mk, ncol)
  n_dropout_all_mk <- sapply(index_dedup1_dropout_all_mk, ncol)
  n_dropin_all_mk <- sapply(index_dedup1_dropin_all_mk, ncol)
  n_ignore_all_mk <- sapply(index_dedup1_ignore_all_mk, ncol)
  n_general_sum <- sum(n_general_all_mk)
  n_dropout_sum <- sum(n_dropout_all_mk)
  n_dropin_sum <- sum(n_dropin_all_mk)
  n_ignore_sum <- sum(n_ignore_all_mk)
  n_dedup1_sum <- n_general_sum + n_dropout_sum + n_dropin_sum + n_ignore_sum
  id_mk_dedup1 <- unlist(mapply(rep, rep(1:n_mk, 4), c(n_general_all_mk, n_dropout_all_mk, n_dropin_all_mk, n_ignore_all_mk)))
  at_dropout <- unlist(mapply(rep, at_all_mk, n_dropout_all_mk))
  loglike_di <- unlist(loglike_di_all_mk)
  if(is.null(loglike_di)){
    loglike_di <- numeric(0)
  }

  # Arrange de-duplicated objects (pattern 1)
  index_dedup1_all_mk <- matrix(NA, nrow = 2 * hnc, ncol = n_dedup1_sum)
  size_dedup1_all_mk <- rep(NA, n_dedup1_sum)
  oph_dedup1_all_mk <- rep(NA, n_dedup1_sum)
  index_general_start <- 1
  index_dropout_start <- n_general_sum + 1
  index_dropin_start <- n_general_sum + n_dropout_sum + 1
  index_ignore_start <- n_general_sum + n_dropout_sum + n_dropin_sum + 1
  sumup_index_dedup1 <- 0
  for(i in 1:n_mk){
    n_general <- n_general_all_mk[i]
    if(n_general != 0){
      index_tmp <- index_general_start:(index_general_start + n_general - 1)
      index_dedup1_all_mk[, index_tmp] <- index_dedup1_general_all_mk[[i]] + sumup_index_dedup1
      size_dedup1_all_mk[index_tmp] <- size_dedup1_general_all_mk[[i]]
      oph_dedup1_all_mk[index_tmp] <- oph_dedup1_general_all_mk[[i]]
    }

    n_dropout <- n_dropout_all_mk[i]
    if(n_dropout != 0){
      index_tmp <- index_dropout_start:(index_dropout_start + n_dropout - 1)
      index_dedup1_all_mk[, index_tmp] <- index_dedup1_dropout_all_mk[[i]] + sumup_index_dedup1
      size_dedup1_all_mk[index_tmp] <- size_dedup1_dropout_all_mk[[i]]
      oph_dedup1_all_mk[index_tmp] <- oph_dedup1_dropout_all_mk[[i]]
    }

    n_dropin <- n_dropin_all_mk[i]
    if(n_dropin != 0){
      index_tmp <- index_dropin_start:(index_dropin_start + n_dropin - 1)
      index_dedup1_all_mk[, index_tmp] <- index_dedup1_dropin_all_mk[[i]] + sumup_index_dedup1
      size_dedup1_all_mk[index_tmp] <- size_dedup1_dropin_all_mk[[i]]
      oph_dedup1_all_mk[index_tmp] <- oph_dedup1_dropin_all_mk[[i]]
    }

    n_ignore <- n_ignore_all_mk[i]
    if(n_ignore != 0){
      index_tmp <- index_ignore_start:(index_ignore_start + n_ignore - 1)
      index_dedup1_all_mk[, index_tmp] <- index_dedup1_ignore_all_mk[[i]] + sumup_index_dedup1
      size_dedup1_all_mk[index_tmp] <- size_dedup1_ignore_all_mk[[i]]
      oph_dedup1_all_mk[index_tmp] <- oph_dedup1_ignore_all_mk[[i]]
    }

    index_general_start <- index_general_start + n_general
    index_dropout_start <- index_dropout_start + n_dropout
    index_dropin_start <- index_dropin_start + n_dropin
    index_ignore_start <- index_ignore_start + n_ignore
    sumup_index_dedup1 <- sumup_index_dedup1 + n_rph_base_all_mk[i]
  }

  # Define de-duplicated objects (pattern 2)
  n_dedup2_all_mk <- n_peak_all_mk * n_gtc_all_mk
  n_dedup2_sum <- sum(n_dedup2_all_mk)
  index_dedup2_all_mk <- rep(NA, n_dedup2_sum)
  index_dedup2_start <- 1
  sumup_index_dedup2 <- 0
  for(i in 1:n_mk){
    index_rph_gtc <- index_rph_gtc_all_mk[[i]] + sumup_index_dedup2
    n_dedup2 <- n_dedup2_all_mk[i]
    index_dedup2_all_mk[index_dedup2_start:(index_dedup2_start + n_dedup2 - 1)] <- make_index_dedup2_mk(index_rph_gtc, index_dedup1_all_mk)

    index_dedup2_start <- index_dedup2_start + n_dedup2
    sumup_index_dedup2 <- sumup_index_dedup2 + n_rph_base_all_mk[i]
  }

  data_list_main <- list(hnc = hnc, n_mk = n_mk, ltph = max_tph / 4, htph = max_tph / 2,
                         n_rph_base_sum = sum(n_rph_base_all_mk), rph_base_all_mk = unlist(rph_base_all_mk), ae_alpha_all_mk = ae_alpha_all_mk, cv = cv,
                         n_dedup1_sum = n_dedup1_sum, index_dedup1_all_mk = index_dedup1_all_mk, size_dedup1_all_mk = size_dedup1_all_mk - min_size, oph_dedup1_all_mk = oph_dedup1_all_mk, id_mk_dedup1 = id_mk_dedup1,
                         n_general_sum = n_general_sum, n_dropout_sum = n_dropout_sum, at_dropout = at_dropout, n_dropin_sum = n_dropin_sum, loglike_di = loglike_di,
                         n_dedup2_sum = n_dedup2_sum, n_dedup2_all_mk = n_dedup2_all_mk, index_dedup2_all_mk = index_dedup2_all_mk,
                         n_peak_all_mk = n_peak_all_mk, n_gtc_all_mk = n_gtc_all_mk, gtc_logprob_all_mk = unlist(gtc_logprob_all_mk))

  # Perform MCMC
  model_main <- stan_package_model(name = "model_main", package = "KongohPlus")
  fit_main <- model_main$sample(
    data_list_main,
    iter_warmup = warmup,
    iter_sampling = sampling,
    chains = chains,
    parallel_chains = detectCores()
  )

  # Trace plot
  draws <- fit_main$draws()
  draws %>% mcmc_trace(pars = c("lp__")) -> traceplot

  # Extract parameters
  fit_main_summary <- fit_main$summary()
  variables <- fit_main_summary$variable
  medians <- fit_main_summary$median
  rhats <- fit_main_summary$rhat
  lp <- medians[which(variables == "lp__")]
  tph <- medians[which(variables == "tph")]
  lambda <- medians[grep("lambda", variables)]
  mp <- medians[grep("mp", variables)]
  deg <- medians[grep("deg", variables)]
  pae <- medians[grep("pae", variables)]
  ae <- n_mk * pae
  rhat_lp <- rhats[which(variables == "lp__")]
  rhat_tph <- rhats[which(variables == "tph")]
  rhat_lambda <- rhats[grep("lambda", variables)]
  rhat_mp <- rhats[grep("mp", variables)]
  rhat_deg <- rhats[grep("deg", variables)]
  rhat_pae <- rhats[grep("pae", variables)]


  index_dedup1_all_mk_cpp <- asplit(index_dedup1_all_mk - 1, 1) # C++ index
  index_dedup2_all_mk_cpp <- index_dedup2_all_mk - 1 # C++ index
  ae_dedup1_all_mk <- ae[id_mk_dedup1]
  lpd_gtc_all_mk <- calc_lpd_gtc_all_mk_cpp(n_mk, mp, deg, hnc, n_dedup1_sum, n_general_sum, n_dropout_sum, n_dropin_sum,
                                            unlist(rph_base_all_mk), index_dedup1_all_mk_cpp, unlist(size_dedup1_all_mk), ae_dedup1_all_mk, oph_dedup1_all_mk,
                                            at_dropout, loglike_di, cv, tph,
                                            index_dedup2_all_mk_cpp, n_dedup2_all_mk, n_peak_all_mk, n_gtc_all_mk)

  dt_pg_list <- make_dt_pg_list(markers_auto, gtc_all_mk, lpd_gtc_all_mk, gtc_logprob_all_mk, paste0("U", 1:hnc))

  log10like_all_mk <- calc_log10like_all_mk(lpd_gtc_all_mk, gtc_logprob_all_mk)

  dt_log10like <- make_dt_log10like(markers_auto, log10like_all_mk)

  dt_mpdeg <- data.table(Contributor = paste0("U", 1:hnc),
                         MP = mp[hnc:1],
                         Deg = deg[hnc:1])

  dt_mcmc <- data.table(Parameter = c("Log probability in stan", "Total peak height", rep("lambda", hnc), rep("Mixture proportion", hnc), rep("Degradation", hnc), rep("Amplification efficiency", n_mk)),
                        Contributor = c("-", "-", paste0("U", 1:hnc), paste0("U", 1:hnc), paste0("U", 1:hnc), rep("-", n_mk)),
                        Marker = c("-", "-", rep("-", 3 * hnc), markers_auto),
                        Estimated_value = c(lp, tph, lambda[hnc:1], mp[hnc:1], deg[hnc:1], ae),
                        Rhat = c(rhat_lp, rhat_tph, rhat_lambda[hnc:1], rhat_mp[hnc:1], rhat_deg[hnc:1], rhat_pae))

  # End computation
  end_time <- proc.time()
  comp_time <- end_time[3] - start_time[3]

  return(list("hnc" = hnc,
              "markers_auto" = markers_auto,
              "pop_al_list" = pop_al_list,
              "pop_prob_list" = pop_prob_list,
              "theta" = theta,
              "gtc_all_mk" = gtc_all_mk,
              "lpd_gtc_all_mk" = lpd_gtc_all_mk,
              "dt_pg_list" = dt_pg_list,
              "dt_log10like" = dt_log10like,
              "dt_mpdeg" = dt_mpdeg,
              "dt_mcmc" = dt_mcmc,
              "traceplot" = traceplot,
              "comp_time" = comp_time))
}
