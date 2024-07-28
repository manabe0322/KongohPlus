perform_deconvo <- function(dt_csp, dt_af, dt_kit_auto, dt_info_mk, dt_par, dt_sr, dt_al_prob, dt_gt_fltr, dt_dropin,
                            hnc, warmup, sampling, chains){
  # Start computation
  start_time <- proc.time()

  # Number of replicates
  n_rep <- length(unique(dt_csp[, SampleName]))

  # Extract used markers and ATs
  markers_use <- dt_info_mk$Marker[!dt_info_mk$Ignore]
  ats_use <- dt_info_mk$AT[!dt_info_mk$Ignore]
  n_mk <- length(markers_use)

  if(n_mk != 0){
    # Extract peak data from crime stain profiles
    tmp <- make_peak_data(dt_csp, dt_kit_auto, markers_use, ats_use)
    peaks_list <- tmp$peaks_list
    sizes_list <- tmp$sizes_list
    heights_list <- tmp$heights_list
    min_size <- min(unlist(sizes_list))

    # Extract parameters for AE
    ae_alpha_all_mk <- dt_par$Alpha[dt_par$Parameter == "AE"]
    ae_alpha_use <- ae_alpha_all_mk[is.element(dt_par$Marker[dt_par$Parameter == "AE"], markers_use)]

    # Extract gamma parameter
    cv <- dt_par$Cv[dt_par$Parameter == "Gamma"]

    # Extract methods to estimate allele probabilities
    al_prob_subcategory <- dt_al_prob[, Subcategory]
    al_prob_value <- dt_al_prob[, Value]
    method_al_prob <- al_prob_value[al_prob_subcategory == "Method"]
    maf <- as.numeric(al_prob_value[al_prob_subcategory == "Maf"])
    theta <- as.numeric(al_prob_value[al_prob_subcategory == "Theta"])

    # Set allele probabilities
    tmp <- set_al_prob(markers_use, dt_af, peaks_list, method_al_prob, maf)
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

    # Setkey of data.table
    setkey(dt_sr, Marker)

    # Define objects
    rph_base_use <- vector(mode = "list", length = n_mk)
    gtc_logprob_use <- vector(mode = "list", length = n_mk)
    index_rph_gtc_use <- vector(mode = "list", length = n_mk)
    index_dedup1_general_use <- vector(mode = "list", length = n_mk)
    index_dedup1_dropout_use <- vector(mode = "list", length = n_mk)
    index_dedup1_dropin_use <- vector(mode = "list", length = n_mk)
    index_dedup1_ignore_use <- vector(mode = "list", length = n_mk)
    size_dedup1_general_use <- vector(mode = "list", length = n_mk)
    size_dedup1_dropout_use <- vector(mode = "list", length = n_mk)
    size_dedup1_dropin_use <- vector(mode = "list", length = n_mk)
    size_dedup1_ignore_use <- vector(mode = "list", length = n_mk)
    oph_dedup1_general_use <- vector(mode = "list", length = n_mk)
    oph_dedup1_dropout_use <- vector(mode = "list", length = n_mk)
    oph_dedup1_dropin_use <- vector(mode = "list", length = n_mk)
    oph_dedup1_ignore_use <- vector(mode = "list", length = n_mk)
    loglike_di_use <- vector(mode = "list", length = n_mk)
    n_peak_use <- rep(NA, n_mk)
    n_gtc_use <- rep(NA, n_mk)
    gtc_use <- vector(mode = "list", length = n_mk)
    max_tph <- 0

    names(gtc_use) <- markers_use

    for(i in 1:n_mk){
      # Extract peak data of one marker
      mk <- markers_use[i]
      peaks <- peaks_list[[mk]]
      sizes <- sizes_list[[mk]]
      heights_all_rep <- heights_list[[mk]]
      heights <- heights_all_rep[[1]]
      n_peak_use[i] <- n_peak <- length(peaks)

      # Update the maximum of the total peak height per marker
      oph_sum <- sum(heights)
      if(oph_sum > max_tph){
        max_tph <- oph_sum
      }

      # Extract parameters of one marker
      at <- ats_use[i]
      dt_sr_mk <- dt_sr[mk]
      sr_alleles_mk <- as.numeric(dt_sr_mk[, Allele])
      bsrs_mk <- dt_sr_mk[, Average_BSR]
      fsrs_mk <- dt_sr_mk[, Average_FSR]
      dsrs_mk <- dt_sr_mk[, Average_DSR]
      m2srs_mk <- dt_sr_mk[, Average_M2SR]
      index_sr <- match(peaks, sr_alleles_mk)
      bsrs <- fsrs <- dsrs <- m2srs <- rep(0, n_peak)
      for(j in 1:length(index_sr)){
        if(!is.na(index_sr[j])){
          bsrs[j] <- bsrs_mk[index_sr[j]]
          fsrs[j] <- fsrs_mk[index_sr[j]]
          dsrs[j] <- dsrs_mk[index_sr[j]]
          m2srs[j] <- m2srs_mk[index_sr[j]]
        }
      }

      # Extract allele probabilities
      pop_al <- pop_al_list[[mk]]
      pop_prob <- pop_prob_list[[mk]]

      # The base of the relative peak heights
      rph_base_use[[i]] <- rph_base_mk <- make_rph_base(peaks, bsrs, fsrs, dsrs, m2srs)

      # All candidate genotypes of one contributor
      cand_gt <- combinations(n = n_peak, r = 2, v = peaks, repeats.allowed = TRUE)
      cand_gt <- asplit(cand_gt, 1)

      # Indices of genotypes of each contributor in all candidate genotype combinations
      index_gtc <- permutations(n = length(cand_gt), r = hnc, repeats.allowed = TRUE)
      index_gtc_cpp <- index_gtc - 1 # C++ index
      index_gtc_cpp <- asplit(index_gtc_cpp, 1)

      # Indices of genotypes of each contributor in each adopted genotype combination
      index_gtc_adopt_cpp <- judge_gtc_mk(cand_gt, index_gtc_cpp, peaks, heights_all_rep,
                                          bsrs, fsrs, dsrs, m2srs,
                                          bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr, n_dropin_fltr, h_dropin_fltr)
      n_gtc_use[i] <- n_gtc <- length(index_gtc_adopt_cpp)

      if(n_gtc != 0){
        # Extract genotype combinations
        gtc_use[[i]] <- gtc_mk <- extract_gtc_adopt(cand_gt, index_gtc_adopt_cpp)

        # Calculate genotype probabilities
        n_obs <- make_n_obs(pop_al, prof_ref = numeric(0))
        gtc_prob_mk <- calc_gtc_prob_mk(gtc_mk, pop_al, pop_prob, theta, prof_ref = numeric(0), hyp_id = rep(0, hnc), n_obs)
        gtc_logprob_use[[i]] <- log(gtc_prob_mk)

        # Connect rph_base and genotype combinations
        index_rph_gtc_cpp <- connect_rph_gtc_mk(cand_gt, index_gtc_adopt_cpp, peaks, bsrs, fsrs, dsrs, m2srs)
        index_rph_gtc_use[[i]] <- index_rph_gtc <- matrix(unlist(index_rph_gtc_cpp) + 1, nrow = 2 * hnc) # R and stan index

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

        index_dedup1_general_use[[i]] <- index_dedup1_tmp[, index_general, drop = FALSE]
        index_dedup1_dropout_use[[i]] <- index_dedup1_tmp[, index_dropout, drop = FALSE]
        index_dedup1_dropin_use[[i]] <- index_dedup1_tmp[, index_dropin, drop = FALSE]
        index_dedup1_ignore_use[[i]] <- index_dedup1_tmp[, index_ignore, drop = FALSE]

        size_dedup1_general_use[[i]] <- size_dedup1_tmp[index_general]
        size_dedup1_dropout_use[[i]] <- size_dedup1_tmp[index_dropout]
        size_dedup1_dropin_use[[i]] <- size_dedup1_tmp[index_dropin]
        size_dedup1_ignore_use[[i]] <- size_dedup1_tmp[index_ignore]

        oph_dedup1_general_use[[i]] <- oph_dedup1_tmp[index_general]
        oph_dedup1_dropout_use[[i]] <- oph_dedup1_tmp[index_dropout]
        oph_dedup1_dropin_use[[i]] <- oph_dedup1_tmp[index_dropin]
        oph_dedup1_ignore_use[[i]] <- oph_dedup1_tmp[index_ignore]

        if(length(index_dropin) != 0){
          loglike_di_use[[i]] <- log(dgamma(oph_dedup1_tmp[index_dropin] - mt_di + 1, shape_di, rate_di)) + logprob_di
        }
      }
    }

    ## check n_gt_comb_all == 0

    n_rph_base_use <- sapply(rph_base_use, length)
    n_general_use <- sapply(index_dedup1_general_use, ncol)
    n_dropout_use <- sapply(index_dedup1_dropout_use, ncol)
    n_dropin_use <- sapply(index_dedup1_dropin_use, ncol)
    n_ignore_use <- sapply(index_dedup1_ignore_use, ncol)
    n_general_sum <- sum(n_general_use)
    n_dropout_sum <- sum(n_dropout_use)
    n_dropin_sum <- sum(n_dropin_use)
    n_ignore_sum <- sum(n_ignore_use)
    n_dedup1_sum <- n_general_sum + n_dropout_sum + n_dropin_sum + n_ignore_sum
    id_mk_dedup1 <- unlist(mapply(rep, rep(1:n_mk, 4), c(n_general_use, n_dropout_use, n_dropin_use, n_ignore_use)))
    at_dropout <- unlist(mapply(rep, ats_use, n_dropout_use))
    loglike_di <- unlist(loglike_di_use)
    if(is.null(loglike_di)){
      loglike_di <- numeric(0)
    }

    # Arrange de-duplicated objects (pattern 1)
    index_dedup1_use <- matrix(NA, nrow = 2 * hnc, ncol = n_dedup1_sum)
    size_dedup1_use <- rep(NA, n_dedup1_sum)
    oph_dedup1_use <- rep(NA, n_dedup1_sum)
    index_general_start <- 1
    index_dropout_start <- n_general_sum + 1
    index_dropin_start <- n_general_sum + n_dropout_sum + 1
    index_ignore_start <- n_general_sum + n_dropout_sum + n_dropin_sum + 1
    sumup_index_dedup1 <- 0
    for(i in 1:n_mk){
      n_general <- n_general_use[i]
      if(n_general != 0){
        index_tmp <- index_general_start:(index_general_start + n_general - 1)
        index_dedup1_use[, index_tmp] <- index_dedup1_general_use[[i]] + sumup_index_dedup1
        size_dedup1_use[index_tmp] <- size_dedup1_general_use[[i]]
        oph_dedup1_use[index_tmp] <- oph_dedup1_general_use[[i]]
      }

      n_dropout <- n_dropout_use[i]
      if(n_dropout != 0){
        index_tmp <- index_dropout_start:(index_dropout_start + n_dropout - 1)
        index_dedup1_use[, index_tmp] <- index_dedup1_dropout_use[[i]] + sumup_index_dedup1
        size_dedup1_use[index_tmp] <- size_dedup1_dropout_use[[i]]
        oph_dedup1_use[index_tmp] <- oph_dedup1_dropout_use[[i]]
      }

      n_dropin <- n_dropin_use[i]
      if(n_dropin != 0){
        index_tmp <- index_dropin_start:(index_dropin_start + n_dropin - 1)
        index_dedup1_use[, index_tmp] <- index_dedup1_dropin_use[[i]] + sumup_index_dedup1
        size_dedup1_use[index_tmp] <- size_dedup1_dropin_use[[i]]
        oph_dedup1_use[index_tmp] <- oph_dedup1_dropin_use[[i]]
      }

      n_ignore <- n_ignore_use[i]
      if(n_ignore != 0){
        index_tmp <- index_ignore_start:(index_ignore_start + n_ignore - 1)
        index_dedup1_use[, index_tmp] <- index_dedup1_ignore_use[[i]] + sumup_index_dedup1
        size_dedup1_use[index_tmp] <- size_dedup1_ignore_use[[i]]
        oph_dedup1_use[index_tmp] <- oph_dedup1_ignore_use[[i]]
      }

      index_general_start <- index_general_start + n_general
      index_dropout_start <- index_dropout_start + n_dropout
      index_dropin_start <- index_dropin_start + n_dropin
      index_ignore_start <- index_ignore_start + n_ignore
      sumup_index_dedup1 <- sumup_index_dedup1 + n_rph_base_use[i]
    }

    # Define de-duplicated objects (pattern 2)
    n_dedup2_use <- n_peak_use * n_gtc_use
    n_dedup2_sum <- sum(n_dedup2_use)
    index_dedup2_use <- rep(NA, n_dedup2_sum)
    index_dedup2_start <- 1
    sumup_index_dedup2 <- 0
    for(i in 1:n_mk){
      index_rph_gtc <- index_rph_gtc_use[[i]] + sumup_index_dedup2
      n_dedup2 <- n_dedup2_use[i]
      index_dedup2_use[index_dedup2_start:(index_dedup2_start + n_dedup2 - 1)] <- make_index_dedup2_mk(index_rph_gtc, index_dedup1_use)

      index_dedup2_start <- index_dedup2_start + n_dedup2
      sumup_index_dedup2 <- sumup_index_dedup2 + n_rph_base_use[i]
    }

    data_list_main <- list(hnc = hnc, n_mk = n_mk, ltph = max_tph / 4, htph = max_tph / 2,
                           n_rph_base_sum = sum(n_rph_base_use), rph_base_use = unlist(rph_base_use), ae_alpha_use = ae_alpha_use, cv = cv,
                           n_dedup1_sum = n_dedup1_sum, index_dedup1_use = index_dedup1_use, size_dedup1_use = size_dedup1_use - min_size, oph_dedup1_use = oph_dedup1_use, id_mk_dedup1 = id_mk_dedup1,
                           n_general_sum = n_general_sum, n_dropout_sum = n_dropout_sum, at_dropout = at_dropout, n_dropin_sum = n_dropin_sum, loglike_di = loglike_di,
                           n_dedup2_sum = n_dedup2_sum, n_dedup2_use = n_dedup2_use, index_dedup2_use = index_dedup2_use,
                           n_peak_use = n_peak_use, n_gtc_use = n_gtc_use, gtc_logprob_use = unlist(gtc_logprob_use))

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


    index_dedup1_use_cpp <- asplit(index_dedup1_use - 1, 1) # C++ index
    index_dedup2_use_cpp <- index_dedup2_use - 1 # C++ index
    ae_dedup1_all_mk <- ae[id_mk_dedup1]
    lpd_gtc_use <- calc_lpd_gtc_use_cpp(n_mk, mp, deg, hnc, n_dedup1_sum, n_general_sum, n_dropout_sum, n_dropin_sum,
                                        unlist(rph_base_use), index_dedup1_use_cpp, unlist(size_dedup1_use), ae_dedup1_all_mk, oph_dedup1_use,
                                        at_dropout, loglike_di, cv, tph,
                                        index_dedup2_use_cpp, n_dedup2_use, n_peak_use, n_gtc_use)

    dt_pg_list <- make_dt_pg_list(markers_use, gtc_use, lpd_gtc_use, gtc_logprob_use, paste0("U", 1:hnc))

    log10like_all_mk <- calc_log10like_all_mk(lpd_gtc_use, gtc_logprob_use)

    dt_log10like <- make_dt_log10like(markers_use, log10like_all_mk)

    dt_mpdeg <- data.table(Contributor = paste0("U", 1:hnc),
                           MP = mp[hnc:1],
                           Deg = deg[hnc:1])

    dt_mcmc <- data.table(Parameter = c("Log probability in stan", "Total peak height", rep("lambda", hnc), rep("Mixture proportion", hnc), rep("Degradation", hnc), rep("Amplification efficiency", n_mk)),
                          Contributor = c("-", "-", paste0("U", 1:hnc), paste0("U", 1:hnc), paste0("U", 1:hnc), rep("-", n_mk)),
                          Marker = c("-", "-", rep("-", 3 * hnc), markers_use),
                          Estimated_value = c(lp, tph, lambda[hnc:1], mp[hnc:1], deg[hnc:1], ae),
                          Rhat = c(rhat_lp, rhat_tph, rhat_lambda[hnc:1], rhat_mp[hnc:1], rhat_deg[hnc:1], rhat_pae))

    # End computation
    end_time <- proc.time()
    comp_time <- end_time[3] - start_time[3]

    return(list("hnc" = hnc,
                "markers_use" = markers_use,
                "pop_al_list" = pop_al_list,
                "pop_prob_list" = pop_prob_list,
                "theta" = theta,
                "gtc_use" = gtc_use,
                "lpd_gtc_use" = lpd_gtc_use,
                "dt_pg_list" = dt_pg_list,
                "dt_log10like" = dt_log10like,
                "dt_mpdeg" = dt_mpdeg,
                "dt_mcmc" = dt_mcmc,
                "traceplot" = traceplot,
                "comp_time" = comp_time))
  }else{
    return(list("hnc" = NULL,
                "markers_use" = NULL,
                "pop_al_list" = NULL,
                "pop_prob_list" = NULL,
                "theta" = NULL,
                "gtc_use" = NULL,
                "lpd_gtc_use" = NULL,
                "dt_pg_list" = NULL,
                "dt_log10like" = NULL,
                "dt_mpdeg" = NULL,
                "dt_mcmc" = NULL,
                "traceplot" = NULL,
                "comp_time" = NULL))
  }
}
