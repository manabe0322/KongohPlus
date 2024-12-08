#' make_index_dedup2_mk
#'
#' @description The function to make 'index_dedup2_mk'
#' @param index_rph_gtc Information on the connection between rph_base and genotype combinations
#' @param index_dedup1_all_mk Index for deduplication (pattern 1)
make_index_dedup2_mk <- function(index_rph_gtc, index_dedup1_all_mk){
  paste_index_rph_gtc <- apply(index_rph_gtc, 2, paste, collapse = "-")
  paste_index_dedup1_all_mk <- apply(index_dedup1_all_mk, 2, paste, collapse = "-")
  return(match(paste_index_rph_gtc, paste_index_dedup1_all_mk))
}

#' perform_deconvo
#'
#' @description The function to perform mixture deconvolution
#' @param dt_csp The data.table of the crime stain profiles
#' @param dt_af The data.table of the allele frequencies
#' @param dt_kit The data.table of the information on the kit
#' @param dt_ignore_mk The data.table of the information on the ignored markers
#' @param dt_at The data.table of the analytical thresholds
#' @param dt_al_prob The data.table of the setting of the allele probabilities
#' @param dt_gt_fltr The data.table of the genotype filters
#' @param dt_dropin The data.table of the setting of the drop-in
#' @param dt_mcmc The data.table of the MCMC parameters
#' @param alpha_auto A vector of the estimated alpha values for Dirichlet distribution of the locus-specific amplification efficiencies
#' @param cv A parameter of the gamma distribution
#' @param list_dt_average_sr A list of the data.table for the average stutter ratios
#' @param hnc The hypothesized number of contributors
perform_deconvo <- function(dt_csp, dt_af, dt_kit, dt_ignore_mk, dt_at, dt_al_prob, dt_gt_fltr, dt_dropin, dt_mcmc,
                            alpha_auto, cv, list_dt_average_sr,
                            hnc){
  # Start computation
  start_time <- proc.time()

  # Number of replicates
  n_rep <- length(unique(dt_csp$SampleName))

  # Extract used markers
  markers_use <- dt_ignore_mk$Marker[!dt_ignore_mk$Ignore]
  n_mk <- length(markers_use)

  if(n_mk != 0){
    dt_kit_auto <- extract_kit_auto(dt_kit)

    # Extract ATs of the used markers
    ats_use <- as.numeric(dt_at$AT[match(markers_use, dt_at$Marker)])

    # Extract peak data from crime stain profiles
    tmp <- make_peak_data(dt_csp, dt_kit_auto, markers_use, ats_use)
    peaks_list <- tmp$peaks_list
    sizes_list <- tmp$sizes_list
    heights_list <- tmp$heights_list
    mean_size <- mean(dt_kit_auto$Size)

    # Extract parameters for AE
    ae_alpha_use <- alpha_auto[match(markers_use, names(alpha_auto))]

    # Extract theta
    theta <- as.numeric(dt_al_prob$Value[dt_al_prob$Subcategory == "Theta"])

    # Set allele probabilities
    tmp <- set_al_prob(markers_use, dt_af, peaks_list)
    pop_al_list <- tmp$pop_al_list
    pop_prob_list <- tmp$pop_prob_list

    # Extract filters to remove unrealistic genotype combinations
    gt_fltr_subcategory <- dt_gt_fltr$Subcategory
    gt_fltr_value <- as.numeric(dt_gt_fltr$Value)
    bsr_fltr <- gt_fltr_value[gt_fltr_subcategory == "Maximum BSR"]
    fsr_fltr <- gt_fltr_value[gt_fltr_subcategory == "Maximum FSR"]
    dsr_fltr <- gt_fltr_value[gt_fltr_subcategory == "Maximum DSR"]
    m2sr_fltr <- gt_fltr_value[gt_fltr_subcategory == "Maximum M2SR"]
    n_dropin_fltr <- gt_fltr_value[gt_fltr_subcategory == "Maximum number of drop-in"]
    h_dropin_fltr <- gt_fltr_value[gt_fltr_subcategory == "Maximum height of drop-in"]

    # Extract drop-in parameters
    dropin_subcategory <- dt_dropin$Subcategory
    dropin_value <- as.numeric(dt_dropin$Value)
    logprob_di <- log(dropin_value[dropin_subcategory == "Probability of drop-in"])
    shape_di <- dropin_value[dropin_subcategory == "Shape parameter"]
    rate_di <- dropin_value[dropin_subcategory == "Rate parameter"]
    mt_di <- dropin_value[dropin_subcategory == "Minimum threshold"]

    # Extract MCMC parameters
    mcmc_subcategory <- dt_mcmc$Subcategory
    mcmc_value <- as.numeric(dt_mcmc$Value)
    warmup <- mcmc_value[mcmc_subcategory == "Warmup"]
    sampling <- mcmc_value[mcmc_subcategory == "Sampling"]
    chains <- mcmc_value[mcmc_subcategory == "Chains"]

    # Define objects
    rph_base_use <- vector(mode = "list", length = n_mk)
    gtc_logprob_use <- vector(mode = "list", length = n_mk)
    index_rph_gtc_use <- vector(mode = "list", length = n_rep)
    index_dedup1_general_use <- vector(mode = "list", length = n_rep)
    index_dedup1_dropout_use <- vector(mode = "list", length = n_rep)
    index_dedup1_dropin_use <- vector(mode = "list", length = n_rep)
    index_dedup1_ignore_use <- vector(mode = "list", length = n_rep)
    size_dedup1_general_use <- vector(mode = "list", length = n_rep)
    size_dedup1_dropout_use <- vector(mode = "list", length = n_rep)
    size_dedup1_dropin_use <- vector(mode = "list", length = n_rep)
    size_dedup1_ignore_use <- vector(mode = "list", length = n_rep)
    oph_dedup1_general_use <- vector(mode = "list", length = n_rep)
    oph_dedup1_dropout_use <- vector(mode = "list", length = n_rep)
    oph_dedup1_dropin_use <- vector(mode = "list", length = n_rep)
    oph_dedup1_ignore_use <- vector(mode = "list", length = n_rep)
    loglike_di_use <- vector(mode = "list", length = n_rep)
    for(i in 1:n_rep){
      index_rph_gtc_use[[i]] <- vector(mode = "list", length = n_mk)
      index_dedup1_general_use[[i]] <- vector(mode = "list", length = n_mk)
      index_dedup1_dropout_use[[i]] <- vector(mode = "list", length = n_mk)
      index_dedup1_dropin_use[[i]] <- vector(mode = "list", length = n_mk)
      index_dedup1_ignore_use[[i]] <- vector(mode = "list", length = n_mk)
      size_dedup1_general_use[[i]] <- vector(mode = "list", length = n_mk)
      size_dedup1_dropout_use[[i]] <- vector(mode = "list", length = n_mk)
      size_dedup1_dropin_use[[i]] <- vector(mode = "list", length = n_mk)
      size_dedup1_ignore_use[[i]] <- vector(mode = "list", length = n_mk)
      oph_dedup1_general_use[[i]] <- vector(mode = "list", length = n_mk)
      oph_dedup1_dropout_use[[i]] <- vector(mode = "list", length = n_mk)
      oph_dedup1_dropin_use[[i]] <- vector(mode = "list", length = n_mk)
      oph_dedup1_ignore_use[[i]] <- vector(mode = "list", length = n_mk)
      loglike_di_use[[i]] <- vector(mode = "list", length = n_mk)
    }

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
      n_peak_use[i] <- n_peak <- length(peaks)

      # Extract parameters of one marker
      at <- ats_use[i]
      dt_ave_sr <- list_dt_average_sr[[mk]]
      sr_alleles_mk <- as.numeric(dt_ave_sr$Allele)
      bsrs_mk <- dt_ave_sr[["Average BSR"]]
      fsrs_mk <- dt_ave_sr[["Average FSR"]]
      dsrs_mk <- dt_ave_sr[["Average DSR"]]
      m2srs_mk <- dt_ave_sr[["Average M2SR"]]
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

      ###########################################
      # List all possible genotype combinations #
      ###########################################

      # All candidate genotypes of one contributor
      cand_gt <- combinations(n = n_peak, r = 2, v = peaks, repeats.allowed = TRUE)
      cand_gt <- asplit(cand_gt, 1)

      # Indices of genotypes of each contributor in all candidate genotype combinations
      index_gtc <- permutations(n = length(cand_gt), r = hnc, repeats.allowed = TRUE)
      index_gtc_cpp <- index_gtc - 1 # C++ index
      index_gtc_cpp <- asplit(index_gtc_cpp, 1)

      # Indices of genotypes of each contributor in each adopted genotype combination
      tmp <- judge_gtc_mk(cand_gt, index_gtc_cpp, peaks, heights_all_rep,
                          bsrs, fsrs, dsrs, m2srs,
                          bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr, n_dropin_fltr, h_dropin_fltr)
      index_gtc_adopt_cpp <- tmp[[1]]
      bool_gtc_di_mk_rep <- tmp[[2]]

      n_gtc_use[i] <- n_gtc <- length(index_gtc_adopt_cpp)

      ##################################################
      # Prepare objects for de-duplication (pattern 1) #
      ##################################################

      if(!any(is.element(n_gtc_use, 0))){
        # Extract genotype combinations
        gtc_use[[i]] <- gtc_mk <- extract_gtc_adopt(cand_gt, index_gtc_adopt_cpp)

        # Calculate genotype probabilities
        n_obs <- make_n_obs(pop_al, prof_ref = numeric(0))
        gtc_prob_mk <- calc_gtc_prob_mk(gtc_mk, pop_al, pop_prob, theta, prof_ref = numeric(0), hyp_id = rep(0, hnc), n_obs)
        gtc_logprob_use[[i]] <- log(gtc_prob_mk)

        for(j in 1:n_rep){
          # Connect rph_base and genotype combinations
          index_rph_gtc_cpp <- connect_rph_gtc_mk(cand_gt, index_gtc_adopt_cpp, peaks, bsrs, fsrs, dsrs, m2srs, bool_gtc_di_mk_rep[[j]])
          index_rph_gtc_use[[j]][[i]] <- index_rph_gtc <- matrix(unlist(index_rph_gtc_cpp) + 1, nrow = 2 * hnc) # R and stan index

          # Delete duplicated indices -> index for deduplication (pattern 1: Before summation of the rph_base)
          bool_duplicated_adopt <- !duplicated(t(index_rph_gtc))
          index_dedup1_tmp <- index_rph_gtc[, bool_duplicated_adopt, drop = FALSE]

          # Sizes for de-duplicated objects (pattern 1)
          size_dedup1_tmp <- rep(sizes, n_gtc)
          size_dedup1_tmp <- size_dedup1_tmp[bool_duplicated_adopt]

          # Update the maximum of the total peak height per marker
          heights <- heights_all_rep[[j]]
          oph_sum <- sum(heights)
          if(oph_sum > max_tph){
            max_tph <- oph_sum
          }

          # Observe peak heights for de-duplicated objects (pattern 1)
          oph_dedup1_tmp <- rep(heights, n_gtc)
          oph_dedup1_tmp <- oph_dedup1_tmp[bool_duplicated_adopt]

          # Indices of each peak type
          rph_dedup1_tmp <- matrix(rph_base_mk[index_dedup1_tmp], nrow = 2 * hnc)
          rph_dedup1_tmp <- apply(rph_dedup1_tmp, 2, sum)
          index_oph0 <- which(oph_dedup1_tmp == 0)
          index_ignore <- c(which(rph_dedup1_tmp == 0), intersect(which(rph_dedup1_tmp < 0), index_oph0)) # intersect(which(rph_dedup1_tmp < 0), index_oph0): for duplicate
          index_dropin <- setdiff(which(rph_dedup1_tmp < 0), index_ignore)
          index_dropout <- setdiff(index_oph0, c(index_ignore, index_dropin))
          index_general <- setdiff(1:length(rph_dedup1_tmp), c(index_ignore, index_dropin, index_dropout))

          index_dedup1_general_use[[j]][[i]] <- index_dedup1_tmp[, index_general, drop = FALSE]
          index_dedup1_dropout_use[[j]][[i]] <- index_dedup1_tmp[, index_dropout, drop = FALSE]
          index_dedup1_dropin_use[[j]][[i]] <- index_dedup1_tmp[, index_dropin, drop = FALSE]
          index_dedup1_ignore_use[[j]][[i]] <- index_dedup1_tmp[, index_ignore, drop = FALSE]

          size_dedup1_general_use[[j]][[i]] <- size_dedup1_tmp[index_general]
          size_dedup1_dropout_use[[j]][[i]] <- size_dedup1_tmp[index_dropout]
          size_dedup1_dropin_use[[j]][[i]] <- size_dedup1_tmp[index_dropin]
          size_dedup1_ignore_use[[j]][[i]] <- size_dedup1_tmp[index_ignore]

          oph_dedup1_general_use[[j]][[i]] <- oph_dedup1_tmp[index_general]
          oph_dedup1_dropout_use[[j]][[i]] <- oph_dedup1_tmp[index_dropout]
          oph_dedup1_dropin_use[[j]][[i]] <- oph_dedup1_tmp[index_dropin]
          oph_dedup1_ignore_use[[j]][[i]] <- oph_dedup1_tmp[index_ignore]

          if(length(index_dropin) != 0){
            loglike_di_use[[j]][[i]] <- log(dgamma(oph_dedup1_tmp[index_dropin] - mt_di + 1, shape_di, rate_di)) + logprob_di
          }
        }
      }
    }

    if(all(n_gtc_use != 0)){

      ##################################
      # Define objects for replicate 1 #
      ##################################

      n_rph_base_use <- sapply(rph_base_use, length)
      n_general_use_1st <- sapply(index_dedup1_general_use[[1]], ncol)
      n_dropout_use_1st <- sapply(index_dedup1_dropout_use[[1]], ncol)
      n_dropin_use_1st <- sapply(index_dedup1_dropin_use[[1]], ncol)
      n_ignore_use_1st <- sapply(index_dedup1_ignore_use[[1]], ncol)
      n_general_sum_1st <- sum(n_general_use_1st)
      n_dropout_sum_1st <- sum(n_dropout_use_1st)
      n_dropin_sum_1st <- sum(n_dropin_use_1st)
      n_ignore_sum_1st <- sum(n_ignore_use_1st)
      n_dedup1_sum_1st <- n_general_sum_1st + n_dropout_sum_1st + n_dropin_sum_1st + n_ignore_sum_1st
      id_mk_dedup1_1st <- unlist(mapply(rep, rep(1:n_mk, 4), c(n_general_use_1st, n_dropout_use_1st, n_dropin_use_1st, n_ignore_use_1st)))
      at_dropout_1st <- unlist(mapply(rep, ats_use, n_dropout_use_1st))
      loglike_di_1st <- unlist(loglike_di_use[[1]])
      if(is.null(loglike_di_1st)){
        loglike_di_1st <- numeric(0)
      }

      ##################################################
      # Arrange objects for de-duplication (pattern 1) #
      ##################################################

      index_dedup1_use_1st <- matrix(NA, nrow = 2 * hnc, ncol = n_dedup1_sum_1st)
      size_dedup1_use_1st <- rep(NA, n_dedup1_sum_1st)
      oph_dedup1_use_1st <- rep(NA, n_dedup1_sum_1st)
      index_general_start <- 1
      index_dropout_start <- n_general_sum_1st + 1
      index_dropin_start <- n_general_sum_1st + n_dropout_sum_1st + 1
      index_ignore_start <- n_general_sum_1st + n_dropout_sum_1st + n_dropin_sum_1st + 1
      sumup_index_dedup1 <- 0
      for(i in 1:n_mk){
        n_general <- n_general_use_1st[i]
        if(n_general != 0){
          index_tmp <- index_general_start:(index_general_start + n_general - 1)
          index_dedup1_use_1st[, index_tmp] <- index_dedup1_general_use[[1]][[i]] + sumup_index_dedup1
          size_dedup1_use_1st[index_tmp] <- size_dedup1_general_use[[1]][[i]]
          oph_dedup1_use_1st[index_tmp] <- oph_dedup1_general_use[[1]][[i]]
        }

        n_dropout <- n_dropout_use_1st[i]
        if(n_dropout != 0){
          index_tmp <- index_dropout_start:(index_dropout_start + n_dropout - 1)
          index_dedup1_use_1st[, index_tmp] <- index_dedup1_dropout_use[[1]][[i]] + sumup_index_dedup1
          size_dedup1_use_1st[index_tmp] <- size_dedup1_dropout_use[[1]][[i]]
          oph_dedup1_use_1st[index_tmp] <- oph_dedup1_dropout_use[[1]][[i]]
        }

        n_dropin <- n_dropin_use_1st[i]
        if(n_dropin != 0){
          index_tmp <- index_dropin_start:(index_dropin_start + n_dropin - 1)
          index_dedup1_use_1st[, index_tmp] <- index_dedup1_dropin_use[[1]][[i]] + sumup_index_dedup1
          size_dedup1_use_1st[index_tmp] <- size_dedup1_dropin_use[[1]][[i]]
          oph_dedup1_use_1st[index_tmp] <- oph_dedup1_dropin_use[[1]][[i]]
        }

        n_ignore <- n_ignore_use_1st[i]
        if(n_ignore != 0){
          index_tmp <- index_ignore_start:(index_ignore_start + n_ignore - 1)
          index_dedup1_use_1st[, index_tmp] <- index_dedup1_ignore_use[[1]][[i]] + sumup_index_dedup1
          size_dedup1_use_1st[index_tmp] <- size_dedup1_ignore_use[[1]][[i]]
          oph_dedup1_use_1st[index_tmp] <- oph_dedup1_ignore_use[[1]][[i]]
        }

        index_general_start <- index_general_start + n_general
        index_dropout_start <- index_dropout_start + n_dropout
        index_dropin_start <- index_dropin_start + n_dropin
        index_ignore_start <- index_ignore_start + n_ignore
        sumup_index_dedup1 <- sumup_index_dedup1 + n_rph_base_use[i]
      }

      ##################################################
      # Prepare objects for de-duplication (pattern 2) #
      ##################################################

      n_dedup2_use <- n_peak_use * n_gtc_use
      n_dedup2_sum <- sum(n_dedup2_use)
      index_dedup2_use_1st <- rep(NA, n_dedup2_sum)
      index_dedup2_start <- 1
      sumup_index_dedup2 <- 0
      for(i in 1:n_mk){
        index_rph_gtc <- index_rph_gtc_use[[1]][[i]] + sumup_index_dedup2
        n_dedup2 <- n_dedup2_use[i]
        index_dedup2_use_1st[index_dedup2_start:(index_dedup2_start + n_dedup2 - 1)] <- make_index_dedup2_mk(index_rph_gtc, index_dedup1_use_1st)

        index_dedup2_start <- index_dedup2_start + n_dedup2
        sumup_index_dedup2 <- sumup_index_dedup2 + n_rph_base_use[i]
      }

      ########################
      # MCMC for one profile #
      ########################

      if(n_rep == 1){

        # Prepare date list
        data_list_main <- list(hnc = hnc,
                               n_mk = n_mk,
                               ltph = max_tph / 4,
                               htph = max_tph / 2,
                               n_rph_base_sum = sum(n_rph_base_use),
                               rph_base_use = unlist(rph_base_use),
                               ae_alpha_use = ae_alpha_use,
                               cv = cv,
                               n_dedup1_sum = n_dedup1_sum_1st,
                               index_dedup1_use = index_dedup1_use_1st,
                               size_dedup1_use = size_dedup1_use_1st,
                               mean_size = mean_size,
                               oph_dedup1_use = oph_dedup1_use_1st,
                               id_mk_dedup1 = id_mk_dedup1_1st,
                               n_general_sum = n_general_sum_1st,
                               n_dropout_sum = n_dropout_sum_1st,
                               at_dropout = at_dropout_1st,
                               n_dropin_sum = n_dropin_sum_1st,
                               loglike_di = loglike_di_1st,
                               n_dedup2_sum = n_dedup2_sum,
                               n_dedup2_use = n_dedup2_use,
                               index_dedup2_use = index_dedup2_use_1st,
                               n_peak_use = n_peak_use,
                               n_gtc_use = n_gtc_use,
                               gtc_logprob_use = unlist(gtc_logprob_use))

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

        # Calculate likelihood under the estimated parameter
        index_dedup1_use_cpp <- asplit(index_dedup1_use_1st - 1, 1) # C++ index
        index_dedup2_use_cpp <- index_dedup2_use_1st - 1 # C++ index
        ae_dedup1_use <- ae[id_mk_dedup1_1st]
        lpd_gtc_use_1st <- calc_lpd_gtc_all_mk_cpp(n_mk, mp, deg, hnc, n_dedup1_sum_1st, n_general_sum_1st, n_dropout_sum_1st, n_dropin_sum_1st,
                                                   unlist(rph_base_use), index_dedup1_use_cpp, unlist(size_dedup1_use_1st), mean_size, ae_dedup1_use, oph_dedup1_use_1st,
                                                   at_dropout_1st, loglike_di_1st, cv, tph,
                                                   index_dedup2_use_cpp, n_dedup2_use, n_peak_use, n_gtc_use)
        lpd_gtc_use_2nd <- NULL
        log10like_use <- calc_log10like_use(gtc_logprob_use, lpd_gtc_use_1st)

        dt_pg_list <- make_dt_pg_list(markers_use, paste0("U", 1:hnc), gtc_use, gtc_logprob_use, lpd_gtc_use_1st)

        dt_diag <- data.table(Parameter = c("Log probability in stan", "Total peak height", rep("lambda", hnc), rep("Mixture proportion", hnc), rep("Degradation", hnc), rep("Amplification efficiency", n_mk)),
                              Contributor = c("-", "-", paste0("U", 1:hnc), paste0("U", 1:hnc), paste0("U", 1:hnc), rep("-", n_mk)),
                              Marker = c("-", "-", rep("-", 3 * hnc), markers_use),
                              Estimated_value = c(lp, tph, lambda[hnc:1], mp[hnc:1], deg[hnc:1], ae),
                              Rhat = c(rhat_lp, rhat_tph, rhat_lambda[hnc:1], rhat_mp[hnc:1], rhat_deg[hnc:1], rhat_pae))

      #######################
      # Duplicated profiles #
      #######################

      }else{

        ##################################
        # Define objects for replicate 2 #
        ##################################

        n_general_use_2nd <- sapply(index_dedup1_general_use[[2]], ncol)
        n_dropout_use_2nd <- sapply(index_dedup1_dropout_use[[2]], ncol)
        n_dropin_use_2nd <- sapply(index_dedup1_dropin_use[[2]], ncol)
        n_ignore_use_2nd <- sapply(index_dedup1_ignore_use[[2]], ncol)
        n_general_sum_2nd <- sum(n_general_use_2nd)
        n_dropout_sum_2nd <- sum(n_dropout_use_2nd)
        n_dropin_sum_2nd <- sum(n_dropin_use_2nd)
        n_ignore_sum_2nd <- sum(n_ignore_use_2nd)
        n_dedup1_sum_2nd <- n_general_sum_2nd + n_dropout_sum_2nd + n_dropin_sum_2nd + n_ignore_sum_2nd
        id_mk_dedup1_2nd <- unlist(mapply(rep, rep(1:n_mk, 4), c(n_general_use_2nd, n_dropout_use_2nd, n_dropin_use_2nd, n_ignore_use_2nd)))
        at_dropout_2nd <- unlist(mapply(rep, ats_use, n_dropout_use_2nd))
        loglike_di_2nd <- unlist(loglike_di_use[[2]])
        if(is.null(loglike_di_2nd)){
          loglike_di_2nd <- numeric(0)
        }

        ##################################################
        # Arrange objects for de-duplication (pattern 1) #
        ##################################################

        index_dedup1_use_2nd <- matrix(NA, nrow = 2 * hnc, ncol = n_dedup1_sum_2nd)
        size_dedup1_use_2nd <- rep(NA, n_dedup1_sum_2nd)
        oph_dedup1_use_2nd <- rep(NA, n_dedup1_sum_2nd)
        index_general_start <- 1
        index_dropout_start <- n_general_sum_2nd + 1
        index_dropin_start <- n_general_sum_2nd + n_dropout_sum_2nd + 1
        index_ignore_start <- n_general_sum_2nd + n_dropout_sum_2nd + n_dropin_sum_2nd + 1
        sumup_index_dedup1 <- 0
        for(i in 1:n_mk){
          n_general <- n_general_use_2nd[i]
          if(n_general != 0){
            index_tmp <- index_general_start:(index_general_start + n_general - 1)
            index_dedup1_use_2nd[, index_tmp] <- index_dedup1_general_use[[2]][[i]] + sumup_index_dedup1
            size_dedup1_use_2nd[index_tmp] <- size_dedup1_general_use[[2]][[i]]
            oph_dedup1_use_2nd[index_tmp] <- oph_dedup1_general_use[[2]][[i]]
          }

          n_dropout <- n_dropout_use_2nd[i]
          if(n_dropout != 0){
            index_tmp <- index_dropout_start:(index_dropout_start + n_dropout - 1)
            index_dedup1_use_2nd[, index_tmp] <- index_dedup1_dropout_use[[2]][[i]] + sumup_index_dedup1
            size_dedup1_use_2nd[index_tmp] <- size_dedup1_dropout_use[[2]][[i]]
            oph_dedup1_use_2nd[index_tmp] <- oph_dedup1_dropout_use[[2]][[i]]
          }

          n_dropin <- n_dropin_use_2nd[i]
          if(n_dropin != 0){
            index_tmp <- index_dropin_start:(index_dropin_start + n_dropin - 1)
            index_dedup1_use_2nd[, index_tmp] <- index_dedup1_dropin_use[[2]][[i]] + sumup_index_dedup1
            size_dedup1_use_2nd[index_tmp] <- size_dedup1_dropin_use[[2]][[i]]
            oph_dedup1_use_2nd[index_tmp] <- oph_dedup1_dropin_use[[2]][[i]]
          }

          n_ignore <- n_ignore_use_2nd[i]
          if(n_ignore != 0){
            index_tmp <- index_ignore_start:(index_ignore_start + n_ignore - 1)
            index_dedup1_use_2nd[, index_tmp] <- index_dedup1_ignore_use[[2]][[i]] + sumup_index_dedup1
            size_dedup1_use_2nd[index_tmp] <- size_dedup1_ignore_use[[2]][[i]]
            oph_dedup1_use_2nd[index_tmp] <- oph_dedup1_ignore_use[[2]][[i]]
          }

          index_general_start <- index_general_start + n_general
          index_dropout_start <- index_dropout_start + n_dropout
          index_dropin_start <- index_dropin_start + n_dropin
          index_ignore_start <- index_ignore_start + n_ignore
          sumup_index_dedup1 <- sumup_index_dedup1 + n_rph_base_use[i]
        }

        ##################################################
        # Prepare objects for de-duplication (pattern 1) #
        ##################################################

        n_dedup2_use <- n_peak_use * n_gtc_use
        n_dedup2_sum <- sum(n_dedup2_use)
        index_dedup2_use_2nd <- rep(NA, n_dedup2_sum)
        index_dedup2_start <- 1
        sumup_index_dedup2 <- 0
        for(i in 1:n_mk){
          index_rph_gtc <- index_rph_gtc_use[[2]][[i]] + sumup_index_dedup2
          n_dedup2 <- n_dedup2_use[i]
          index_dedup2_use_2nd[index_dedup2_start:(index_dedup2_start + n_dedup2 - 1)] <- make_index_dedup2_mk(index_rph_gtc, index_dedup1_use_2nd)

          index_dedup2_start <- index_dedup2_start + n_dedup2
          sumup_index_dedup2 <- sumup_index_dedup2 + n_rph_base_use[i]
        }

        ################################
        # MCMC for duplicated profiles #
        ################################

        # Prepare date list
        data_list_main <- list(hnc = hnc,
                               n_mk = n_mk,
                               ltph = max_tph / 4,
                               htph = max_tph / 2,
                               n_rph_base_sum = sum(n_rph_base_use),
                               rph_base_use = unlist(rph_base_use),
                               ae_alpha_use = ae_alpha_use,
                               cv = cv,
                               mean_size = mean_size,
                               n_dedup1_sum_1st = n_dedup1_sum_1st,
                               index_dedup1_use_1st = index_dedup1_use_1st,
                               size_dedup1_use_1st = size_dedup1_use_1st,
                               oph_dedup1_use_1st = oph_dedup1_use_1st,
                               id_mk_dedup1_1st = id_mk_dedup1_1st,
                               n_general_sum_1st = n_general_sum_1st,
                               n_dropout_sum_1st = n_dropout_sum_1st,
                               at_dropout_1st = at_dropout_1st,
                               n_dropin_sum_1st = n_dropin_sum_1st,
                               loglike_di_1st = loglike_di_1st,
                               n_dedup1_sum_2nd = n_dedup1_sum_2nd,
                               index_dedup1_use_2nd = index_dedup1_use_2nd,
                               size_dedup1_use_2nd = size_dedup1_use_2nd,
                               oph_dedup1_use_2nd = oph_dedup1_use_2nd,
                               id_mk_dedup1_2nd = id_mk_dedup1_2nd,
                               n_general_sum_2nd = n_general_sum_2nd,
                               n_dropout_sum_2nd = n_dropout_sum_2nd,
                               at_dropout_2nd = at_dropout_2nd,
                               n_dropin_sum_2nd = n_dropin_sum_2nd,
                               loglike_di_2nd = loglike_di_2nd,
                               n_dedup2_sum = n_dedup2_sum,
                               n_dedup2_use = n_dedup2_use,
                               index_dedup2_use_1st = index_dedup2_use_1st,
                               index_dedup2_use_2nd = index_dedup2_use_2nd,
                               n_peak_use = n_peak_use,
                               n_gtc_use = n_gtc_use,
                               gtc_logprob_use = unlist(gtc_logprob_use))

        # Perform MCMC
        model_main_duplicate <- stan_package_model(name = "model_main_duplicate", package = "KongohPlus")
        fit_main <- model_main_duplicate$sample(
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
        tph_1st <- medians[which(variables == "tph_1st")]
        tph_2nd <- medians[which(variables == "tph_2nd")]
        lambda <- medians[grep("lambda", variables)]
        mp <- medians[grep("mp", variables)]
        deg <- medians[grep("deg", variables)]
        pae <- medians[grep("pae", variables)]
        ae <- n_mk * pae
        rhat_lp <- rhats[which(variables == "lp__")]
        rhat_tph_1st <- rhats[which(variables == "tph_1st")]
        rhat_tph_2nd <- rhats[which(variables == "tph_2nd")]
        rhat_lambda <- rhats[grep("lambda", variables)]
        rhat_mp <- rhats[grep("mp", variables)]
        rhat_deg <- rhats[grep("deg", variables)]
        rhat_pae <- rhats[grep("pae", variables)]

        # Calculate likelihood under the estimated parameter
        index_dedup1_use_cpp_1st <- asplit(index_dedup1_use_1st - 1, 1) # C++ index
        index_dedup2_use_cpp_1st <- index_dedup2_use_1st - 1 # C++ index
        ae_dedup1_use_1st <- ae[id_mk_dedup1_1st]
        lpd_gtc_use_1st <- calc_lpd_gtc_all_mk_cpp(n_mk, mp, deg, hnc, n_dedup1_sum_1st, n_general_sum_1st, n_dropout_sum_1st, n_dropin_sum_1st,
                                                   unlist(rph_base_use), index_dedup1_use_cpp_1st, unlist(size_dedup1_use_1st), mean_size, ae_dedup1_use_1st, oph_dedup1_use_1st,
                                                   at_dropout_1st, loglike_di_1st, cv, tph_1st,
                                                   index_dedup2_use_cpp_1st, n_dedup2_use, n_peak_use, n_gtc_use)
        index_dedup1_use_cpp_2nd <- asplit(index_dedup1_use_2nd - 1, 1) # C++ index
        index_dedup2_use_cpp_2nd <- index_dedup2_use_2nd - 1 # C++ index
        ae_dedup1_use_2nd <- ae[id_mk_dedup1_2nd]
        lpd_gtc_use_2nd <- calc_lpd_gtc_all_mk_cpp(n_mk, mp, deg, hnc, n_dedup1_sum_2nd, n_general_sum_2nd, n_dropout_sum_2nd, n_dropin_sum_2nd,
                                                   unlist(rph_base_use), index_dedup1_use_cpp_2nd, unlist(size_dedup1_use_2nd), mean_size, ae_dedup1_use_2nd, oph_dedup1_use_2nd,
                                                   at_dropout_2nd, loglike_di_2nd, cv, tph_2nd,
                                                   index_dedup2_use_cpp_2nd, n_dedup2_use, n_peak_use, n_gtc_use)
        log10like_use <- calc_log10like_use(gtc_logprob_use, lpd_gtc_use_1st, lpd_gtc_use_2nd)

        dt_pg_list <- make_dt_pg_list(markers_use, paste0("U", 1:hnc), gtc_use, gtc_logprob_use, lpd_gtc_use_1st, lpd_gtc_use_2nd)

        dt_diag <- data.table(Parameter = c("Log probability in stan", "Total peak height (1st)", "Total peak height (2nd)", rep("lambda", hnc), rep("Mixture proportion", hnc), rep("Degradation", hnc), rep("Amplification efficiency", n_mk)),
                              Contributor = c("-", "-", "-", paste0("U", 1:hnc), paste0("U", 1:hnc), paste0("U", 1:hnc), rep("-", n_mk)),
                              Marker = c("-", "-", "-", rep("-", 3 * hnc), markers_use),
                              Estimated_value = c(lp, tph_1st, tph_2nd, lambda[hnc:1], mp[hnc:1], deg[hnc:1], ae),
                              Rhat = c(rhat_lp, rhat_tph_1st, rhat_tph_2nd, rhat_lambda[hnc:1], rhat_mp[hnc:1], rhat_deg[hnc:1], rhat_pae))
      }

      dt_log10like <- data.table(Marker = c(markers_use, "Total"),
                                 Log10_Likelihood = c(log10like_use, sum(log10like_use)))

      dt_mpdeg <- data.table(Contributor = paste0("U", 1:hnc),
                             MP = mp[hnc:1],
                             Deg = deg[hnc:1])

      # End computation
      end_time <- proc.time()
      comp_time <- end_time[3] - start_time[3]

      return(list("markers_use" = markers_use,
                  "peaks_list" = peaks_list,
                  "pop_al_list" = pop_al_list,
                  "pop_prob_list" = pop_prob_list,
                  "gtc_use" = gtc_use,
                  "lpd_gtc_use_1st" = lpd_gtc_use_1st,
                  "lpd_gtc_use_2nd" = lpd_gtc_use_2nd,
                  "dt_pg_list" = dt_pg_list,
                  "dt_log10like" = dt_log10like,
                  "dt_mpdeg" = dt_mpdeg,
                  "dt_diag" = dt_diag,
                  "traceplot" = traceplot,
                  "comp_time" = comp_time,
                  "error_message" = NULL))
    }else{
      markers_unexplainable <- markers_use[n_gtc_use == 0]
      if(hnc == 1){
        error_message <- paste0("The following marker(s) cannot be explained under ", hnc, " contributor.<br>")
      }else{
        error_message <- paste0("The following marker(s) cannot be explained under ", hnc, " contributors.<br>")
      }
      for(i in 1:length(markers_unexplainable)){
        error_message <- paste0(error_message, "  - ", markers_unexplainable[i], "<br>")
      }

      return(list("markers_use" = NULL,
                  "peaks_list" = NULL,
                  "pop_al_list" = NULL,
                  "pop_prob_list" = NULL,
                  "gtc_use" = NULL,
                  "lpd_gtc_use_1st" = NULL,
                  "lpd_gtc_use_2nd" = NULL,
                  "dt_pg_list" = NULL,
                  "dt_log10like" = NULL,
                  "dt_mpdeg" = NULL,
                  "dt_diag" = NULL,
                  "traceplot" = NULL,
                  "comp_time" = NULL,
                  "error_message" = error_message))
    }
  }else{
    return(list("markers_use" = NULL,
                "peaks_list" = NULL,
                "pop_al_list" = NULL,
                "pop_prob_list" = NULL,
                "gtc_use" = NULL,
                "lpd_gtc_use_1st" = NULL,
                "lpd_gtc_use_2nd" = NULL,
                "dt_pg_list" = NULL,
                "dt_log10like" = NULL,
                "dt_mpdeg" = NULL,
                "dt_diag" = NULL,
                "traceplot" = NULL,
                "comp_time" = NULL,
                "error_message" = "All markers are ignored!"))
  }
}
