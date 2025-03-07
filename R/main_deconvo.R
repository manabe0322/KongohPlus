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
#' @param alpha_auto A vector of the estimated alpha values for Dirichlet distribution of the locus-specific amplification efficiencies
#' @param cv A parameter of the gamma distribution
#' @param list_dt_average_sr A list of the data.table for the average stutter ratios
#' @param hnc The number of contributors
perform_deconvo <- function(hnc, n_rep, markers_use, ats_use, peaks_list, sizes_list, heights_list, mean_size,
                            pop_al_list, pop_prob_list, theta,
                            bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr, n_dropin_fltr, h_dropin_fltr,
                            logprob_di, shape_di, rate_di,
                            warmup, sampling, chains,
                            ae_alpha_use, cv, list_dt_average_sr){
  # Start computation
  start_time <- proc.time()

  # Define objects
  n_mk <- length(markers_use)
  n_peak_use <- rep(NA, n_mk)
  n_gtc_use <- rep(NA, n_mk)
  gtc_use <- vector(mode = "list", length = n_mk)
  names(gtc_use) <- markers_use
  max_tph <- 0
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
      obs_al <- as.numeric(names(n_obs))
      gtc_prob_mk <- calc_gtc_prob_mk(gtc_mk, pop_al, pop_prob, theta, prof_ref = numeric(0), hyp_id = rep(0, hnc), n_obs, obs_al)
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
          loglike_di_use[[j]][[i]] <- log(dgamma(oph_dedup1_tmp[index_dropin], shape_di, rate_di)) + logprob_di
        }
      }
    }
  }

  if(all(n_gtc_use != 0)){

    ###########################
    # Define objects for MCMC #
    ###########################

    n_rph_base_use <- sapply(rph_base_use, length)

    n_general_use <- rep(NA, n_rep * n_mk)
    n_dropout_use <- rep(NA, n_rep * n_mk)
    n_dropin_use <- rep(NA, n_rep * n_mk)
    n_ignore_use <- rep(NA, n_rep * n_mk)
    n_general_sum_per_rep <- rep(NA, n_rep)
    n_dropout_sum_per_rep <- rep(NA, n_rep)
    n_dropin_sum_per_rep <- rep(NA, n_rep)
    n_ignore_sum_per_rep <- rep(NA, n_rep)
    n_dedup1_sum_per_rep <- rep(NA, n_rep)
    id_mk_dedup1 <- vector(mode = "list", length = n_rep)
    at_dropout <- vector(mode = "list", length = n_rep)
    for(i in 1:n_rep){
      n_general_use[(n_mk * (i - 1) + 1):(n_mk * i)] <- n_general_rep1 <- sapply(index_dedup1_general_use[[i]], ncol)
      n_dropout_use[(n_mk * (i - 1) + 1):(n_mk * i)] <- n_dropout_rep1 <- sapply(index_dedup1_dropout_use[[i]], ncol)
      n_dropin_use[(n_mk * (i - 1) + 1):(n_mk * i)] <- n_dropin_rep1 <- sapply(index_dedup1_dropin_use[[i]], ncol)
      n_ignore_use[(n_mk * (i - 1) + 1):(n_mk * i)] <- n_ignore_rep1 <- sapply(index_dedup1_ignore_use[[i]], ncol)
      n_general_sum_per_rep[i] <- sum(n_general_rep1)
      n_dropout_sum_per_rep[i] <- sum(n_dropout_rep1)
      n_dropin_sum_per_rep[i] <- sum(n_dropin_rep1)
      n_ignore_sum_per_rep[i] <- sum(n_ignore_rep1)
      n_dedup1_sum_per_rep[i] <- n_general_sum_per_rep[i] + n_dropout_sum_per_rep[i] + n_dropin_sum_per_rep[i] + n_ignore_sum_per_rep[i]
      id_mk_dedup1[[i]] <- unlist(mapply(rep, rep(1:n_mk, 4), c(n_general_rep1, n_dropout_rep1, n_dropin_rep1, n_ignore_rep1)))
      at_dropout[[i]] <- unlist(mapply(rep, ats_use, n_dropout_rep1))
    }
    id_mk_dedup1 <- unlist(id_mk_dedup1)
    at_dropout <- unlist(at_dropout)

    n_general_sum <- sum(n_general_sum_per_rep)
    n_dropout_sum <- sum(n_dropout_sum_per_rep)
    n_dropin_sum <- sum(n_dropin_sum_per_rep)
    n_ignore_sum <- sum(n_dedup1_sum_per_rep)
    n_dedup1_sum <- sum(n_dedup1_sum_per_rep)

    loglike_di <- unlist(loglike_di_use)
    if(is.null(loglike_di)){
      loglike_di <- numeric(0)
    }

    ##################################################
    # Arrange objects for de-duplication (pattern 1) #
    ##################################################

    index_dedup1_use <- matrix(NA, nrow = 2 * hnc, ncol = n_dedup1_sum)
    size_dedup1_use <- rep(NA, n_dedup1_sum)
    oph_dedup1_use <- rep(NA, n_dedup1_sum)

    index_general_start <- 1
    for(i in 1:n_rep){
      n_general_rep1 <- n_general_use[(n_mk * (i - 1) + 1):(n_mk * i)]
      n_dropout_rep1 <- n_dropout_use[(n_mk * (i - 1) + 1):(n_mk * i)]
      n_dropin_rep1 <- n_dropin_use[(n_mk * (i - 1) + 1):(n_mk * i)]
      n_ignore_rep1 <- n_ignore_use[(n_mk * (i - 1) + 1):(n_mk * i)]

      index_dropout_start <- index_general_start + n_general_sum_per_rep[i]
      index_dropin_start <- index_dropout_start + n_dropout_sum_per_rep[i]
      index_ignore_start <- index_dropin_start + n_dropin_sum_per_rep[i]
      sumup_index_dedup1 <- 0

      for(j in 1:n_mk){
        n_general <- n_general_rep1[j]
        if(n_general != 0){
          index_tmp <- index_general_start:(index_general_start + n_general - 1)
          index_dedup1_use[, index_tmp] <- index_dedup1_general_use[[i]][[j]] + sumup_index_dedup1
          size_dedup1_use[index_tmp] <- size_dedup1_general_use[[i]][[j]]
          oph_dedup1_use[index_tmp] <- oph_dedup1_general_use[[i]][[j]]
        }

        n_dropout <- n_dropout_rep1[j]
        if(n_dropout != 0){
          index_tmp <- index_dropout_start:(index_dropout_start + n_dropout - 1)
          index_dedup1_use[, index_tmp] <- index_dedup1_dropout_use[[i]][[j]] + sumup_index_dedup1
          size_dedup1_use[index_tmp] <- size_dedup1_dropout_use[[i]][[j]]
          oph_dedup1_use[index_tmp] <- oph_dedup1_dropout_use[[i]][[j]]
        }

        n_dropin <- n_dropin_rep1[j]
        if(n_dropin != 0){
          index_tmp <- index_dropin_start:(index_dropin_start + n_dropin - 1)
          index_dedup1_use[, index_tmp] <- index_dedup1_dropin_use[[i]][[j]] + sumup_index_dedup1
          size_dedup1_use[index_tmp] <- size_dedup1_dropin_use[[i]][[j]]
          oph_dedup1_use[index_tmp] <- oph_dedup1_dropin_use[[i]][[j]]
        }

        n_ignore <- n_ignore_rep1[j]
        if(n_ignore != 0){
          index_tmp <- index_ignore_start:(index_ignore_start + n_ignore - 1)
          index_dedup1_use[, index_tmp] <- index_dedup1_ignore_use[[i]][[j]] + sumup_index_dedup1
          size_dedup1_use[index_tmp] <- size_dedup1_ignore_use[[i]][[j]]
          oph_dedup1_use[index_tmp] <- oph_dedup1_ignore_use[[i]][[j]]
        }

        index_general_start <- index_general_start + n_general
        index_dropout_start <- index_dropout_start + n_dropout
        index_dropin_start <- index_dropin_start + n_dropin
        index_ignore_start <- index_ignore_start + n_ignore
        sumup_index_dedup1 <- sumup_index_dedup1 + n_rph_base_use[j]
      }

      index_general_start <- 1 + sum(n_dedup1_sum_per_rep[1:i])
    }

    ##################################################
    # Prepare objects for de-duplication (pattern 2) #
    ##################################################

    n_dedup2_use <- n_peak_use * n_gtc_use
    n_dedup2_sum_per_rep <- rep(sum(n_dedup2_use), n_rep)
    n_dedup2_sum <- sum(n_dedup2_sum_per_rep)
    index_dedup2_use <- rep(NA, n_dedup2_sum)
    index_dedup1_start <- 1
    index_dedup2_start <- 1

    for(i in 1:n_rep){
      index_dedup1_use_rep1 <- index_dedup1_use[, index_dedup1_start:(index_dedup1_start + n_dedup1_sum_per_rep[i] - 1)]
      sumup_index_dedup2 <- 0
      for(j in 1:n_mk){
        index_rph_gtc <- index_rph_gtc_use[[i]][[j]] + sumup_index_dedup2
        n_dedup2 <- n_dedup2_use[j]
        index_dedup2_use[index_dedup2_start:(index_dedup2_start + n_dedup2 - 1)] <- make_index_dedup2_mk(index_rph_gtc, index_dedup1_use_rep1)

        index_dedup2_start <- index_dedup2_start + n_dedup2
        sumup_index_dedup2 <- sumup_index_dedup2 + n_rph_base_use[j]
      }
      index_dedup1_start <- index_dedup1_start + n_dedup1_sum_per_rep[i]
    }

    ################
    # Perform MCMC #
    ################

    # Prepare date list
    data_list_main <- list(hnc = hnc,
                           n_rep = n_rep,
                           n_mk = n_mk,
                           ltph = max_tph / 4,
                           htph = max_tph / 2,
                           n_rph_base_sum = sum(n_rph_base_use),
                           rph_base_use = unlist(rph_base_use),
                           ae_alpha_use = ae_alpha_use,
                           cv = cv,
                           mean_size = mean_size,
                           n_dedup1_sum_per_rep = n_dedup1_sum_per_rep,
                           n_dedup1_sum = n_dedup1_sum,
                           index_dedup1_use = index_dedup1_use,
                           size_dedup1_use = size_dedup1_use,
                           oph_dedup1_use = oph_dedup1_use,
                           id_mk_dedup1 = id_mk_dedup1,
                           n_general_sum_per_rep = n_general_sum_per_rep,
                           n_dropout_sum_per_rep = n_dropout_sum_per_rep,
                           n_dropout_sum = sum(n_dropout_sum_per_rep),
                           at_dropout = at_dropout,
                           n_dropin_sum = n_dropin_sum,
                           n_dropin_sum_per_rep = n_dropin_sum_per_rep,
                           loglike_di = loglike_di,
                           n_dedup2_sum = n_dedup2_sum,
                           n_dedup2_sum_per_rep = n_dedup2_sum_per_rep,
                           n_dedup2_use = n_dedup2_use,
                           index_dedup2_use = index_dedup2_use,
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
    tph <- medians[grep("tph", variables)]
    mp <- medians[grep("mp", variables)]
    deg <- medians[grep("deg", variables)]
    pae <- medians[grep("pae", variables)]
    ae <- n_mk * pae
    rhat_lp <- rhats[which(variables == "lp__")]
    rhat_tph <- rhats[grep("tph", variables)]
    rhat_mp <- rhats[grep("mp", variables)]
    rhat_deg <- rhats[grep("deg", variables)]
    rhat_pae <- rhats[grep("pae", variables)]

    ######################################################
    # Calculate likelihood under the estimated parameter #
    ######################################################

    index_dedup1_start <- 1
    index_dedup2_start <- 1
    index_dropout_start <- 1
    index_dropin_start <- 1
    lpd_gtc_use <- vector(mode = "list", length = n_rep)
    for(i in 1:n_rep){
      n_dedup1_sum_rep1 <- n_dedup1_sum_per_rep[i]
      index_dedup1_use_rep1 <- index_dedup1_use[, index_dedup1_start:(index_dedup1_start + n_dedup1_sum_rep1 - 1)]
      index_dedup1_use_cpp <- asplit(index_dedup1_use_rep1 - 1, 1) # C++ index

      size_dedup1_use_rep1 <- size_dedup1_use[index_dedup1_start:(index_dedup1_start + n_dedup1_sum_rep1 - 1)]
      oph_dedup1_use_rep1 <- oph_dedup1_use[index_dedup1_start:(index_dedup1_start + n_dedup1_sum_rep1 - 1)]

      id_mk_dedup1_rep1 <- id_mk_dedup1[index_dedup1_start:(index_dedup1_start + n_dedup1_sum_rep1 - 1)]
      ae_dedup1_use_rep1 <- ae[id_mk_dedup1_rep1]

      n_dropout_sum_rep1 <- n_dropout_sum_per_rep[i]
      if(n_dropout_sum_rep1 != 0){
        at_dropout_rep1 <- at_dropout[index_dropout_start:(index_dropout_start + n_dropout_sum_rep1 - 1)]
      }else{
        at_dropout_rep1 <- numeric(0)
      }

      n_dropin_sum_rep1 <- n_dropin_sum_per_rep[i]
      if(n_dropin_sum_rep1 != 0){
        loglike_di_rep1 <- loglike_di[index_dropin_start:(index_dropin_start + n_dropin_sum_rep1 - 1)]
      }else{
        loglike_di_rep1 <- numeric(0)
      }

      n_dedup2_sum_rep1 <- n_dedup2_sum_per_rep[i]
      index_dedup2_use_rep1 <- index_dedup2_use[index_dedup2_start:(index_dedup2_start + n_dedup2_sum_rep1 - 1)]
      index_dedup2_use_cpp <- index_dedup2_use_rep1 - 1 # C++ index

      lpd_gtc_use[[i]] <- calc_lpd_gtc_all_mk_cpp(n_mk, mp, deg, hnc, n_dedup1_sum_rep1, n_general_sum_per_rep[i], n_dropout_sum_per_rep[i], n_dropin_sum_per_rep[i],
                                                  unlist(rph_base_use), index_dedup1_use_cpp, size_dedup1_use_rep1, mean_size, ae_dedup1_use_rep1, oph_dedup1_use_rep1,
                                                  at_dropout_rep1, loglike_di_rep1, cv, tph[i],
                                                  index_dedup2_use_cpp, n_dedup2_use, n_peak_use, n_gtc_use)

      index_dedup1_start <- index_dedup1_start + n_dedup1_sum_rep1
      index_dropout_start <- index_dropout_start + n_dropout_sum_rep1
      index_dropin_start <- index_dropin_start + n_dropin_sum_rep1
      index_dedup2_start <- index_dedup2_start + n_dedup2_sum_rep1
    }

    dt_pg_list <- make_dt_pg_list(markers_use, paste0("U", 1:hnc), gtc_use, gtc_logprob_use, lpd_gtc_use)

    log10like_use <- calc_log10like_use(dt_pg_list)

    dt_log10like <- data.table(Marker = c(markers_use, "Overall"),
                               Log10_Likelihood = log10like_use)

    dt_diag <- data.table("Parameter" = c("lp__", rep("Total peak height", n_rep), rep("Mixture proportion", hnc), rep("Degradation", hnc), rep("Amplification efficiency", n_mk)),
                          "Contributor" = c(NA, rep(NA, n_rep), paste0("U", 1:hnc), paste0("U", 1:hnc), rep(NA, n_mk)),
                          "Replicate" = c(NA, 1:n_rep, rep(NA, hnc), rep(NA, hnc), rep(NA, n_mk)),
                          "Marker" = c(NA, rep(NA, n_rep), rep(NA, hnc), rep(NA, hnc), markers_use),
                          "Estimated value" = c(lp, tph, mp[hnc:1], deg[hnc:1], ae),
                          "Rhat" = c(rhat_lp, rhat_tph, rhat_mp[hnc:1], rhat_deg[hnc:1], rhat_pae))

    # End computation
    end_time <- proc.time()
    comp_time <- end_time[3] - start_time[3]

    return(list("dt_pg_list" = dt_pg_list,
                "dt_log10like" = dt_log10like,
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

    return(list("dt_pg_list" = NULL,
                "dt_log10like" = NULL,
                "dt_diag" = NULL,
                "traceplot" = NULL,
                "comp_time" = NULL,
                "error_message" = error_message))
  }
}

#' perform_deconvo_all_noc
#'
#' @description The function to perform mixture deconvolution for all assumed number of contributors
#' @param dt_csp The data.table of the crime stain profiles
#' @param dt_af The data.table of the allele frequencies
#' @param dt_kit_auto The data.table of the information on the autosomal markers in a kit
#' @param dt_ignore_mk The data.table of the information on the ignored markers
#' @param dt_anameth The data.table of the analysis method
#' @param alpha_auto A vector of the estimated alpha values for Dirichlet distribution of the locus-specific amplification efficiencies
#' @param cv A parameter of the gamma distribution
#' @param list_dt_average_sr A list of the data.table for the average stutter ratios
#' @param hnc_all A vector of the assumed number of contributors
perform_deconvo_all_noc <- function(dt_csp, dt_af, dt_kit_auto, dt_ignore_mk, dt_anameth, alpha_auto, cv, list_dt_average_sr, hnc_all){
  markers_use <- dt_ignore_mk$Marker[!dt_ignore_mk$Ignore]
  n_mk <- length(markers_use)
  if(n_mk != 0){
    n_rep <- length(unique(dt_csp$SampleName))

    # Extract ATs of the used markers
    markers_at <- dt_anameth$Subcategory[dt_anameth$Category == "AT"]
    ats_auto <- dt_anameth$Value[dt_anameth$Category == "AT"]
    ats_use <- as.numeric(ats_auto[match(markers_use, markers_at)])
    rm(markers_at, ats_auto)

    # Extract peak data from crime stain profiles
    tmp <- make_peak_data(dt_csp, dt_kit_auto, markers_use, ats_use)
    peaks_list <- tmp$peaks_list
    sizes_list <- tmp$sizes_list
    heights_list <- tmp$heights_list
    mean_size <- mean(dt_kit_auto$Size)

    # Set allele probabilities
    tmp <- set_al_prob(markers_use, dt_af, peaks_list)
    pop_al_list <- tmp$pop_al_list
    pop_prob_list <- tmp$pop_prob_list

    # Extract the theta value
    theta <- as.numeric(dt_anameth$Value[dt_anameth$Subcategory == "Theta"])

    # Extract filters to remove unrealistic genotype combinations
    bsr_fltr <- dt_anameth$Value[dt_anameth$Subcategory == "Maximum BSR"]
    fsr_fltr <- dt_anameth$Value[dt_anameth$Subcategory == "Maximum FSR"]
    dsr_fltr <- dt_anameth$Value[dt_anameth$Subcategory == "Maximum DSR"]
    m2sr_fltr <- dt_anameth$Value[dt_anameth$Subcategory == "Maximum M2SR"]
    n_dropin_fltr <- dt_anameth$Value[dt_anameth$Subcategory == "Maximum number of drop-in"]
    h_dropin_fltr <- dt_anameth$Value[dt_anameth$Subcategory == "Maximum height of drop-in"]

    # Extract drop-in parameters
    logprob_di <- log(dt_anameth$Value[dt_anameth$Subcategory == "Probability of drop-in"])
    shape_di <- dt_anameth$Value[dt_anameth$Subcategory == "Shape parameter"]
    rate_di <- dt_anameth$Value[dt_anameth$Subcategory == "Rate parameter"]

    # Extract MCMC parameters
    warmup <- dt_anameth$Value[dt_anameth$Subcategory == "Warmup"]
    sampling <- dt_anameth$Value[dt_anameth$Subcategory == "Sampling"]
    chains <- dt_anameth$Value[dt_anameth$Subcategory == "Chains"]

    # Set LSAE
    ae_alpha_use <- alpha_auto[match(markers_use, names(alpha_auto))]

    n_hnc <- length(hnc_all)
    data_hnc_list <- vector(mode = "list", length = n_hnc)
    for(i in 1:n_hnc){
      hnc <- hnc_all[i]
      data_hnc_list[[i]] <- perform_deconvo(hnc, n_rep, markers_use, ats_use, peaks_list, sizes_list, heights_list, mean_size,
                                            pop_al_list, pop_prob_list, theta,
                                            bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr, n_dropin_fltr, h_dropin_fltr,
                                            logprob_di, shape_di, rate_di,
                                            warmup, sampling, chains,
                                            ae_alpha_use, cv, list_dt_average_sr)
    }
  }else{
    peaks_list <- NULL
    heights_list <- NULL
    pop_al_list <- NULL
    pop_prob_list <- NULL
    mean_size <- NULL
    data_hnc_list <- NULL
  }

  return(list("hnc_all" = hnc_all,
              "markers_use" = markers_use,
              "peaks_list" = peaks_list,
              "heights_list" = heights_list,
              "pop_al_list" = pop_al_list,
              "pop_prob_list" = pop_prob_list,
              "mean_size" = mean_size,
              "data_hnc_list" = data_hnc_list))
}
