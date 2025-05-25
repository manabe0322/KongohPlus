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

#' generate_init
#'
#' @description The function to generate initial values
#' @param n_rep The number of replicates
#' @param hnc The hypothesized number of contributors
#' @param tph_mean The mean of the total peak height
#' @param tph_sd The SD of the total peak height
generate_init <- function(n_rep, hnc, tph_mean, tph_sd){
  tph <- rnorm(n_rep, tph_mean, tph_sd)
  lambda <- rnorm(hnc, 0, 1)
  deg <- rbeta(hnc, 1, 1)
  return(list(list(tph = tph, lambda = lambda, deg = deg)))
}

#' perform_mle
#'
#' @description The function to perform MLE
#' @param mdoel_main The stan model for deconvolution
#' @param data_list_main The input data for stan
#' @param n_rep The number of replicates
#' @param hnc The hypothesized number of contributors
#' @param tph_mean The mean of the total peak height
#' @param tph_sd The SD of the total peak height
#' @param n_conv The number of successful convergences
#' @param n_max_init The maximum number of initializations
perform_mle <- function(model_main, data_list_main, n_rep, hnc, tph_mean, tph_sd, n_conv, n_max_init){
  count_all <- 0
  count_conv <- 0
  max_lp <- -Inf

  if(hnc == 1){
    est_param <- rep(NA, n_rep + 1)
  }else{
    est_param <- rep(NA, n_rep + 3 * hnc)
  }

  # Repeat random initialization and optimization until the condition is satisfied
  repeat{
    count_all <- count_all + 1
    init_vals <- generate_init(n_rep, hnc, tph_mean, tph_sd)
    fit <- model_main$optimize(data_list_main, init = init_vals)

    # If convergence is successful, perform the following steps
    if(fit$return_codes() == 0){
      lp <- fit$lp()

      # If the difference between the current likelihood and the maximum likelihood is less than 1e-6, the current run is considered convergent.
      if(abs(max_lp - lp) < 1e-6){
        count_conv <- count_conv + 1
        if(lp - max_lp > 0){
          max_lp <- lp
          est_param <- fit$mle()
        }
      }

      # Reset convergence count to 1 if a new maximum likelihood is found.
      else if(lp - max_lp > 0){
        count_conv <- 1
        max_lp <- lp
        est_param <- fit$mle()
      }

      # Stop the simulation if convergence is successful for the required number of times.
      if(count_conv == n_conv){
        break
      }

      # Stop the simulation if the maximum number of simulations is reached, treating it as a convergence failure.
      if(count_all == n_max_init){
        if(hnc == 1){
          est_param <- rep(NA, n_rep + 1)
        }else{
          est_param <- rep(NA, n_rep + 3 * hnc)
        }
        break
      }
    }
  }

  if(hnc == 1){
    names(est_param) <- c(paste0("tph[", 1:n_rep, "]"), "deg[1]")
  }else{
    names(est_param) <- c(paste0("tph[", 1:n_rep, "]"), paste0("lambda[", 1:hnc, "]"), paste0("deg[", 1:hnc, "]"), paste0("mp[", 1:hnc, "]"))
  }
  return(est_param)
}

#' perform_deconvo
#'
#' @description The function to perform mixture deconvolution
#' @param hnc The number of contributors
#' @param cv_allele The coefficient of variation for allelels
#' @param cv_stutter The coefficient of variation for stutters
#' @param list_dt_average_sr A list of the data.table for the average stutter ratios
#' @param n_conv The number of successful convergences
#' @param n_max_init The maximum number of initializations
perform_deconvo <- function(hnc, n_rep, markers_use, ats_use, peaks_list, sizes_list, heights_list, min_size,
                            pop_al_list, pop_prob_list, theta,
                            bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr, n_dropin_fltr, h_dropin_fltr,
                            freq_di, shape_di, rate_di,
                            ae_use, cv_allele, cv_stutter, list_dt_average_sr,
                            n_conv, n_max_init){
  # Start computation
  start_time <- proc.time()

  # Define objects
  n_mk <- length(markers_use)
  n_peak_use <- rep(NA, n_mk)
  n_gtc_use <- rep(NA, n_mk)
  gtc_use <- vector(mode = "list", length = n_mk)
  names(gtc_use) <- markers_use
  median_sizes_use <- rep(NA, n_mk)
  tph_per_mk_use <- matrix(0, n_mk, n_rep)
  rph_base_use <- vector(mode = "list", length = n_mk)
  gtc_logprob_use <- vector(mode = "list", length = n_mk)
  index_rph_gtc_use <- vector(mode = "list", length = n_rep)
  index_dedup1_general_use <- vector(mode = "list", length = n_rep)
  index_dedup1_dropout_use <- vector(mode = "list", length = n_rep)
  index_dedup1_dropin_use <- vector(mode = "list", length = n_rep)
  index_dedup1_ignore_use <- vector(mode = "list", length = n_rep)
  cv_square_dedup1_general_use <- vector(mode = "list", length = n_rep)
  cv_square_dedup1_dropout_use <- vector(mode = "list", length = n_rep)
  size_dedup1_general_use <- vector(mode = "list", length = n_rep)
  size_dedup1_dropout_use <- vector(mode = "list", length = n_rep)
  oph_dedup1_general_use <- vector(mode = "list", length = n_rep)
  oph_dedup1_dropout_use <- vector(mode = "list", length = n_rep)
  loglike_di_use <- vector(mode = "list", length = n_rep)
  for(i in 1:n_rep){
    index_rph_gtc_use[[i]] <- vector(mode = "list", length = n_mk)
    index_dedup1_general_use[[i]] <- vector(mode = "list", length = n_mk)
    index_dedup1_dropout_use[[i]] <- vector(mode = "list", length = n_mk)
    index_dedup1_dropin_use[[i]] <- vector(mode = "list", length = n_mk)
    index_dedup1_ignore_use[[i]] <- vector(mode = "list", length = n_mk)
    cv_square_dedup1_general_use[[i]] <- vector(mode = "list", length = n_mk)
    cv_square_dedup1_dropout_use[[i]] <- vector(mode = "list", length = n_mk)
    size_dedup1_general_use[[i]] <- vector(mode = "list", length = n_mk)
    size_dedup1_dropout_use[[i]] <- vector(mode = "list", length = n_mk)
    oph_dedup1_general_use[[i]] <- vector(mode = "list", length = n_mk)
    oph_dedup1_dropout_use[[i]] <- vector(mode = "list", length = n_mk)
    loglike_di_use[[i]] <- vector(mode = "list", length = n_mk)
  }

  for(i in 1:n_mk){
    # Extract peak data of one marker
    mk <- markers_use[i]
    peaks <- peaks_list[[mk]]
    sizes <- sizes_list[[mk]]
    heights_all_rep <- heights_list[[mk]]
    n_peak_use[i] <- n_peak <- length(peaks)
    median_sizes_use[i] <- median(sizes)
    tph_per_mk_use[i, ] <- sapply(heights_all_rep, sum) / 2

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

    # Calculate the probability of drop-in in the current marker
    auc_below_at <- pgamma(q = at, shape = shape_di, rate = rate_di)
    logprob_di <- log(freq_di / (1 - auc_below_at))

    # The base of the relative peak heights
    rph_base_use[[i]] <- rph_base_mk <- make_rph_base(peaks, bsrs, fsrs, dsrs, m2srs, ae_use[mk])

    cv_square_base_mk <- c(rep(cv_allele^2, n_peak), rep(cv_stutter^2, 4 * n_peak), 0, -1)

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

        # CV for de-duplicated objects (pattern 1)
        cv_square_dedup1_tmp <- matrix(cv_square_base_mk[index_dedup1_tmp], nrow = 2 * hnc)

        # Sizes for de-duplicated objects (pattern 1)
        size_dedup1_tmp <- rep(sizes, n_gtc)
        size_dedup1_tmp <- size_dedup1_tmp[bool_duplicated_adopt]

        # Observe peak heights for de-duplicated objects (pattern 1)
        heights <- heights_all_rep[[j]]
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

        cv_square_dedup1_general_use[[j]][[i]] <- cv_square_dedup1_tmp[, index_general, drop = FALSE]
        cv_square_dedup1_dropout_use[[j]][[i]] <- cv_square_dedup1_tmp[, index_dropout, drop = FALSE]

        size_dedup1_general_use[[j]][[i]] <- size_dedup1_tmp[index_general]
        size_dedup1_dropout_use[[j]][[i]] <- size_dedup1_tmp[index_dropout]

        oph_dedup1_general_use[[j]][[i]] <- oph_dedup1_tmp[index_general]
        oph_dedup1_dropout_use[[j]][[i]] <- oph_dedup1_tmp[index_dropout]

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
    n_dedup1_sum_per_rep <- rep(NA, n_rep)
    for(i in 1:n_rep){
      n_general_use[(n_mk * (i - 1) + 1):(n_mk * i)] <- n_general_rep1 <- sapply(index_dedup1_general_use[[i]], ncol)
      n_dropout_use[(n_mk * (i - 1) + 1):(n_mk * i)] <- n_dropout_rep1 <- sapply(index_dedup1_dropout_use[[i]], ncol)
      n_dropin_use[(n_mk * (i - 1) + 1):(n_mk * i)] <- n_dropin_rep1 <- sapply(index_dedup1_dropin_use[[i]], ncol)
      n_ignore_use[(n_mk * (i - 1) + 1):(n_mk * i)] <- n_ignore_rep1 <- sapply(index_dedup1_ignore_use[[i]], ncol)
      n_general_sum_per_rep[i] <- sum(n_general_rep1)
      n_dropout_sum_per_rep[i] <- sum(n_dropout_rep1)
      n_dropin_sum_per_rep[i] <- sum(n_dropin_rep1)
      n_dedup1_sum_per_rep[i] <- n_general_sum_per_rep[i] + n_dropout_sum_per_rep[i] + n_dropin_sum_per_rep[i] + sum(n_ignore_rep1)
    }

    n_general_sum <- sum(n_general_sum_per_rep)
    n_dropout_sum <- sum(n_dropout_sum_per_rep)
    n_dropin_sum <- sum(n_dropin_sum_per_rep)
    n_dedup1_sum <- sum(n_dedup1_sum_per_rep)
    n_rph_sum_per_rep <- n_general_sum_per_rep + n_dropout_sum_per_rep
    n_rph_sum <- n_general_sum + n_dropout_sum

    loglike_di <- unlist(loglike_di_use)
    if(is.null(loglike_di)){
      loglike_di <- numeric(0)
    }

    ##################################################
    # Arrange objects for de-duplication (pattern 1) #
    ##################################################

    index_dedup1_use <- matrix(NA, nrow = 2 * hnc, ncol = n_dedup1_sum)
    index_dedup1_mcmc <- matrix(NA, nrow = 2 * hnc, ncol = n_rph_sum)
    cv_square_mcmc <- matrix(NA, nrow = 2 * hnc, ncol = n_rph_sum)
    size_mcmc <- rep(NA, n_rph_sum)
    oph_mcmc <- rep(NA, n_rph_sum)
    at_dropout <- rep(NA, n_dropout_sum)

    index_start_rph_general <- 1
    index_start_dedup1_use_general <- 1
    index_at_start <- 1
    for(i in 1:n_rep){
      n_general_rep1 <- n_general_use[(n_mk * (i - 1) + 1):(n_mk * i)]
      n_dropout_rep1 <- n_dropout_use[(n_mk * (i - 1) + 1):(n_mk * i)]
      n_dropin_rep1 <- n_dropin_use[(n_mk * (i - 1) + 1):(n_mk * i)]
      n_ignore_rep1 <- n_ignore_use[(n_mk * (i - 1) + 1):(n_mk * i)]

      index_start_rph_dropout <- index_start_rph_general + n_general_sum_per_rep[i]

      index_start_dedup1_use_dropout <- index_start_dedup1_use_general + n_general_sum_per_rep[i]
      index_start_dedup1_use_dropin <- index_start_dedup1_use_dropout + n_dropout_sum_per_rep[i]
      index_start_dedup1_use_ignore <- index_start_dedup1_use_dropin + n_dropin_sum_per_rep[i]
      sumup_index_dedup1 <- 0

      for(j in 1:n_mk){
        n_general <- n_general_rep1[j]
        if(n_general != 0){
          index_tmp <- index_start_rph_general:(index_start_rph_general + n_general - 1)
          index_dedup1_mcmc[, index_tmp] <- index_dedup1_general_use[[i]][[j]] + sumup_index_dedup1
          cv_square_mcmc[, index_tmp] <- cv_square_dedup1_general_use[[i]][[j]]
          size_mcmc[index_tmp] <- size_dedup1_general_use[[i]][[j]]
          oph_mcmc[index_tmp] <- oph_dedup1_general_use[[i]][[j]]

          index_tmp <- index_start_dedup1_use_general:(index_start_dedup1_use_general + n_general - 1)
          index_dedup1_use[, index_tmp] <- index_dedup1_general_use[[i]][[j]] + sumup_index_dedup1
        }

        n_dropout <- n_dropout_rep1[j]
        if(n_dropout != 0){
          index_tmp <- index_start_rph_dropout:(index_start_rph_dropout + n_dropout - 1)
          index_dedup1_mcmc[, index_tmp] <- index_dedup1_dropout_use[[i]][[j]] + sumup_index_dedup1
          cv_square_mcmc[, index_tmp] <- cv_square_dedup1_dropout_use[[i]][[j]]
          size_mcmc[index_tmp] <- size_dedup1_dropout_use[[i]][[j]]
          oph_mcmc[index_tmp] <- oph_dedup1_dropout_use[[i]][[j]]

          at_dropout[index_at_start:(index_at_start + n_dropout - 1)] <- ats_use[j]

          index_tmp <- index_start_dedup1_use_dropout:(index_start_dedup1_use_dropout + n_dropout - 1)
          index_dedup1_use[, index_tmp] <- index_dedup1_dropout_use[[i]][[j]] + sumup_index_dedup1
        }

        n_dropin <- n_dropin_rep1[j]
        if(n_dropin != 0){
          index_tmp <- index_start_dedup1_use_dropin:(index_start_dedup1_use_dropin + n_dropin - 1)
          index_dedup1_use[, index_tmp] <- index_dedup1_dropin_use[[i]][[j]] + sumup_index_dedup1
        }

        n_ignore <- n_ignore_rep1[j]
        if(n_ignore != 0){
          index_tmp <- index_start_dedup1_use_ignore:(index_start_dedup1_use_ignore + n_ignore - 1)
          index_dedup1_use[, index_tmp] <- index_dedup1_ignore_use[[i]][[j]] + sumup_index_dedup1
        }

        index_start_rph_general <- index_start_rph_general + n_general
        index_start_rph_dropout <- index_start_rph_dropout + n_dropout

        index_start_dedup1_use_general <- index_start_dedup1_use_general + n_general
        index_start_dedup1_use_dropout <- index_start_dedup1_use_dropout + n_dropout
        index_start_dedup1_use_dropin <- index_start_dedup1_use_dropin + n_dropin
        index_start_dedup1_use_ignore <- index_start_dedup1_use_ignore + n_ignore
        index_at_start <- index_at_start + n_dropout
        sumup_index_dedup1 <- sumup_index_dedup1 + n_rph_base_use[j]
      }

      index_start_rph_general <- n_rph_sum_per_rep[i] + 1
      index_start_dedup1_use_general <- n_dedup1_sum_per_rep[i] + 1
    }

    deg_size_mcmc <- (size_mcmc - min_size) / min_size

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

    ##########################################################
    # Determine the prior distribution for the parameter tph #
    ##########################################################

    median_sizes_use <- rep(median_sizes_use, n_rep)
    tph_per_mk_use <- as.vector(tph_per_mk_use)
    fit_exp <- lm(log(tph_per_mk_use) ~ median_sizes_use)
    tph_mean <- exp(coef(fit_exp)[1]) * exp(coef(fit_exp)[2] * min_size)

    min_tph <- min(tph_per_mk_use)
    max_tph <- max(tph_per_mk_use)
    if(min_tph == max_tph){
      min_tph <- 0
    }
    tph_sd <- (max_tph - min_tph) / 2

    ###############
    # Perform MLE #
    ###############

    if(n_rep == 1){
      if(hnc == 1){
        model_main <- stan_package_model(name = "model_main_rep1_single", package = "KongohPlus")

        data_list_main <- list(n_mk = n_mk,
                               tph_mean = tph_mean,
                               tph_sd = tph_sd,
                               n_rph_base_sum = sum(n_rph_base_use),
                               rph_base_use = unlist(rph_base_use),
                               n_dedup1_sum = n_dedup1_sum,
                               n_rph_sum = n_rph_sum,
                               index_dedup1_rep1 = index_dedup1_mcmc,
                               cv_square_rep1 = t(cv_square_mcmc),
                               deg_size_rep1 = deg_size_mcmc,
                               oph_rep1 = oph_mcmc,
                               n_general_sum = n_general_sum,
                               n_dropout_sum = n_dropout_sum,
                               at_dropout_rep1 = at_dropout,
                               n_dropin_sum = n_dropin_sum,
                               loglike_di_rep1 = loglike_di,
                               n_dedup2_sum = n_dedup2_sum,
                               n_dedup2_use = n_dedup2_use,
                               index_dedup2_use_rep1 = index_dedup2_use,
                               n_peak_use = n_peak_use,
                               n_gtc_use = n_gtc_use,
                               gtc_logprob_use = unlist(gtc_logprob_use))
      }else{
        model_main <- stan_package_model(name = "model_main_rep1", package = "KongohPlus")

        data_list_main <- list(hnc = hnc,
                               n_mk = n_mk,
                               tph_mean = tph_mean,
                               tph_sd = tph_sd,
                               n_rph_base_sum = sum(n_rph_base_use),
                               rph_base_use = unlist(rph_base_use),
                               n_dedup1_sum = n_dedup1_sum,
                               n_rph_sum = n_rph_sum,
                               index_dedup1_rep1 = index_dedup1_mcmc,
                               cv_square_rep1 = t(cv_square_mcmc),
                               deg_size_rep1 = deg_size_mcmc,
                               oph_rep1 = oph_mcmc,
                               n_general_sum = n_general_sum,
                               n_dropout_sum = n_dropout_sum,
                               at_dropout_rep1 = at_dropout,
                               n_dropin_sum = n_dropin_sum,
                               loglike_di_rep1 = loglike_di,
                               n_dedup2_sum = n_dedup2_sum,
                               n_dedup2_use = n_dedup2_use,
                               index_dedup2_use_rep1 = index_dedup2_use,
                               n_peak_use = n_peak_use,
                               n_gtc_use = n_gtc_use,
                               gtc_logprob_use = unlist(gtc_logprob_use))
      }
    }else{
      if(hnc == 1){
        model_main <- stan_package_model(name = "model_main_single", package = "KongohPlus")

        data_list_main <- list(n_rep = n_rep,
                               n_mk = n_mk,
                               tph_mean = tph_mean,
                               tph_sd = tph_sd,
                               n_rph_base_sum = sum(n_rph_base_use),
                               rph_base_use = unlist(rph_base_use),
                               n_dedup1_sum_per_rep = n_dedup1_sum_per_rep,
                               n_rph_sum_per_rep = n_rph_sum_per_rep,
                               n_rph_sum = n_rph_sum,
                               index_dedup1_mcmc = index_dedup1_mcmc,
                               cv_square_mcmc = t(cv_square_mcmc),
                               deg_size_mcmc = deg_size_mcmc,
                               oph_mcmc = oph_mcmc,
                               n_general_sum_per_rep = n_general_sum_per_rep,
                               n_dropout_sum_per_rep = n_dropout_sum_per_rep,
                               n_dropout_sum = n_dropout_sum,
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
      }else{
        model_main <- stan_package_model(name = "model_main", package = "KongohPlus")

        data_list_main <- list(hnc = hnc,
                               n_rep = n_rep,
                               n_mk = n_mk,
                               tph_mean = tph_mean,
                               tph_sd = tph_sd,
                               n_rph_base_sum = sum(n_rph_base_use),
                               rph_base_use = unlist(rph_base_use),
                               n_dedup1_sum_per_rep = n_dedup1_sum_per_rep,
                               n_rph_sum_per_rep = n_rph_sum_per_rep,
                               n_rph_sum = n_rph_sum,
                               index_dedup1_mcmc = index_dedup1_mcmc,
                               cv_square_mcmc = t(cv_square_mcmc),
                               deg_size_mcmc = deg_size_mcmc,
                               oph_mcmc = oph_mcmc,
                               n_general_sum_per_rep = n_general_sum_per_rep,
                               n_dropout_sum_per_rep = n_dropout_sum_per_rep,
                               n_dropout_sum = n_dropout_sum,
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
      }
    }

    # Perform MLE
    est_param <- perform_mle(model_main, data_list_main, n_rep, hnc, tph_mean, tph_sd, n_conv, n_max_init)

    if(!any(is.na(est_param))){
      tph <- est_param[paste0("tph[", 1:n_rep, "]")]
      deg <- est_param[paste0("deg[", 1:hnc, "]")]
      if(hnc == 1){
        mp <- 1
      }else{
        mp <- est_param[paste0("mp[", 1:hnc, "]")]
      }

      ######################################################
      # Calculate likelihood under the estimated parameter #
      ######################################################

      index_rph_start <- 1
      index_dropout_start <- 1
      index_dropin_start <- 1
      index_dedup2_start <- 1
      lpd_gtc_use <- vector(mode = "list", length = n_rep)
      for(i in 1:n_rep){
        n_rph_sum_rep1 <- n_rph_sum_per_rep[i]
        index_dedup1_rep1 <- index_dedup1_mcmc[, index_rph_start:(index_rph_start + n_rph_sum_rep1 - 1)]
        index_dedup1_rep1_cpp <- asplit(index_dedup1_rep1 - 1, 1) # C++ index
        cv_square_rep1 <- cv_square_mcmc[, index_rph_start:(index_rph_start + n_rph_sum_rep1 - 1)]
        cv_square_rep1 <- asplit(cv_square_rep1, 1)
        size_rep1 <- size_mcmc[index_rph_start:(index_rph_start + n_rph_sum_rep1 - 1)]
        oph_rep1 <- oph_mcmc[index_rph_start:(index_rph_start + n_rph_sum_rep1 - 1)]

        n_dedup1_sum_rep1 <- n_dedup1_sum_per_rep[i]
        n_general_sum_rep1 <- n_general_sum_per_rep[i]

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
        index_dedup2_use_rep1_cpp <- index_dedup2_use_rep1 - 1 # C++ index

        lpd_gtc_use[[i]] <- calc_lpd_gtc_all_mk_cpp(mp, deg, hnc, n_rph_sum_rep1, unlist(rph_base_use), index_dedup1_rep1_cpp, size_rep1, min_size,
                                                    cv_square_rep1, tph[i],
                                                    n_dedup1_sum_rep1, n_general_sum_rep1, n_dropout_sum_rep1, n_dropin_sum_rep1, oph_rep1, at_dropout_rep1, loglike_di_rep1,
                                                    n_mk, index_dedup2_use_rep1_cpp, n_dedup2_use, n_peak_use, n_gtc_use)

        index_rph_start <- index_rph_start + n_rph_sum_rep1
        index_dropout_start <- index_dropout_start + n_dropout_sum_rep1
        index_dropin_start <- index_dropin_start + n_dropin_sum_rep1
        index_dedup2_start <- index_dedup2_start + n_dedup2_sum_rep1
      }

      dt_pg_list <- make_dt_pg_list(markers_use, paste0("U", 1:hnc), gtc_use, gtc_logprob_use, lpd_gtc_use)

      log10like_use <- calc_log10like_use(dt_pg_list)

      dt_log10like <- data.table(Marker = c(markers_use, "Overall"),
                                 Log10_Likelihood = c(log10like_use, sum(log10like_use)))

      dt_param <- data.table("Parameter" = c(rep("Template amount", n_rep), rep("Mixture proportion", hnc), rep("Degradation", hnc)),
                             "Contributor" = c(rep(NA, n_rep), paste0("U", 1:hnc), paste0("U", 1:hnc)),
                             "Replicate" = c(1:n_rep, rep(NA, hnc), rep(NA, hnc)),
                             "Estimated value" = c(tph, mp[hnc:1], deg[hnc:1]))

      # End computation
      end_time <- proc.time()
      comp_time <- end_time[3] - start_time[3]

      return(list("dt_pg_list" = dt_pg_list,
                  "dt_log10like" = dt_log10like,
                  "dt_param" = dt_param,
                  "comp_time" = comp_time,
                  "error_message" = NULL))
    }else{
      error_message <- "Convergence failure."
      return(list("dt_pg_list" = NULL,
                  "dt_log10like" = NULL,
                  "dt_param" = NULL,
                  "comp_time" = NULL,
                  "error_message" = error_message))
    }
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
                "dt_param" = NULL,
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
#' @param expect_ae A vector of the expected locus-specific amplification efficiencies
#' @param cv_allele The coefficient of variation for allelels
#' @param cv_stutter The coefficient of variation for stutters
#' @param list_dt_average_sr A list of the data.table for the average stutter ratios
#' @param hnc_all A vector of the assumed number of contributors
#' @param n_conv The number of successful convergences
#' @param n_max_init The maximum number of initializations
perform_deconvo_all_noc <- function(dt_csp, dt_af, dt_kit_auto, dt_ignore_mk, dt_anameth, expect_ae, cv_allele, cv_stutter, list_dt_average_sr, hnc_all, n_conv, n_max_init){
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
    min_size <- min(dt_kit_auto$Size)

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
    freq_di <- dt_anameth$Value[dt_anameth$Subcategory == "Frequency of drop-in"]
    shape_di <- dt_anameth$Value[dt_anameth$Subcategory == "Shape parameter"]
    rate_di <- dt_anameth$Value[dt_anameth$Subcategory == "Rate parameter"]

    # Extract LSAE
    ae_use <- expect_ae[match(markers_use, names(expect_ae))]

    n_hnc <- length(hnc_all)
    data_hnc_list <- vector(mode = "list", length = n_hnc)
    for(i in 1:n_hnc){
      hnc <- hnc_all[i]
      data_hnc <- try(perform_deconvo(hnc, n_rep, markers_use, ats_use, peaks_list, sizes_list, heights_list, min_size,
                                      pop_al_list, pop_prob_list, theta,
                                      bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr, n_dropin_fltr, h_dropin_fltr,
                                      freq_di, shape_di, rate_di,
                                      ae_use, cv_allele, cv_stutter, list_dt_average_sr,
                                      n_conv, n_max_init))
      if(class(data_hnc) == "try-error"){
        data_hnc_list[[i]] <- list("dt_pg_list" = NULL,
                                   "dt_log10like" = NULL,
                                   "dt_diag" = NULL,
                                   "traceplot" = NULL,
                                   "comp_time" = NULL,
                                   "error_message" = "Unexpected error")
      }else{
        data_hnc_list[[i]] <- data_hnc
      }
    }
  }else{
    peaks_list <- NULL
    heights_list <- NULL
    pop_al_list <- NULL
    pop_prob_list <- NULL
    min_size <- NULL
    data_hnc_list <- NULL
  }

  return(list("hnc_all" = hnc_all,
              "markers_use" = markers_use,
              "peaks_list" = peaks_list,
              "heights_list" = heights_list,
              "pop_al_list" = pop_al_list,
              "pop_prob_list" = pop_prob_list,
              "ae_use" = ae_use,
              "min_size" = min_size,
              "data_hnc_list" = data_hnc_list))
}
