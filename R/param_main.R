#' set_consider_sr_model
#'
#' @description The function to set consideration of each stutter ratio model
#' @param model_sr The stutter model
#' @param reps_one_mk Repeat motifs of one marker
#' @param motif_lus The minimum number of repeat motifs for considering the LUS model
#' @param motif_multi The minimum number of repeat motifs for considering the multi-sequence model
set_consider_sr_model <- function(model_sr, reps_one_mk, motif_lus, motif_multi){
  consider_sr_model <- rep(FALSE, 3)
  names(consider_sr_model) <- c("Allele", "LUS", "Multi-sequence")
  if(model_sr == "Allele"){
    consider_sr_model["Allele"] <- TRUE
  }else if(model_sr == "LUS"){
    consider_sr_model["LUS"] <- TRUE
  }else if(model_sr == "Multi-sequence"){
    consider_sr_model["Multi-sequence"] <- TRUE
  }else{
    consider_sr_model["Allele"] <- TRUE
    max_motif <- max(sapply(reps_one_mk, length))
    if(max_motif >= motif_lus){
      consider_sr_model["LUS"] <- TRUE
    }
    if(max_motif >= motif_multi){
      consider_sr_model["Multi-sequence"] <- TRUE
    }
  }
  return(consider_sr_model)
}

#' calc_aic
#'
#' @description The function to calculate the AIC
#' @param log_like A logarithmic likelihood
#' @param n_par The number of parameters
calc_aic <- function(log_like, n_par){
  return(-2 * log_like + 2 * n_par)
}

#' estimate_parameters
#'
#' @description The function to estimate parameters
#' @param dt_exp The data.table of the experimental data
#' @param dt_gt The data.table of the genotype data
#' @param dt_seq The data.table of the sequence data
#' @param dt_kit_set The data.table of the kit-setting
#' @param n_parallel_chains The number of parallel chains
#' @param model_ae The stan model for AE
#' @param model_sr_gen The stan model for the general stutter ratio
#' @param model_sr_together The stan model for the stutter ratio when considering multiple loci together
#' @param model_sr_uni The stan model for the uniform stutter ratio
#' @param model_gamma The stan model for gamma distribution
#' @param n_warmup The number of simulations for warm-up
#' @param n_sampling The number of simulations for sampling
#' @param n_chains The number of chains
#' @param motif_lus The minimum number of repeat motifs for considering the LUS model
#' @param motif_multi The minimum number of repeat motifs for considering the multi-sequence model
estimate_parameters <- function(dt_exp, dt_gt, dt_seq, dt_kit_set,
                                model_ae, model_sr_gen, model_sr_together, model_sr_uni, model_gamma,
                                n_warmup, n_sampling, n_chains,
                                motif_lus, motif_multi, max_cand_us = 10,
                                n_parallel_chains = detectCores()){
  # Extract information on setting
  markers_par <- dt_kit_set[, Marker]
  n_mk <- length(markers_par)
  at_all_mk <- dt_kit_set[, AT]
  dye_all_mk <- dt_kit_set[, Dye]
  model_bsr_all_mk <- dt_kit_set[, Model_BSR]
  model_fsr_all_mk <- dt_kit_set[, Model_FSR]
  model_dsr_all_mk <- dt_kit_set[, Model_DSR]
  model_m2sr_all_mk <- dt_kit_set[, Model_M2SR]

  # Information on considering stutter ratios
  sr_consider_all_mk <- matrix(TRUE, n_mk, 4)
  rownames(sr_consider_all_mk) <- markers_par
  colnames(sr_consider_all_mk) <- c("BS", "FS", "DS", "M2S")
  sr_consider_all_mk[, 1] <- model_bsr_all_mk != "Not consider"
  sr_consider_all_mk[, 2] <- model_fsr_all_mk != "Not consider"
  sr_consider_all_mk[, 3] <- model_dsr_all_mk != "Not consider"
  sr_consider_all_mk[, 4] <- model_m2sr_all_mk != "Not consider"

  # Arrange experimental data
  tmp <- arrange_exp_data(dt_exp, dt_gt, markers_par, at_all_mk, dye_all_mk, sr_consider_all_mk, max_dist_pu = 1)
  peak_list <- tmp[[1]]
  size_list <- tmp[[2]]
  height_list <- tmp[[3]]
  gt_list <- tmp[[4]]

  # Prepare data for the locus-specific amplification efficiency
  tmp <- get_ae_data(height_list)
  ae_data <- tmp[[1]]
  height_ae_data <- tmp[[2]]

  # Prepare data for the stutter ratios
  tmp <- get_sr_data(peak_list, height_list, gt_list, sr_consider_all_mk)
  bsr_data <- tmp[[1]]
  fsr_data <- tmp[[2]]
  dsr_data <- tmp[[3]]
  m2sr_data <- tmp[[4]]
  allele_sr_data <- tmp[[5]]
  height_sr_data <- tmp[[6]]

  # Prepare data related to the sequence data
  tmp <- analyze_seq(dt_seq, markers_par)
  lus_all_mk <- tmp[[1]]
  reps_all_mk <- tmp[[2]]
  seq_al_all_mk <- tmp[[3]]
  seq_count_all_mk <- tmp[[4]]

  # Prepare data for the gamma distribution
  tph_all_data <- calc_tph_all_data(height_list)
  oph_all_data <- unlist(height_list)
  n_oph <- length(oph_all_data)

  # Perform the MCMC for the locus-specific amplification efficiency
  n_ae_data <- nrow(ae_data)
  ae_median_myu_all_mk <- ae_median_sigma_all_mk <- ae_rhat_myu_all_mk <- ae_rhat_sigma_all_mk <- ae_rhat_lp_all_mk <- rep(NA, n_mk)
  ae_trace_all_mk <- vector(mode = "list", length = n_mk)
  for(i in 1:n_mk){
    cat(paste0("MCMC for the locus-specific amplification efficiency of the ", markers_par[i], "\n"))
    data_list_ae <- list(N = n_ae_data, AE = ae_data[, i], Height = height_ae_data[, i])
    tmp <- run_mcmc_ae(model_ae, data_list_ae, n_warmup, n_sampling, n_chains, n_parallel_chains)
    ae_median_myu_all_mk[i] <- tmp$median_myu
    ae_median_sigma_all_mk[i] <- tmp$median_sigma
    ae_rhat_myu_all_mk[i] <- tmp$rhat_myu
    ae_rhat_sigma_all_mk[i] <- tmp$rhat_sigma
    ae_rhat_lp_all_mk[i] <- tmp$rhat_lp
    ae_trace_all_mk[[i]] <- tmp$traceplot
  }

  # Perform the MCMC for the back stutter ratio

  ## Define objects
  bsr_data_all_mk <- vector(mode = "list", length = n_mk)
  names(bsr_data_all_mk) <- markers_par
  bsr_alcor_all_mk <- vector(mode = "list", length = n_mk)
  names(bsr_alcor_all_mk) <- markers_par
  aics_bsr_all_mk <- vector(mode = "list", length = n_mk)
  names(aics_bsr_all_mk) <- markers_par
  bsr_allele_median_slope_all_mk <- bsr_allele_median_intercept_all_mk <- bsr_allele_median_sigma_all_mk <- rep(NA, n_mk)
  bsr_allele_rhat_slope_all_mk <- bsr_allele_rhat_intercept_all_mk <- bsr_allele_rhat_sigma_all_mk <- bsr_allele_rhat_lp_all_mk <- rep(NA, n_mk)
  bsr_allele_trace_all_mk <- vector(mode = "list", length = n_mk)
  bsr_lus_median_slope_all_mk <- bsr_lus_median_intercept_all_mk <- bsr_lus_median_sigma_all_mk <- rep(NA, n_mk)
  bsr_lus_rhat_slope_all_mk <- bsr_lus_rhat_intercept_all_mk <- bsr_lus_rhat_sigma_all_mk <- bsr_lus_rhat_lp_all_mk <- rep(NA, n_mk)
  bsr_lus_trace_all_mk <- vector(mode = "list", length = n_mk)
  bsr_multi_median_slope_all_mk <- bsr_multi_median_intercept_all_mk <- bsr_multi_median_sigma_all_mk <- bsr_multi_min_us_all_mk <- rep(NA, n_mk)
  bsr_multi_rhat_slope_all_mk <- bsr_multi_rhat_intercept_all_mk <- bsr_multi_rhat_sigma_all_mk <- bsr_multi_rhat_lp_all_mk <- rep(NA, n_mk)
  bsr_multi_trace_all_mk <- vector(mode = "list", length = n_mk)
  bsr_best_model_all_mk <- rep(NA, n_mk)
  bsr_best_median_myu_all_mk <- bsr_best_median_slope_all_mk <- bsr_best_median_intercept_all_mk <- bsr_best_median_sigma_all_mk <- bsr_best_min_us_all_mk <- rep(NA, n_mk)
  bsr_best_rhat_myu_all_mk <- bsr_best_rhat_slope_all_mk <- bsr_best_rhat_intercept_all_mk <- bsr_best_rhat_sigma_all_mk <- bsr_best_rhat_lp_all_mk <- rep(NA, n_mk)
  bsr_best_trace_all_mk <- vector(mode = "list", length = n_mk)

  ## BSR: locus-specific
  index_mk_specific_bsr <- which(is.element(model_bsr_all_mk, c("Best", "Allele", "LUS", "Multi-sequence")))
  n_mk_specific_bsr <- length(index_mk_specific_bsr)
  if(n_mk_specific_bsr > 0){
    for(i in 1:n_mk_specific_bsr){
      id_mk <- index_mk_specific_bsr[i]
      cat(paste0("MCMC for the back stutter ratio of the ", markers_par[id_mk], "\n"))

      # Extract data of one marker
      info_lus_one_mk <- lus_all_mk[[id_mk]]
      tmp <- extract_sr_use(bsr_data[, id_mk], allele_sr_data[, id_mk], height_sr_data[, id_mk], info_lus_one_mk)
      bsr_one_mk <- tmp[[1]]
      bsr_al_one_mk <- tmp[[2]]
      bsr_h_one_mk <- tmp[[3]]
      bsr_lus_one_mk <- tmp[[4]]
      reps_one_mk <- reps_all_mk[[id_mk]]
      seq_al_one_mk <- seq_al_all_mk[[id_mk]]
      seq_count_one_mk <- seq_count_all_mk[[id_mk]]

      # Set consideration of each model
      model_bsr <- model_bsr_all_mk[id_mk]
      consider_bsr <- set_consider_sr_model(model_bsr, reps_one_mk, motif_lus, motif_multi)

      # Define an object to assign log-likelihoods of each model
      log_likes <- rep(NA, 3)
      names(log_likes) <- c("Allele", "LUS", "Multi-sequence")

      # Allele model
      if(consider_bsr["Allele"]){
        cat("Allele model\n")
        data_list_bsr_allele <- list(N = length(bsr_one_mk), SR = bsr_one_mk, XVal = bsr_al_one_mk, Height = bsr_h_one_mk)
        tmp <- run_mcmc_sr_gen(model_sr_gen, data_list_bsr_allele, n_warmup, n_sampling, n_chains, n_parallel_chains)
        log_likes[1] <- tmp$median_loglike
        bsr_allele_median_slope_all_mk[id_mk] <- tmp$median_slope
        bsr_allele_median_intercept_all_mk[id_mk] <- tmp$median_intercept
        bsr_allele_median_sigma_all_mk[id_mk] <- tmp$median_sigma
        bsr_allele_rhat_slope_all_mk[id_mk] <- tmp$rhat_slope
        bsr_allele_rhat_intercept_all_mk[id_mk] <- tmp$rhat_intercept
        bsr_allele_rhat_sigma_all_mk[id_mk] <- tmp$rhat_sigma
        bsr_allele_rhat_lp_all_mk[id_mk] <- tmp$rhat_lp
        bsr_allele_trace_all_mk[[id_mk]] <- tmp$traceplot
      }

      # LUS model
      if(consider_bsr["LUS"]){
        cat("LUS model\n")
        data_list_bsr_lus <- list(N = length(bsr_one_mk), SR = bsr_one_mk, XVal = bsr_lus_one_mk, Height = bsr_h_one_mk)
        tmp <- run_mcmc_sr_gen(model_sr_gen, data_list_bsr_lus, n_warmup, n_sampling, n_chains, n_parallel_chains)
        log_likes[2] <- tmp$median_loglike
        bsr_lus_median_slope_all_mk[id_mk] <- tmp$median_slope
        bsr_lus_median_intercept_all_mk[id_mk] <- tmp$median_intercept
        bsr_lus_median_sigma_all_mk[id_mk] <- tmp$median_sigma
        bsr_lus_rhat_slope_all_mk[id_mk] <- tmp$rhat_slope
        bsr_lus_rhat_intercept_all_mk[id_mk] <- tmp$rhat_intercept
        bsr_lus_rhat_sigma_all_mk[id_mk] <- tmp$rhat_sigma
        bsr_lus_rhat_lp_all_mk[id_mk] <- tmp$rhat_lp
        bsr_lus_trace_all_mk[[id_mk]] <- tmp$traceplot
      }

      # Multi-sequence model
      if(consider_bsr["Multi-sequence"]){
        cat("Multi-sequence model\n")
        cand_min_us <- sort(unique(unlist(reps_one_mk)))
        cand_min_us <- cand_min_us[cand_min_us <= max_cand_us]
        n_cand_min_us <- length(cand_min_us)
        loglike_multi <- rep(NA, n_cand_min_us)
        bsr_ca_one_mk_cand <- result_bsr_multi_cand <- list()
        for(j in 1:n_cand_min_us){
          bsr_ca_one_mk_cand[[j]] <- bsr_ca_one_mk <- make_sr_ca_one_mk(cand_min_us[j], reps_one_mk, seq_al_one_mk, seq_count_one_mk, bsr_al_one_mk)
          data_list_bsr_multi <- list(N = length(bsr_one_mk), SR = bsr_one_mk, XVal = bsr_ca_one_mk, Height = bsr_h_one_mk)
          result_bsr_multi_cand[[j]] <- run_mcmc_sr_gen(model_sr_gen, data_list_bsr_multi, n_warmup, n_sampling, n_chains, n_parallel_chains)
          loglike_multi[j] <- result_bsr_multi_cand[[j]]$median_loglike
        }
        index_max_like_multi <- which.max(loglike_multi)
        bsr_ca_one_mk <- bsr_ca_one_mk_cand[[index_max_like_multi]]
        tmp <- result_bsr_multi_cand[[index_max_like_multi]]
        log_likes[3] <- tmp$median_loglike
        bsr_multi_median_slope_all_mk[id_mk] <- tmp$median_slope
        bsr_multi_median_intercept_all_mk[id_mk] <- tmp$median_intercept
        bsr_multi_median_sigma_all_mk[id_mk] <- tmp$median_sigma
        bsr_multi_min_us_all_mk[id_mk] <- cand_min_us[index_max_like_multi]
        bsr_multi_rhat_slope_all_mk[id_mk] <- tmp$rhat_slope
        bsr_multi_rhat_intercept_all_mk[id_mk] <- tmp$rhat_intercept
        bsr_multi_rhat_sigma_all_mk[id_mk] <- tmp$rhat_sigma
        bsr_multi_rhat_lp_all_mk[id_mk] <- tmp$rhat_lp
        bsr_multi_trace_all_mk[[id_mk]] <- tmp$traceplot
      }

      # Best model
      index_cand <- which(!is.na(log_likes))
      log_likes <- log_likes[index_cand]
      n_pars <- c(3, 3, 4)
      n_pars <- n_pars[index_cand]
      aics <- mapply(calc_aic, log_likes, n_pars)
      names(aics) <- names(log_likes)
      index_best_model <- which.min(aics)
      if(index_best_model == 1){
        bsr_best_model_all_mk[id_mk] <- "Allele"
        bsr_best_median_slope_all_mk[id_mk] <- bsr_allele_median_slope_all_mk[id_mk]
        bsr_best_median_intercept_all_mk[id_mk] <- bsr_allele_median_intercept_all_mk[id_mk]
        bsr_best_median_sigma_all_mk[id_mk] <- bsr_allele_median_sigma_all_mk[id_mk]
        bsr_best_rhat_slope_all_mk[id_mk] <- bsr_allele_rhat_slope_all_mk[id_mk]
        bsr_best_rhat_intercept_all_mk[id_mk] <- bsr_allele_rhat_intercept_all_mk[id_mk]
        bsr_best_rhat_sigma_all_mk[id_mk] <- bsr_allele_rhat_sigma_all_mk[id_mk]
        bsr_best_rhat_lp_all_mk[id_mk] <- bsr_allele_rhat_lp_all_mk[id_mk]
        bsr_best_trace_all_mk[[id_mk]] <- bsr_allele_trace_all_mk[[id_mk]]
      }else if(index_best_model == 2){
        bsr_best_model_all_mk[id_mk] <- "LUS"
        bsr_best_median_slope_all_mk[id_mk] <- bsr_lus_median_slope_all_mk[id_mk]
        bsr_best_median_intercept_all_mk[id_mk] <- bsr_lus_median_intercept_all_mk[id_mk]
        bsr_best_median_sigma_all_mk[id_mk] <- bsr_lus_median_sigma_all_mk[id_mk]
        bsr_best_rhat_slope_all_mk[id_mk] <- bsr_lus_rhat_slope_all_mk[id_mk]
        bsr_best_rhat_intercept_all_mk[id_mk] <- bsr_lus_rhat_intercept_all_mk[id_mk]
        bsr_best_rhat_sigma_all_mk[id_mk] <- bsr_lus_rhat_sigma_all_mk[id_mk]
        bsr_best_rhat_lp_all_mk[id_mk] <- bsr_lus_rhat_lp_all_mk[id_mk]
        bsr_best_trace_all_mk[[id_mk]] <- bsr_lus_trace_all_mk[[id_mk]]
      }else{
        bsr_best_model_all_mk[id_mk] <- "Multi-sequence"
        bsr_best_median_slope_all_mk[id_mk] <- bsr_multi_median_slope_all_mk[id_mk]
        bsr_best_median_intercept_all_mk[id_mk] <- bsr_multi_median_intercept_all_mk[id_mk]
        bsr_best_median_sigma_all_mk[id_mk] <- bsr_multi_median_sigma_all_mk[id_mk]
        bsr_best_min_us_all_mk[id_mk] <- bsr_multi_min_us_all_mk[id_mk]
        bsr_best_rhat_slope_all_mk[id_mk] <- bsr_multi_rhat_slope_all_mk[id_mk]
        bsr_best_rhat_intercept_all_mk[id_mk] <- bsr_multi_rhat_intercept_all_mk[id_mk]
        bsr_best_rhat_sigma_all_mk[id_mk] <- bsr_multi_rhat_sigma_all_mk[id_mk]
        bsr_best_rhat_lp_all_mk[id_mk] <- bsr_multi_rhat_lp_all_mk[id_mk]
        bsr_best_trace_all_mk[[id_mk]] <- bsr_multi_trace_all_mk[[id_mk]]
      }

      # Assign data
      bsr_data_one_mk <- vector(mode = "list", length = 4)
      names(bsr_data_one_mk) <- c("bsr", "bsr_al", "bsr_lus", "bsr_ca")
      bsr_data_one_mk[[1]] <- bsr_one_mk
      if(consider_bsr["Allele"]){
        bsr_data_one_mk[[2]] <- bsr_al_one_mk
      }
      if(consider_bsr["LUS"]){
        bsr_data_one_mk[[3]] <- bsr_lus_one_mk
      }
      if(consider_bsr["Multi-sequence"]){
        bsr_data_one_mk[[4]] <- bsr_ca_one_mk
      }
      bsr_data_all_mk[[id_mk]] <- bsr_data_one_mk
      aics_bsr_all_mk[[id_mk]] <- aics

      # Allele repeat correction
      info_al_one_mk <- as.numeric(names(info_lus_one_mk))
      if(consider_bsr["Multi-sequence"]){
        info_ca_one_mk <- make_sr_ca_one_mk(cand_min_us[index_max_like_multi], reps_one_mk, seq_al_one_mk, seq_count_one_mk, info_al_one_mk)
      }else{
        info_ca_one_mk <- NULL
      }
      bsr_alcor_all_mk[[id_mk]] <- list(info_al = info_al_one_mk,
                                        info_lus = info_lus_one_mk,
                                        info_ca = info_ca_one_mk)
    }
  }

  ## BSR: multiple loci together
  index_mk_together_bsr <- which(model_bsr_all_mk == "Multiple loci together")
  n_mk_together_bsr <- length(index_mk_together_bsr)
  bsr_together <- bsr_h_together <- id_mk_together <- numeric(0)
  id_mk_count <- 0
  if(n_mk_together_bsr > 0){
    cat("MCMC for the back stutter ratio of the multiple loci together\n")

    # Prepare data of the target loci
    for(i in 1:n_mk_together_bsr){
      # Extract data of one marker
      id_mk <- index_mk_together_bsr[i]
      tmp <- extract_sr_use(bsr_data[, id_mk], allele_sr_data[, id_mk], height_sr_data[, id_mk], lus_all_mk[[id_mk]])
      bsr_one_mk <- tmp[[1]]
      bsr_al_one_mk <- tmp[[2]]
      bsr_h_one_mk <- tmp[[3]]

      # Combine data
      bsr_together <- c(bsr_together, bsr_one_mk)
      bsr_h_together <- c(bsr_h_together, bsr_h_one_mk)

      # Make id_mk_together
      id_mk_count <- id_mk_count + 1
      id_mk_together <- c(id_mk_together, rep(id_mk_count, length(bsr_one_mk)))

      # Assign data
      bsr_data_one_mk <- vector(mode = "list", length = 4)
      names(bsr_data_one_mk) <- c("bsr", "bsr_al", "bsr_lus", "bsr_ca")
      bsr_data_one_mk[[1]] <- bsr_one_mk
      bsr_data_one_mk[[2]] <- bsr_al_one_mk
      bsr_data_all_mk[[id_mk]] <- bsr_data_one_mk
    }
    data_list_bsr_together <- list(N = length(bsr_together), NMk = n_mk_together_bsr, SR = bsr_together, Height = bsr_h_together, IndexMk = id_mk_together)
    tmp <- run_mcmc_sr_together(model_sr_together, data_list_bsr_together, n_warmup, n_sampling, n_chains, n_parallel_chains)
    bsr_best_model_all_mk[index_mk_together_bsr] <- "Multiple loci together"
    bsr_best_median_myu_all_mk[index_mk_together_bsr] <- tmp$median_myus
    bsr_best_median_sigma_all_mk[index_mk_together_bsr] <- tmp$median_sigma
    bsr_best_rhat_myu_all_mk[index_mk_together_bsr] <- tmp$rhat_myus
    bsr_best_rhat_sigma_all_mk[index_mk_together_bsr] <- tmp$rhat_sigma
    bsr_best_rhat_lp_all_mk[index_mk_together_bsr] <- tmp$rhat_lp
    bsr_best_trace_all_mk[index_mk_together_bsr] <- tmp$traceplot
  }

  ## BSR: not consider
  index_mk_not_consider_bsr <- which(model_bsr_all_mk == "Not consider")
  if(length(index_mk_not_consider_bsr) > 0){
    bsr_best_model_all_mk[index_mk_not_consider_bsr] <- "Not consider"
  }

  # Perform the MCMC for the forward stutter ratio

  ## Define objects
  fsr_data_all_mk <- vector(mode = "list", length = n_mk)
  names(fsr_data_all_mk) <- markers_par
  fsr_alcor_all_mk <- vector(mode = "list", length = n_mk)
  names(fsr_alcor_all_mk) <- markers_par
  aics_fsr_all_mk <- vector(mode = "list", length = n_mk)
  names(aics_fsr_all_mk) <- markers_par
  fsr_allele_median_slope_all_mk <- fsr_allele_median_intercept_all_mk <- fsr_allele_median_sigma_all_mk <- rep(NA, n_mk)
  fsr_allele_rhat_slope_all_mk <- fsr_allele_rhat_intercept_all_mk <- fsr_allele_rhat_sigma_all_mk <- fsr_allele_rhat_lp_all_mk <- rep(NA, n_mk)
  fsr_allele_trace_all_mk <- vector(mode = "list", length = n_mk)
  fsr_lus_median_slope_all_mk <- fsr_lus_median_intercept_all_mk <- fsr_lus_median_sigma_all_mk <- rep(NA, n_mk)
  fsr_lus_rhat_slope_all_mk <- fsr_lus_rhat_intercept_all_mk <- fsr_lus_rhat_sigma_all_mk <- fsr_lus_rhat_lp_all_mk <- rep(NA, n_mk)
  fsr_lus_trace_all_mk <- vector(mode = "list", length = n_mk)
  fsr_multi_median_slope_all_mk <- fsr_multi_median_intercept_all_mk <- fsr_multi_median_sigma_all_mk <- fsr_multi_min_us_all_mk <- rep(NA, n_mk)
  fsr_multi_rhat_slope_all_mk <- fsr_multi_rhat_intercept_all_mk <- fsr_multi_rhat_sigma_all_mk <- fsr_multi_rhat_lp_all_mk <- rep(NA, n_mk)
  fsr_multi_trace_all_mk <- vector(mode = "list", length = n_mk)
  fsr_best_model_all_mk <- rep(NA, n_mk)
  fsr_best_median_myu_all_mk <- fsr_best_median_slope_all_mk <- fsr_best_median_intercept_all_mk <- fsr_best_median_sigma_all_mk <- fsr_best_min_us_all_mk <- rep(NA, n_mk)
  fsr_best_rhat_myu_all_mk <- fsr_best_rhat_slope_all_mk <- fsr_best_rhat_intercept_all_mk <- fsr_best_rhat_sigma_all_mk <- fsr_best_rhat_lp_all_mk <- rep(NA, n_mk)
  fsr_best_trace_all_mk <- vector(mode = "list", length = n_mk)

  ## FSR: locus-specific
  index_mk_specific_fsr <- which(is.element(model_fsr_all_mk, c("Best", "Allele", "LUS", "Multi-sequence")))
  n_mk_specific_fsr <- length(index_mk_specific_fsr)
  if(n_mk_specific_fsr > 0){
    for(i in 1:n_mk_specific_fsr){
      id_mk <- index_mk_specific_fsr[i]
      cat(paste0("MCMC for the forward stutter ratio of the ", markers_par[id_mk], "\n"))

      # Extract data of one marker
      info_lus_one_mk <- lus_all_mk[[id_mk]]
      tmp <- extract_sr_use(fsr_data[, id_mk], allele_sr_data[, id_mk], height_sr_data[, id_mk], info_lus_one_mk)
      fsr_one_mk <- tmp[[1]]
      fsr_al_one_mk <- tmp[[2]]
      fsr_h_one_mk <- tmp[[3]]
      fsr_lus_one_mk <- tmp[[4]]
      reps_one_mk <- reps_all_mk[[id_mk]]
      seq_al_one_mk <- seq_al_all_mk[[id_mk]]
      seq_count_one_mk <- seq_count_all_mk[[id_mk]]

      # Set consideration of each model
      model_fsr <- model_fsr_all_mk[id_mk]
      consider_fsr <- set_consider_sr_model(model_fsr, reps_one_mk, motif_lus, motif_multi)

      # Define an object to assign log-likelihoods of each model
      log_likes <- rep(NA, 3)
      names(log_likes) <- c("Allele", "LUS", "Multi-sequence")

      # Allele model
      if(consider_fsr["Allele"]){
        cat("Allele model\n")
        data_list_fsr_allele <- list(N = length(fsr_one_mk), SR = fsr_one_mk, XVal = fsr_al_one_mk, Height = fsr_h_one_mk)
        tmp <- run_mcmc_sr_gen(model_sr_gen, data_list_fsr_allele, n_warmup, n_sampling, n_chains, n_parallel_chains)
        log_likes[1] <- tmp$median_loglike
        fsr_allele_median_slope_all_mk[id_mk] <- tmp$median_slope
        fsr_allele_median_intercept_all_mk[id_mk] <- tmp$median_intercept
        fsr_allele_median_sigma_all_mk[id_mk] <- tmp$median_sigma
        fsr_allele_rhat_slope_all_mk[id_mk] <- tmp$rhat_slope
        fsr_allele_rhat_intercept_all_mk[id_mk] <- tmp$rhat_intercept
        fsr_allele_rhat_sigma_all_mk[id_mk] <- tmp$rhat_sigma
        fsr_allele_rhat_lp_all_mk[id_mk] <- tmp$rhat_lp
        fsr_allele_trace_all_mk[[id_mk]] <- tmp$traceplot
      }

      # LUS model
      if(consider_fsr["LUS"]){
        cat("LUS model\n")
        data_list_fsr_lus <- list(N = length(fsr_one_mk), SR = fsr_one_mk, XVal = fsr_lus_one_mk, Height = fsr_h_one_mk)
        tmp <- run_mcmc_sr_gen(model_sr_gen, data_list_fsr_lus, n_warmup, n_sampling, n_chains, n_parallel_chains)
        log_likes[2] <- tmp$median_loglike
        fsr_lus_median_slope_all_mk[id_mk] <- tmp$median_slope
        fsr_lus_median_intercept_all_mk[id_mk] <- tmp$median_intercept
        fsr_lus_median_sigma_all_mk[id_mk] <- tmp$median_sigma
        fsr_lus_rhat_slope_all_mk[id_mk] <- tmp$rhat_slope
        fsr_lus_rhat_intercept_all_mk[id_mk] <- tmp$rhat_intercept
        fsr_lus_rhat_sigma_all_mk[id_mk] <- tmp$rhat_sigma
        fsr_lus_rhat_lp_all_mk[id_mk] <- tmp$rhat_lp
        fsr_lus_trace_all_mk[[id_mk]] <- tmp$traceplot
      }

      # Multi-sequence model
      if(consider_fsr["Multi-sequence"]){
        cat("Multi-sequence model\n")
        cand_min_us <- sort(unique(unlist(reps_one_mk)))
        cand_min_us <- cand_min_us[cand_min_us <= max_cand_us]
        n_cand_min_us <- length(cand_min_us)
        loglike_multi <- rep(NA, n_cand_min_us)
        fsr_ca_one_mk_cand <- result_fsr_multi_cand <- list()
        for(j in 1:n_cand_min_us){
          fsr_ca_one_mk_cand[[j]] <- fsr_ca_one_mk <- make_sr_ca_one_mk(cand_min_us[j], reps_one_mk, seq_al_one_mk, seq_count_one_mk, fsr_al_one_mk)
          data_list_fsr_multi <- list(N = length(fsr_one_mk), SR = fsr_one_mk, XVal = fsr_ca_one_mk, Height = fsr_h_one_mk)
          result_fsr_multi_cand[[j]] <- run_mcmc_sr_gen(model_sr_gen, data_list_fsr_multi, n_warmup, n_sampling, n_chains, n_parallel_chains)
          loglike_multi[j] <- result_fsr_multi_cand[[j]]$median_loglike
        }
        index_max_like_multi <- which.max(loglike_multi)
        fsr_ca_one_mk <- fsr_ca_one_mk_cand[[index_max_like_multi]]
        tmp <- result_fsr_multi_cand[[index_max_like_multi]]
        log_likes[3] <- tmp$median_loglike
        fsr_multi_median_slope_all_mk[id_mk] <- tmp$median_slope
        fsr_multi_median_intercept_all_mk[id_mk] <- tmp$median_intercept
        fsr_multi_median_sigma_all_mk[id_mk] <- tmp$median_sigma
        fsr_multi_min_us_all_mk[id_mk] <- cand_min_us[index_max_like_multi]
        fsr_multi_rhat_slope_all_mk[id_mk] <- tmp$rhat_slope
        fsr_multi_rhat_intercept_all_mk[id_mk] <- tmp$rhat_intercept
        fsr_multi_rhat_sigma_all_mk[id_mk] <- tmp$rhat_sigma
        fsr_multi_rhat_lp_all_mk[id_mk] <- tmp$rhat_lp
        fsr_multi_trace_all_mk[[id_mk]] <- tmp$traceplot
      }

      # Best model
      index_cand <- which(!is.na(log_likes))
      log_likes <- log_likes[index_cand]
      n_pars <- c(3, 3, 4)
      n_pars <- n_pars[index_cand]
      aics <- mapply(calc_aic, log_likes, n_pars)
      names(aics) <- names(log_likes)
      index_best_model <- which.min(aics)
      if(index_best_model == 1){
        fsr_best_model_all_mk[id_mk] <- "Allele"
        fsr_best_median_slope_all_mk[id_mk] <- fsr_allele_median_slope_all_mk[id_mk]
        fsr_best_median_intercept_all_mk[id_mk] <- fsr_allele_median_intercept_all_mk[id_mk]
        fsr_best_median_sigma_all_mk[id_mk] <- fsr_allele_median_sigma_all_mk[id_mk]
        fsr_best_rhat_slope_all_mk[id_mk] <- fsr_allele_rhat_slope_all_mk[id_mk]
        fsr_best_rhat_intercept_all_mk[id_mk] <- fsr_allele_rhat_intercept_all_mk[id_mk]
        fsr_best_rhat_sigma_all_mk[id_mk] <- fsr_allele_rhat_sigma_all_mk[id_mk]
        fsr_best_rhat_lp_all_mk[id_mk] <- fsr_allele_rhat_lp_all_mk[id_mk]
        fsr_best_trace_all_mk[[id_mk]] <- fsr_allele_trace_all_mk[[id_mk]]
      }else if(index_best_model == 2){
        fsr_best_model_all_mk[id_mk] <- "LUS"
        fsr_best_median_slope_all_mk[id_mk] <- fsr_lus_median_slope_all_mk[id_mk]
        fsr_best_median_intercept_all_mk[id_mk] <- fsr_lus_median_intercept_all_mk[id_mk]
        fsr_best_median_sigma_all_mk[id_mk] <- fsr_lus_median_sigma_all_mk[id_mk]
        fsr_best_rhat_slope_all_mk[id_mk] <- fsr_lus_rhat_slope_all_mk[id_mk]
        fsr_best_rhat_intercept_all_mk[id_mk] <- fsr_lus_rhat_intercept_all_mk[id_mk]
        fsr_best_rhat_sigma_all_mk[id_mk] <- fsr_lus_rhat_sigma_all_mk[id_mk]
        fsr_best_rhat_lp_all_mk[id_mk] <- fsr_lus_rhat_lp_all_mk[id_mk]
        fsr_best_trace_all_mk[[id_mk]] <- fsr_lus_trace_all_mk[[id_mk]]
      }else{
        fsr_best_model_all_mk[id_mk] <- "Multi-sequence"
        fsr_best_median_slope_all_mk[id_mk] <- fsr_multi_median_slope_all_mk[id_mk]
        fsr_best_median_intercept_all_mk[id_mk] <- fsr_multi_median_intercept_all_mk[id_mk]
        fsr_best_median_sigma_all_mk[id_mk] <- fsr_multi_median_sigma_all_mk[id_mk]
        fsr_best_min_us_all_mk[id_mk] <- fsr_multi_min_us_all_mk[id_mk]
        fsr_best_rhat_slope_all_mk[id_mk] <- fsr_multi_rhat_slope_all_mk[id_mk]
        fsr_best_rhat_intercept_all_mk[id_mk] <- fsr_multi_rhat_intercept_all_mk[id_mk]
        fsr_best_rhat_sigma_all_mk[id_mk] <- fsr_multi_rhat_sigma_all_mk[id_mk]
        fsr_best_rhat_lp_all_mk[id_mk] <- fsr_multi_rhat_lp_all_mk[id_mk]
        fsr_best_trace_all_mk[[id_mk]] <- fsr_multi_trace_all_mk[[id_mk]]
      }

      # Assign data
      fsr_data_one_mk <- vector(mode = "list", length = 4)
      names(fsr_data_one_mk) <- c("fsr", "fsr_al", "fsr_lus", "fsr_ca")
      fsr_data_one_mk[[1]] <- fsr_one_mk
      if(consider_fsr["Allele"]){
        fsr_data_one_mk[[2]] <- fsr_al_one_mk
      }
      if(consider_fsr["LUS"]){
        fsr_data_one_mk[[3]] <- fsr_lus_one_mk
      }
      if(consider_fsr["Multi-sequence"]){
        fsr_data_one_mk[[4]] <- fsr_ca_one_mk
      }
      fsr_data_all_mk[[id_mk]] <- fsr_data_one_mk
      aics_fsr_all_mk[[id_mk]] <- aics

      # Allele repeat correction
      info_al_one_mk <- as.numeric(names(info_lus_one_mk))
      if(consider_fsr["Multi-sequence"]){
        info_ca_one_mk <- make_sr_ca_one_mk(cand_min_us[index_max_like_multi], reps_one_mk, seq_al_one_mk, seq_count_one_mk, info_al_one_mk)
      }else{
        info_ca_one_mk <- NULL
      }
      fsr_alcor_all_mk[[id_mk]] <- list(info_al = info_al_one_mk,
                                        info_lus = info_lus_one_mk,
                                        info_ca = info_ca_one_mk)
    }
  }

  ## FSR: multiple loci together
  index_mk_together_fsr <- which(model_fsr_all_mk == "Multiple loci together")
  n_mk_together_fsr <- length(index_mk_together_fsr)
  fsr_together <- fsr_h_together <- id_mk_together <- numeric(0)
  id_mk_count <- 0
  if(n_mk_together_fsr > 0){
    cat("MCMC for the forward stutter ratio of the multiple loci together\n")

    # Prepare data of the target loci
    for(i in 1:n_mk_together_fsr){
      # Extract data of one marker
      id_mk <- index_mk_together_fsr[i]
      tmp <- extract_sr_use(fsr_data[, id_mk], allele_sr_data[, id_mk], height_sr_data[, id_mk], lus_all_mk[[id_mk]])
      fsr_one_mk <- tmp[[1]]
      fsr_al_one_mk <- tmp[[2]]
      fsr_h_one_mk <- tmp[[3]]

      # Combine data
      fsr_together <- c(fsr_together, fsr_one_mk)
      fsr_h_together <- c(fsr_h_together, fsr_h_one_mk)

      # Make id_mk_together
      id_mk_count <- id_mk_count + 1
      id_mk_together <- c(id_mk_together, rep(id_mk_count, length(fsr_one_mk)))

      # Assign data
      fsr_data_one_mk <- vector(mode = "list", length = 4)
      names(fsr_data_one_mk) <- c("fsr", "fsr_al", "fsr_lus", "fsr_ca")
      fsr_data_one_mk[[1]] <- fsr_one_mk
      fsr_data_one_mk[[2]] <- fsr_al_one_mk
      fsr_data_all_mk[[id_mk]] <- fsr_data_one_mk
    }
    data_list_fsr_together <- list(N = length(fsr_together), NMk = n_mk_together_fsr, SR = fsr_together, Height = fsr_h_together, IndexMk = id_mk_together)
    tmp <- run_mcmc_sr_together(model_sr_together, data_list_fsr_together, n_warmup, n_sampling, n_chains, n_parallel_chains)
    fsr_best_model_all_mk[index_mk_together_fsr] <- "Multiple loci together"
    fsr_best_median_myu_all_mk[index_mk_together_fsr] <- tmp$median_myus
    fsr_best_median_sigma_all_mk[index_mk_together_fsr] <- tmp$median_sigma
    fsr_best_rhat_myu_all_mk[index_mk_together_fsr] <- tmp$rhat_myus
    fsr_best_rhat_sigma_all_mk[index_mk_together_fsr] <- tmp$rhat_sigma
    fsr_best_rhat_lp_all_mk[index_mk_together_fsr] <- tmp$rhat_lp
    fsr_best_trace_all_mk[index_mk_together_fsr] <- tmp$traceplot
  }

  ## FSR: not consider
  index_mk_not_consider_fsr <- which(model_fsr_all_mk == "Not consider")
  if(length(index_mk_not_consider_fsr) > 0){
    fsr_best_model_all_mk[index_mk_not_consider_fsr] <- "Not consider"
  }

  # Perform the MCMC for the double-back stutter ratio

  ## Define objects
  dsr_data_all_mk <- vector(mode = "list", length = n_mk)
  names(dsr_data_all_mk) <- markers_par
  dsr_alcor_all_mk <- vector(mode = "list", length = n_mk)
  names(dsr_alcor_all_mk) <- markers_par
  aics_dsr_all_mk <- vector(mode = "list", length = n_mk)
  names(aics_dsr_all_mk) <- markers_par
  dsr_allele_median_slope_all_mk <- dsr_allele_median_intercept_all_mk <- dsr_allele_median_sigma_all_mk <- rep(NA, n_mk)
  dsr_allele_rhat_slope_all_mk <- dsr_allele_rhat_intercept_all_mk <- dsr_allele_rhat_sigma_all_mk <- dsr_allele_rhat_lp_all_mk <- rep(NA, n_mk)
  dsr_allele_trace_all_mk <- vector(mode = "list", length = n_mk)
  dsr_lus_median_slope_all_mk <- dsr_lus_median_intercept_all_mk <- dsr_lus_median_sigma_all_mk <- rep(NA, n_mk)
  dsr_lus_rhat_slope_all_mk <- dsr_lus_rhat_intercept_all_mk <- dsr_lus_rhat_sigma_all_mk <- dsr_lus_rhat_lp_all_mk <- rep(NA, n_mk)
  dsr_lus_trace_all_mk <- vector(mode = "list", length = n_mk)
  dsr_multi_median_slope_all_mk <- dsr_multi_median_intercept_all_mk <- dsr_multi_median_sigma_all_mk <- dsr_multi_min_us_all_mk <- rep(NA, n_mk)
  dsr_multi_rhat_slope_all_mk <- dsr_multi_rhat_intercept_all_mk <- dsr_multi_rhat_sigma_all_mk <- dsr_multi_rhat_lp_all_mk <- rep(NA, n_mk)
  dsr_multi_trace_all_mk <- vector(mode = "list", length = n_mk)
  dsr_best_model_all_mk <- rep(NA, n_mk)
  dsr_best_median_myu_all_mk <- dsr_best_median_slope_all_mk <- dsr_best_median_intercept_all_mk <- dsr_best_median_sigma_all_mk <- dsr_best_min_us_all_mk <- rep(NA, n_mk)
  dsr_best_rhat_myu_all_mk <- dsr_best_rhat_slope_all_mk <- dsr_best_rhat_intercept_all_mk <- dsr_best_rhat_sigma_all_mk <- dsr_best_rhat_lp_all_mk <- rep(NA, n_mk)
  dsr_best_trace_all_mk <- vector(mode = "list", length = n_mk)

  ## DSR: locus-specific
  index_mk_specific_dsr <- which(is.element(model_dsr_all_mk, c("Best", "Allele", "LUS", "Multi-sequence")))
  n_mk_specific_dsr <- length(index_mk_specific_dsr)
  if(n_mk_specific_dsr > 0){
    for(i in 1:n_mk_specific_dsr){
      id_mk <- index_mk_specific_dsr[i]
      cat(paste0("MCMC for the double-back stutter ratio of the ", markers_par[id_mk], "\n"))

      # Extract data of one marker
      info_lus_one_mk <- lus_all_mk[[id_mk]]
      tmp <- extract_sr_use(dsr_data[, id_mk], allele_sr_data[, id_mk], height_sr_data[, id_mk], info_lus_one_mk)
      dsr_one_mk <- tmp[[1]]
      dsr_al_one_mk <- tmp[[2]]
      dsr_h_one_mk <- tmp[[3]]
      dsr_lus_one_mk <- tmp[[4]]
      reps_one_mk <- reps_all_mk[[id_mk]]
      seq_al_one_mk <- seq_al_all_mk[[id_mk]]
      seq_count_one_mk <- seq_count_all_mk[[id_mk]]

      # Set consideration of each model
      model_dsr <- model_dsr_all_mk[id_mk]
      consider_dsr <- set_consider_sr_model(model_dsr, reps_one_mk, motif_lus, motif_multi)

      # Define an object to assign log-likelihoods of each model
      log_likes <- rep(NA, 3)
      names(log_likes) <- c("Allele", "LUS", "Multi-sequence")

      # Allele model
      if(model_dsr["Allele"]){
        cat("Allele model\n")
        data_list_dsr_allele <- list(N = length(dsr_one_mk), SR = dsr_one_mk, XVal = dsr_al_one_mk, Height = dsr_h_one_mk)
        tmp <- run_mcmc_sr_gen(model_sr_gen, data_list_dsr_allele, n_warmup, n_sampling, n_chains, n_parallel_chains)
        log_likes[1] <- tmp$median_loglike
        dsr_allele_median_slope_all_mk[id_mk] <- tmp$median_slope
        dsr_allele_median_intercept_all_mk[id_mk] <- tmp$median_intercept
        dsr_allele_median_sigma_all_mk[id_mk] <- tmp$median_sigma
        dsr_allele_rhat_slope_all_mk[id_mk] <- tmp$rhat_slope
        dsr_allele_rhat_intercept_all_mk[id_mk] <- tmp$rhat_intercept
        dsr_allele_rhat_sigma_all_mk[id_mk] <- tmp$rhat_sigma
        dsr_allele_rhat_lp_all_mk[id_mk] <- tmp$rhat_lp
        dsr_allele_trace_all_mk[[id_mk]] <- tmp$traceplot
      }

      # LUS model
      if(model_dsr["LUS"]){
        cat("LUS model\n")
        data_list_dsr_lus <- list(N = length(dsr_one_mk), SR = dsr_one_mk, XVal = dsr_lus_one_mk, Height = dsr_h_one_mk)
        tmp <- run_mcmc_sr_gen(model_sr_gen, data_list_dsr_lus, n_warmup, n_sampling, n_chains, n_parallel_chains)
        log_likes[2] <- tmp$median_loglike
        dsr_lus_median_slope_all_mk[id_mk] <- tmp$median_slope
        dsr_lus_median_intercept_all_mk[id_mk] <- tmp$median_intercept
        dsr_lus_median_sigma_all_mk[id_mk] <- tmp$median_sigma
        dsr_lus_rhat_slope_all_mk[id_mk] <- tmp$rhat_slope
        dsr_lus_rhat_intercept_all_mk[id_mk] <- tmp$rhat_intercept
        dsr_lus_rhat_sigma_all_mk[id_mk] <- tmp$rhat_sigma
        dsr_lus_rhat_lp_all_mk[id_mk] <- tmp$rhat_lp
        dsr_lus_trace_all_mk[[id_mk]] <- tmp$traceplot
      }

      # Multi-sequence model
      if(model_dsr["Multi-sequence"]){
        cat("Multi-sequence model\n")
        cand_min_us <- sort(unique(unlist(reps_one_mk)))
        cand_min_us <- cand_min_us[cand_min_us <= max_cand_us]
        n_cand_min_us <- length(cand_min_us)
        loglike_multi <- rep(NA, n_cand_min_us)
        dsr_ca_one_mk_cand <- result_dsr_multi_cand <- list()
        for(j in 1:n_cand_min_us){
          dsr_ca_one_mk_cand[[j]] <- dsr_ca_one_mk <- make_sr_ca_one_mk(cand_min_us[j], reps_one_mk, seq_al_one_mk, seq_count_one_mk, dsr_al_one_mk)
          data_list_dsr_multi <- list(N = length(dsr_one_mk), SR = dsr_one_mk, XVal = dsr_ca_one_mk, Height = dsr_h_one_mk)
          result_dsr_multi_cand[[j]] <- run_mcmc_sr_gen(model_sr_gen, data_list_dsr_multi, n_warmup, n_sampling, n_chains, n_parallel_chains)
          loglike_multi[j] <- result_dsr_multi_cand[[j]]$median_loglike
        }
        index_max_like_multi <- which.max(loglike_multi)
        dsr_ca_one_mk <- dsr_ca_one_mk_cand[[index_max_like_multi]]
        tmp <- result_dsr_multi_cand[[index_max_like_multi]]
        log_likes[3] <- tmp$median_loglike
        dsr_multi_median_slope_all_mk[id_mk] <- tmp$median_slope
        dsr_multi_median_intercept_all_mk[id_mk] <- tmp$median_intercept
        dsr_multi_median_sigma_all_mk[id_mk] <- tmp$median_sigma
        dsr_multi_min_us_all_mk[id_mk] <- cand_min_us[index_max_like_multi]
        dsr_multi_rhat_slope_all_mk[id_mk] <- tmp$rhat_slope
        dsr_multi_rhat_intercept_all_mk[id_mk] <- tmp$rhat_intercept
        dsr_multi_rhat_sigma_all_mk[id_mk] <- tmp$rhat_sigma
        dsr_multi_rhat_lp_all_mk[id_mk] <- tmp$rhat_lp
        dsr_multi_trace_all_mk[[id_mk]] <- tmp$traceplot
      }

      # Best model
      index_cand <- which(!is.na(log_likes))
      log_likes <- log_likes[index_cand]
      n_pars <- c(3, 3, 4)
      n_pars <- n_pars[index_cand]
      aics <- mapply(calc_aic, log_likes, n_pars)
      names(aics) <- names(log_likes)
      index_best_model <- which.min(aics)
      if(index_best_model == 1){
        dsr_best_model_all_mk[id_mk] <- "Allele"
        dsr_best_median_slope_all_mk[id_mk] <- dsr_allele_median_slope_all_mk[id_mk]
        dsr_best_median_intercept_all_mk[id_mk] <- dsr_allele_median_intercept_all_mk[id_mk]
        dsr_best_median_sigma_all_mk[id_mk] <- dsr_allele_median_sigma_all_mk[id_mk]
        dsr_best_rhat_slope_all_mk[id_mk] <- dsr_allele_rhat_slope_all_mk[id_mk]
        dsr_best_rhat_intercept_all_mk[id_mk] <- dsr_allele_rhat_intercept_all_mk[id_mk]
        dsr_best_rhat_sigma_all_mk[id_mk] <- dsr_allele_rhat_sigma_all_mk[id_mk]
        dsr_best_rhat_lp_all_mk[id_mk] <- dsr_allele_rhat_lp_all_mk[id_mk]
        dsr_best_trace_all_mk[[id_mk]] <- dsr_allele_trace_all_mk[[id_mk]]
      }else if(index_best_model == 2){
        dsr_best_model_all_mk[id_mk] <- "LUS"
        dsr_best_median_slope_all_mk[id_mk] <- dsr_lus_median_slope_all_mk[id_mk]
        dsr_best_median_intercept_all_mk[id_mk] <- dsr_lus_median_intercept_all_mk[id_mk]
        dsr_best_median_sigma_all_mk[id_mk] <- dsr_lus_median_sigma_all_mk[id_mk]
        dsr_best_rhat_slope_all_mk[id_mk] <- dsr_lus_rhat_slope_all_mk[id_mk]
        dsr_best_rhat_intercept_all_mk[id_mk] <- dsr_lus_rhat_intercept_all_mk[id_mk]
        dsr_best_rhat_sigma_all_mk[id_mk] <- dsr_lus_rhat_sigma_all_mk[id_mk]
        dsr_best_rhat_lp_all_mk[id_mk] <- dsr_lus_rhat_lp_all_mk[id_mk]
        dsr_best_trace_all_mk[[id_mk]] <- dsr_lus_trace_all_mk[[id_mk]]
      }else{
        dsr_best_model_all_mk[id_mk] <- "Multi-sequence"
        dsr_best_median_slope_all_mk[id_mk] <- dsr_multi_median_slope_all_mk[id_mk]
        dsr_best_median_intercept_all_mk[id_mk] <- dsr_multi_median_intercept_all_mk[id_mk]
        dsr_best_median_sigma_all_mk[id_mk] <- dsr_multi_median_sigma_all_mk[id_mk]
        dsr_best_min_us_all_mk[id_mk] <- dsr_multi_min_us_all_mk[id_mk]
        dsr_best_rhat_slope_all_mk[id_mk] <- dsr_multi_rhat_slope_all_mk[id_mk]
        dsr_best_rhat_intercept_all_mk[id_mk] <- dsr_multi_rhat_intercept_all_mk[id_mk]
        dsr_best_rhat_sigma_all_mk[id_mk] <- dsr_multi_rhat_sigma_all_mk[id_mk]
        dsr_best_rhat_lp_all_mk[id_mk] <- dsr_multi_rhat_lp_all_mk[id_mk]
        dsr_best_trace_all_mk[[id_mk]] <- dsr_multi_trace_all_mk[[id_mk]]
      }

      # Assign data
      dsr_data_one_mk <- vector(mode = "list", length = 4)
      names(dsr_data_one_mk) <- c("dsr", "dsr_al", "dsr_lus", "dsr_ca")
      dsr_data_one_mk[[1]] <- dsr_one_mk
      if(consider_dsr["Allele"]){
        dsr_data_one_mk[[2]] <- dsr_al_one_mk
      }
      if(consider_dsr["LUS"]){
        dsr_data_one_mk[[3]] <- dsr_lus_one_mk
      }
      if(consider_dsr["Multi-sequence"]){
        dsr_data_one_mk[[4]] <- dsr_ca_one_mk
      }
      dsr_data_all_mk[[id_mk]] <- dsr_data_one_mk
      aics_dsr_all_mk[[id_mk]] <- aics

      # Allele repeat correction
      info_al_one_mk <- as.numeric(names(info_lus_one_mk))
      if(consider_dsr["Multi-sequence"]){
        info_ca_one_mk <- make_sr_ca_one_mk(cand_min_us[index_max_like_multi], reps_one_mk, seq_al_one_mk, seq_count_one_mk, info_al_one_mk)
      }else{
        info_ca_one_mk <- NULL
      }
      dsr_alcor_all_mk[[id_mk]] <- list(info_al = info_al_one_mk,
                                        info_lus = info_lus_one_mk,
                                        info_ca = info_ca_one_mk)
    }
  }

  ## DSR: multiple loci together
  index_mk_together_dsr <- which(model_dsr_all_mk == "Multiple loci together")
  n_mk_together_dsr <- length(index_mk_together_dsr)
  dsr_together <- dsr_h_together <- id_mk_together <- numeric(0)
  id_mk_count <- 0
  if(n_mk_together_dsr > 0){
    cat("MCMC for the double-back stutter ratio of the multiple loci together\n")

    # Prepare data of the target loci
    for(i in 1:n_mk_together_dsr){
      # Extract data of one marker
      id_mk <- index_mk_together_dsr[i]
      tmp <- extract_sr_use(dsr_data[, id_mk], allele_sr_data[, id_mk], height_sr_data[, id_mk], lus_all_mk[[id_mk]])
      dsr_one_mk <- tmp[[1]]
      dsr_al_one_mk <- tmp[[2]]
      dsr_h_one_mk <- tmp[[3]]

      # Combine data
      dsr_together <- c(dsr_together, dsr_one_mk)
      dsr_h_together <- c(dsr_h_together, dsr_h_one_mk)

      # Make id_mk_together
      id_mk_count <- id_mk_count + 1
      id_mk_together <- c(id_mk_together, rep(id_mk_count, length(dsr_one_mk)))

      # Assign data
      dsr_data_one_mk <- vector(mode = "list", length = 4)
      names(dsr_data_one_mk) <- c("dsr", "dsr_al", "dsr_lus", "dsr_ca")
      dsr_data_one_mk[[1]] <- dsr_one_mk
      dsr_data_one_mk[[2]] <- dsr_al_one_mk
      dsr_data_all_mk[[id_mk]] <- dsr_data_one_mk
    }
    data_list_dsr_together <- list(N = length(dsr_together), NMk = n_mk_together_dsr, SR = dsr_together, Height = dsr_h_together, IndexMk = id_mk_together)
    tmp <- run_mcmc_sr_together(model_sr_together, data_list_dsr_together, n_warmup, n_sampling, n_chains, n_parallel_chains)
    dsr_best_model_all_mk[index_mk_together_dsr] <- "Multiple loci together"
    dsr_best_median_myu_all_mk[index_mk_together_dsr] <- tmp$median_myus
    dsr_best_median_sigma_all_mk[index_mk_together_dsr] <- tmp$median_sigma
    dsr_best_rhat_myu_all_mk[index_mk_together_dsr] <- tmp$rhat_myus
    dsr_best_rhat_sigma_all_mk[index_mk_together_dsr] <- tmp$rhat_sigma
    dsr_best_rhat_lp_all_mk[index_mk_together_dsr] <- tmp$rhat_lp
    dsr_best_trace_all_mk[index_mk_together_dsr] <- tmp$traceplot
  }

  ## DSR: not consider
  index_mk_not_consider_dsr <- which(model_dsr_all_mk == "Not consider")
  if(length(index_mk_not_consider_dsr) > 0){
    dsr_best_model_all_mk[index_mk_not_consider_dsr] <- "Not consider"
  }

  # Perform the MCMC for the minus 2-nt stutter ratio

  ## Define objects
  m2sr_best_model_all_mk <- rep(NA, n_mk)
  m2sr_best_median_myu_all_mk <- m2sr_best_median_sigma_all_mk <- rep(NA, n_mk)
  m2sr_best_rhat_myu_all_mk <- m2sr_best_rhat_sigma_all_mk <- m2sr_best_rhat_lp_all_mk <- rep(NA, n_mk)
  m2sr_best_trace_all_mk <- vector(mode = "list", length = n_mk)

  ## M2SR: locus-specific
  m2sr_data_all_mk <- vector(mode = "list", length = n_mk)
  names(m2sr_data_all_mk) <- markers_par
  index_mk_specific_m2sr <- which(model_m2sr_all_mk == "Uniform")
  n_mk_specific_m2sr <- length(index_mk_specific_m2sr)
  if(n_mk_specific_m2sr > 0){
    for(i in 1:n_mk_specific_m2sr){
      id_mk <- index_mk_specific_m2sr[i]
      cat(paste0("MCMC for the minus 2-nt stutter ratio of the ", markers_par[id_mk], "\n"))

      # Extract data of one marker
      tmp <- extract_sr_use(m2sr_data[, id_mk], allele_sr_data[, id_mk], height_sr_data[, id_mk], lus_all_mk[[id_mk]])
      m2sr_one_mk <- tmp[[1]]
      m2sr_al_one_mk <- tmp[[2]]
      m2sr_h_one_mk <- tmp[[3]]

      data_list_m2sr_uniform <- list(N = length(m2sr_one_mk), SR = m2sr_one_mk, Height = m2sr_h_one_mk)
      tmp <- run_mcmc_sr_uni(model_sr_uni, data_list_m2sr_uniform, n_warmup, n_sampling, n_chains, n_parallel_chains)
      m2sr_best_model_all_mk[id_mk] <- "Uniform"
      m2sr_best_median_myu_all_mk[id_mk] <- tmp$median_myu
      m2sr_best_median_sigma_all_mk[id_mk] <- tmp$median_sigma
      m2sr_best_rhat_myu_all_mk[id_mk] <- tmp$rhat_myu
      m2sr_best_rhat_sigma_all_mk[id_mk] <- tmp$rhat_sigma
      m2sr_best_rhat_lp_all_mk[id_mk] <- tmp$rhat_lp
      m2sr_best_trace_all_mk[[id_mk]] <- tmp$traceplot

      # Assign data
      m2sr_data_one_mk <- vector(mode = "list", length = 4)
      names(m2sr_data_one_mk) <- c("m2sr", "m2sr_al", "m2sr_lus", "m2sr_ca")
      m2sr_data_one_mk[[1]] <- m2sr_one_mk
      m2sr_data_one_mk[[2]] <- m2sr_al_one_mk
      m2sr_data_all_mk[[id_mk]] <- m2sr_data_one_mk
    }
  }

  ## M2SR: multiple loci together
  index_mk_together_m2sr <- which(model_m2sr_all_mk == "Multiple loci together")
  n_mk_together_m2sr <- length(index_mk_together_m2sr)
  m2sr_together <- m2sr_h_together <- id_mk_together <- numeric(0)
  id_mk_count <- 0
  if(n_mk_together_m2sr > 0){
    cat("MCMC for the minus 2-nt stutter ratio of the multiple loci together\n")

    # Prepare data of the target loci
    for(i in 1:n_mk_together_m2sr){
      # Extract data of one marker
      id_mk <- index_mk_together_m2sr[i]
      tmp <- extract_sr_use(m2sr_data[, id_mk], allele_sr_data[, id_mk], height_sr_data[, id_mk], lus_all_mk[[id_mk]])
      m2sr_one_mk <- tmp[[1]]
      m2sr_al_one_mk <- tmp[[2]]
      m2sr_h_one_mk <- tmp[[3]]

      # Combine data
      m2sr_together <- c(m2sr_together, m2sr_one_mk)
      m2sr_h_together <- c(m2sr_h_together, m2sr_h_one_mk)

      # Make id_mk_together
      id_mk_count <- id_mk_count + 1
      id_mk_together <- c(id_mk_together, rep(id_mk_count, length(m2sr_one_mk)))

      # Assign data
      m2sr_data_one_mk <- vector(mode = "list", length = 4)
      names(m2sr_data_one_mk) <- c("m2sr", "m2sr_al", "m2sr_lus", "m2sr_ca")
      m2sr_data_one_mk[[1]] <- m2sr_one_mk
      m2sr_data_one_mk[[2]] <- m2sr_al_one_mk
      m2sr_data_all_mk[[id_mk]] <- m2sr_data_one_mk
    }
    data_list_m2sr_together <- list(N = length(m2sr_together), NMk = n_mk_together_m2sr, SR = m2sr_together, Height = m2sr_h_together, IndexMk = id_mk_together)
    tmp <- run_mcmc_sr_together(model_sr_together, data_list_m2sr_together, n_warmup, n_sampling, n_chains, n_parallel_chains)
    m2sr_best_model_all_mk[index_mk_together_m2sr] <- "Multiple loci together"
    m2sr_best_median_myu_all_mk[index_mk_together_m2sr] <- tmp$median_myus
    m2sr_best_median_sigma_all_mk[index_mk_together_m2sr] <- tmp$median_sigma
    m2sr_best_rhat_myu_all_mk[index_mk_together_m2sr] <- tmp$rhat_myus
    m2sr_best_rhat_sigma_all_mk[index_mk_together_m2sr] <- tmp$rhat_sigma
    m2sr_best_rhat_lp_all_mk[index_mk_together_m2sr] <- tmp$rhat_lp
    m2sr_best_trace_all_mk[index_mk_together_m2sr] <- tmp$traceplot
  }

  ## M2SR: not consider
  index_mk_not_consider_m2sr <- which(model_m2sr_all_mk == "Not consider")
  if(length(index_mk_not_consider_m2sr) > 0){
    m2sr_best_model_all_mk[index_mk_not_consider_m2sr] <- "Not consider"
  }

  #### Gamma parameter ####

  cat("MCMC for the gamma parameter\n")

  # Calculate expected peak heights
  eph_all_data <- calc_eph_all_data(n_oph, tph_all_data, gt_list, peak_list,
                                    ae_median_myu_all_mk,
                                    bsr_alcor_all_mk, bsr_best_model_all_mk, bsr_best_median_myu_all_mk, bsr_best_median_slope_all_mk, bsr_best_median_intercept_all_mk, bsr_best_min_us_all_mk,
                                    fsr_alcor_all_mk, fsr_best_model_all_mk, fsr_best_median_myu_all_mk, fsr_best_median_slope_all_mk, fsr_best_median_intercept_all_mk, fsr_best_min_us_all_mk,
                                    dsr_alcor_all_mk, dsr_best_model_all_mk, dsr_best_median_myu_all_mk, dsr_best_median_slope_all_mk, dsr_best_median_intercept_all_mk, dsr_best_min_us_all_mk,
                                    m2sr_best_model_all_mk, m2sr_best_median_myu_all_mk)

  # Remove data without sequence information (i.e., NA data)
  index_adopt <- which(!is.na(eph_all_data))
  oph_all_data <- oph_all_data[index_adopt]
  eph_all_data <- eph_all_data[index_adopt]

  # Perform MCMC
  data_list_gamma <- list(N = length(oph_all_data), OPH = oph_all_data, EPH = eph_all_data)
  tmp <- run_mcmc_gamma(model_gamma, data_list_gamma, n_warmup, n_sampling, n_chains, n_parallel_chains)
  gamma_median_beta <- tmp$median_beta
  gamma_rhat_beta <- tmp$rhat_beta
  gamma_rhat_lp <- tmp$rhat_lp
  gamma_traceplot <- tmp$traceplot

  return(list("markers_par" = markers_par,
              "model_bsr_all_mk" = model_bsr_all_mk,
              "model_fsr_all_mk" = model_fsr_all_mk,
              "model_dsr_all_mk" = model_dsr_all_mk,
              "ae_data" = ae_data,
              "ae_median_myu_all_mk" = ae_median_myu_all_mk,
              "ae_median_sigma_all_mk" = ae_median_sigma_all_mk,
              "ae_rhat_myu_all_mk" = ae_rhat_myu_all_mk,
              "ae_rhat_sigma_all_mk" = ae_rhat_sigma_all_mk,
              "ae_rhat_lp_all_mk" = ae_rhat_lp_all_mk,
              "ae_trace_all_mk" = ae_trace_all_mk,
              "bsr_data_all_mk" = bsr_data_all_mk,
              "bsr_alcor_all_mk" = bsr_alcor_all_mk,
              "aics_bsr_all_mk" = aics_bsr_all_mk,
              "bsr_allele_median_slope_all_mk" = bsr_allele_median_slope_all_mk,
              "bsr_allele_median_intercept_all_mk" = bsr_allele_median_intercept_all_mk,
              "bsr_allele_median_sigma_all_mk" = bsr_allele_median_sigma_all_mk,
              "bsr_allele_rhat_slope_all_mk" = bsr_allele_rhat_slope_all_mk,
              "bsr_allele_rhat_intercept_all_mk" = bsr_allele_rhat_intercept_all_mk,
              "bsr_allele_rhat_sigma_all_mk" = bsr_allele_rhat_sigma_all_mk,
              "bsr_allele_rhat_lp_all_mk" = bsr_allele_rhat_lp_all_mk,
              "bsr_lus_median_slope_all_mk" = bsr_lus_median_slope_all_mk,
              "bsr_lus_median_intercept_all_mk" = bsr_lus_median_intercept_all_mk,
              "bsr_lus_median_sigma_all_mk" = bsr_lus_median_sigma_all_mk,
              "bsr_lus_rhat_slope_all_mk" = bsr_lus_rhat_slope_all_mk,
              "bsr_lus_rhat_intercept_all_mk" = bsr_lus_rhat_intercept_all_mk,
              "bsr_lus_rhat_sigma_all_mk" = bsr_lus_rhat_sigma_all_mk,
              "bsr_lus_rhat_lp_all_mk" = bsr_lus_rhat_lp_all_mk,
              "bsr_multi_median_slope_all_mk" = bsr_multi_median_slope_all_mk,
              "bsr_multi_median_intercept_all_mk" = bsr_multi_median_intercept_all_mk,
              "bsr_multi_median_sigma_all_mk" = bsr_multi_median_sigma_all_mk,
              "bsr_multi_min_us_all_mk" = bsr_multi_min_us_all_mk,
              "bsr_multi_rhat_slope_all_mk" = bsr_multi_rhat_slope_all_mk,
              "bsr_multi_rhat_intercept_all_mk" = bsr_multi_rhat_intercept_all_mk,
              "bsr_multi_rhat_sigma_all_mk" = bsr_multi_rhat_sigma_all_mk,
              "bsr_multi_rhat_lp_all_mk" = bsr_multi_rhat_lp_all_mk,
              "bsr_best_model_all_mk" = bsr_best_model_all_mk,
              "bsr_best_median_myu_all_mk" = bsr_best_median_myu_all_mk,
              "bsr_best_median_slope_all_mk" = bsr_best_median_slope_all_mk,
              "bsr_best_median_intercept_all_mk" = bsr_best_median_intercept_all_mk,
              "bsr_best_median_sigma_all_mk" = bsr_best_median_sigma_all_mk,
              "bsr_best_min_us_all_mk" = bsr_best_min_us_all_mk,
              "bsr_best_rhat_myu_all_mk" = bsr_best_rhat_myu_all_mk,
              "bsr_best_rhat_slope_all_mk" = bsr_best_rhat_slope_all_mk,
              "bsr_best_rhat_intercept_all_mk" = bsr_best_rhat_intercept_all_mk,
              "bsr_best_rhat_sigma_all_mk" = bsr_best_rhat_sigma_all_mk,
              "bsr_best_rhat_lp_all_mk" = bsr_best_rhat_lp_all_mk,
              "bsr_best_trace_all_mk" = bsr_best_trace_all_mk,
              "fsr_data_all_mk" = fsr_data_all_mk,
              "fsr_alcor_all_mk" = fsr_alcor_all_mk,
              "aics_fsr_all_mk" = aics_fsr_all_mk,
              "fsr_allele_median_slope_all_mk" = fsr_allele_median_slope_all_mk,
              "fsr_allele_median_intercept_all_mk" = fsr_allele_median_intercept_all_mk,
              "fsr_allele_median_sigma_all_mk" = fsr_allele_median_sigma_all_mk,
              "fsr_allele_rhat_slope_all_mk" = fsr_allele_rhat_slope_all_mk,
              "fsr_allele_rhat_intercept_all_mk" = fsr_allele_rhat_intercept_all_mk,
              "fsr_allele_rhat_sigma_all_mk" = fsr_allele_rhat_sigma_all_mk,
              "fsr_allele_rhat_lp_all_mk" = fsr_allele_rhat_lp_all_mk,
              "fsr_lus_median_slope_all_mk" = fsr_lus_median_slope_all_mk,
              "fsr_lus_median_intercept_all_mk" = fsr_lus_median_intercept_all_mk,
              "fsr_lus_median_sigma_all_mk" = fsr_lus_median_sigma_all_mk,
              "fsr_lus_rhat_slope_all_mk" = fsr_lus_rhat_slope_all_mk,
              "fsr_lus_rhat_intercept_all_mk" = fsr_lus_rhat_intercept_all_mk,
              "fsr_lus_rhat_sigma_all_mk" = fsr_lus_rhat_sigma_all_mk,
              "fsr_lus_rhat_lp_all_mk" = fsr_lus_rhat_lp_all_mk,
              "fsr_multi_median_slope_all_mk" = fsr_multi_median_slope_all_mk,
              "fsr_multi_median_intercept_all_mk" = fsr_multi_median_intercept_all_mk,
              "fsr_multi_median_sigma_all_mk" = fsr_multi_median_sigma_all_mk,
              "fsr_multi_min_us_all_mk" = fsr_multi_min_us_all_mk,
              "fsr_multi_rhat_slope_all_mk" = fsr_multi_rhat_slope_all_mk,
              "fsr_multi_rhat_intercept_all_mk" = fsr_multi_rhat_intercept_all_mk,
              "fsr_multi_rhat_sigma_all_mk" = fsr_multi_rhat_sigma_all_mk,
              "fsr_multi_rhat_lp_all_mk" = fsr_multi_rhat_lp_all_mk,
              "fsr_best_model_all_mk" = fsr_best_model_all_mk,
              "fsr_best_median_myu_all_mk" = fsr_best_median_myu_all_mk,
              "fsr_best_median_slope_all_mk" = fsr_best_median_slope_all_mk,
              "fsr_best_median_intercept_all_mk" = fsr_best_median_intercept_all_mk,
              "fsr_best_median_sigma_all_mk" = fsr_best_median_sigma_all_mk,
              "fsr_best_min_us_all_mk" = fsr_best_min_us_all_mk,
              "fsr_best_rhat_myu_all_mk" = fsr_best_rhat_myu_all_mk,
              "fsr_best_rhat_slope_all_mk" = fsr_best_rhat_slope_all_mk,
              "fsr_best_rhat_intercept_all_mk" = fsr_best_rhat_intercept_all_mk,
              "fsr_best_rhat_sigma_all_mk" = fsr_best_rhat_sigma_all_mk,
              "fsr_best_rhat_lp_all_mk" = fsr_best_rhat_lp_all_mk,
              "fsr_best_trace_all_mk" = fsr_best_trace_all_mk,
              "dsr_data_all_mk" = dsr_data_all_mk,
              "dsr_alcor_all_mk" = dsr_alcor_all_mk,
              "aics_dsr_all_mk" = aics_dsr_all_mk,
              "dsr_allele_median_slope_all_mk" = dsr_allele_median_slope_all_mk,
              "dsr_allele_median_intercept_all_mk" = dsr_allele_median_intercept_all_mk,
              "dsr_allele_median_sigma_all_mk" = dsr_allele_median_sigma_all_mk,
              "dsr_allele_rhat_slope_all_mk" = dsr_allele_rhat_slope_all_mk,
              "dsr_allele_rhat_intercept_all_mk" = dsr_allele_rhat_intercept_all_mk,
              "dsr_allele_rhat_sigma_all_mk" = dsr_allele_rhat_sigma_all_mk,
              "dsr_allele_rhat_lp_all_mk" = dsr_allele_rhat_lp_all_mk,
              "dsr_lus_median_slope_all_mk" = dsr_lus_median_slope_all_mk,
              "dsr_lus_median_intercept_all_mk" = dsr_lus_median_intercept_all_mk,
              "dsr_lus_median_sigma_all_mk" = dsr_lus_median_sigma_all_mk,
              "dsr_lus_rhat_slope_all_mk" = dsr_lus_rhat_slope_all_mk,
              "dsr_lus_rhat_intercept_all_mk" = dsr_lus_rhat_intercept_all_mk,
              "dsr_lus_rhat_sigma_all_mk" = dsr_lus_rhat_sigma_all_mk,
              "dsr_lus_rhat_lp_all_mk" = dsr_lus_rhat_lp_all_mk,
              "dsr_multi_median_slope_all_mk" = dsr_multi_median_slope_all_mk,
              "dsr_multi_median_intercept_all_mk" = dsr_multi_median_intercept_all_mk,
              "dsr_multi_median_sigma_all_mk" = dsr_multi_median_sigma_all_mk,
              "dsr_multi_min_us_all_mk" = dsr_multi_min_us_all_mk,
              "dsr_multi_rhat_slope_all_mk" = dsr_multi_rhat_slope_all_mk,
              "dsr_multi_rhat_intercept_all_mk" = dsr_multi_rhat_intercept_all_mk,
              "dsr_multi_rhat_sigma_all_mk" = dsr_multi_rhat_sigma_all_mk,
              "dsr_multi_rhat_lp_all_mk" = dsr_multi_rhat_lp_all_mk,
              "dsr_best_model_all_mk" = dsr_best_model_all_mk,
              "dsr_best_median_myu_all_mk" = dsr_best_median_myu_all_mk,
              "dsr_best_median_slope_all_mk" = dsr_best_median_slope_all_mk,
              "dsr_best_median_intercept_all_mk" = dsr_best_median_intercept_all_mk,
              "dsr_best_median_sigma_all_mk" = dsr_best_median_sigma_all_mk,
              "dsr_best_min_us_all_mk" = dsr_best_min_us_all_mk,
              "dsr_best_rhat_myu_all_mk" = dsr_best_rhat_myu_all_mk,
              "dsr_best_rhat_slope_all_mk" = dsr_best_rhat_slope_all_mk,
              "dsr_best_rhat_intercept_all_mk" = dsr_best_rhat_intercept_all_mk,
              "dsr_best_rhat_sigma_all_mk" = dsr_best_rhat_sigma_all_mk,
              "dsr_best_rhat_lp_all_mk" = dsr_best_rhat_lp_all_mk,
              "dsr_best_trace_all_mk" = dsr_best_trace_all_mk,
              "m2sr_data_all_mk" = m2sr_data_all_mk,
              "m2sr_best_model_all_mk" = m2sr_best_model_all_mk,
              "m2sr_best_median_myu_all_mk" = m2sr_best_median_myu_all_mk,
              "m2sr_best_median_sigma_all_mk" = m2sr_best_median_sigma_all_mk,
              "m2sr_best_rhat_myu_all_mk" = m2sr_best_rhat_myu_all_mk,
              "m2sr_best_rhat_sigma_all_mk" = m2sr_best_rhat_sigma_all_mk,
              "m2sr_best_rhat_lp_all_mk" = m2sr_best_rhat_lp_all_mk,
              "m2sr_best_trace_all_mk" = m2sr_best_trace_all_mk,
              "gamma_median_beta" = gamma_median_beta,
              "gamma_rhat_beta" = gamma_rhat_beta,
              "gamma_rhat_lp" = gamma_rhat_lp,
              "gamma_traceplot" = gamma_traceplot))
}
