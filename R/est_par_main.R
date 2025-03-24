#' create_dt_ae_data
#'
#' @description The function to create a data.table of the locus-specific amplification efficiencies
#' @param peak_list A list of peaks
#' @param height_list A list of heights
#' @param gt_list A list of genotypes
create_dt_ae_data <- function(peak_list, height_list, gt_list){
  sample_names <- names(peak_list)
  n_sample <- length(sample_names)
  markers <- names(peak_list[[1]])
  n_mk <- length(markers)

  n_data <- n_sample * n_mk
  sample_name_all <- rep("", n_data)
  marker_all <- rep("", n_data)
  temp_amount_all <- rep(NA, n_data)
  prop_ae_all <- rep(NA, n_data)
  ae_all <- rep(NA, n_data)
  bool_ae_use_all <- rep(FALSE, n_data)

  current_index <- 1
  for(i in 1:n_sample){
    peak_list_one <- peak_list[[i]]
    height_list_one <- height_list[[i]]
    gt_list_one <- gt_list[[i]]
    bool_detect <- TRUE
    for(j in 1:n_mk){
      peaks <- peak_list_one[[j]]
      heights <- height_list_one[[j]]
      gt <- unique(gt_list_one[j, ])
      if(!all(is.element(gt, peaks))){
        bool_detect <- FALSE
        break
      }
    }

    index_assign <- current_index:(current_index + n_mk - 1)
    sample_name_all[index_assign] <- sample_names[i]
    marker_all[index_assign] <- markers

    if(bool_detect){
      temp_amount <- sapply(height_list_one, sum)
      sum_temp_amount <- sum(temp_amount)
      temp_amount_all[index_assign] <- temp_amount
      prop_ae_all[index_assign] <- temp_amount / sum_temp_amount
      ae_all[index_assign] <- n_mk * temp_amount / sum_temp_amount
      bool_ae_use_all[index_assign] <- TRUE
    }

    current_index <- current_index + n_mk
  }

  dt_ae_data <- data.table(SampleName = sample_name_all,
                           Marker = marker_all,
                           TempAmount = temp_amount_all,
                           PropAE = prop_ae_all,
                           AE = ae_all,
                           Bool_Use = bool_ae_use_all)

  return(dt_ae_data)
}

#' calc_ave_het_ph_all_data
#'
#' @description The function to calculate the average heterozygote allele peak heights of all marker
#' @param height_list The list of peak heights
calc_ave_het_ph_all_data <- function(height_list){
  n_data <- length(height_list)
  ave_het_ph_all_data <- rep(0, n_data)
  for(i in 1:n_data){
    height_list_one_data <- height_list[[i]]
    ave_het_ph_per_mk <- sapply(height_list_one_data, sum)
    ave_het_ph_all_data[i] <- mean(ave_het_ph_per_mk) / 2
  }
  return(ave_het_ph_all_data)
}

#' run_mcmc_ae
#'
#' @description The function to perform MCMC for the locus-specific amplification efficiency
#' @param model_ae The stan model of the locus-specific amplification efficiency
#' @param data_list_ae The data list of the locus-specific amplification efficiency
#' @param n_mk The number of markers
#' @param warmup The number of warm-up iterations for each chain
#' @param sampling The number of sampling iterations for each chain
#' @param chains The number of chains
#' @param n_parallel_chains The number of parallel chains
run_mcmc_ae <- function(model_ae, data_list_ae, n_mk, warmup, sampling, chains, n_parallel_chains){
  # Perform MCMC
  fit_ae <- model_ae$sample(
    data_list_ae,
    iter_warmup = warmup,
    iter_sampling = sampling,
    chains = chains,
    parallel_chains = n_parallel_chains
  )

  # Summary data
  fit_ae_summary <- fit_ae$summary()
  variables <- fit_ae_summary$variable
  medians <- fit_ae_summary$median
  rhats <- fit_ae_summary$rhat
  median_lp <- medians[which(variables == "lp__")]
  median_alpha <- medians[which(is.element(variables, paste0("alpha[", 1:n_mk, "]")))]
  rhat_lp <- rhats[which(variables == "lp__")]
  rhat_alpha <- rhats[which(is.element(variables, paste0("alpha[", 1:n_mk, "]")))]

  # Trace plot
  fit_ae$draws() %>% mcmc_trace(pars = paste0("alpha[", 1:n_mk, "]")) -> traceplot_ae

  return(list("median_lp" = median_lp,
              "median_alpha" = median_alpha,
              "rhat_lp" = rhat_lp,
              "rhat_alpha" = rhat_alpha,
              "traceplot_ae" = traceplot_ae))
}

#' perform_mcmc_ae_all_mk
#'
#' @description The function to perform MCMC for amplification efficiencies of all loci
#' @param prop_ae_data_use A matrix of the proportion of the amplification efficiencies (only data used for MCMC)
#' @param markers_auto A vector of the autosomal markers
#' @param model_ae The stan model of the locus-specific amplification efficiency
#' @param warmup The number of warm-up iterations for each chain
#' @param sampling The number of sampling iterations for each chain
#' @param chains The number of chains
#' @param n_parallel_chains The number of parallel chains
perform_mcmc_ae_all_mk <- function(prop_ae_data_use, markers_auto,
                                   model_ae,
                                   warmup, sampling, chains, n_parallel_chains){
  n_mk <- length(markers_auto)
  data_list_ae <- list(N = ncol(prop_ae_data_use), N_Mk = nrow(prop_ae_data_use), PAE = prop_ae_data_use)

  tmp <- run_mcmc_ae(model_ae, data_list_ae, n_mk, warmup, sampling, chains, n_parallel_chains)
  median_alpha <- tmp$median_alpha
  median_lp <- tmp$median_lp
  rhat_alpha <- tmp$rhat_alpha
  rhat_lp <- tmp$rhat_lp
  traceplot_ae <- tmp$traceplot_ae

  names(median_alpha) <- markers_auto
  names(rhat_alpha) <- markers_auto

  return(list("median_alpha" = median_alpha,
              "median_lp" = median_lp,
              "rhat_alpha" = rhat_alpha,
              "rhat_lp" = rhat_lp,
              "traceplot_ae" = traceplot_ae))
}

#' calc_relative_ph
#'
#' @description The function to calculate the relative peak heights in one locus
#' @param peaks A vector of the peaks in one locus
#' @param gt A vector of the genotype in one locus
#' @param ae The amplification efficiency in one locus
#' @param allele_sr A vector of the alleles for stutter ratios
#' @param bsr_average A vector of the average back stutter ratios
#' @param fsr_average A vector of the average forward stutter ratios
#' @param dsr_average A vector of the average double-back stutter ratios
#' @param m2sr_average A vector of the average minus 2-nt stutter ratios
calc_relative_ph <- function(peaks, gt, ae,
                             allele_sr, bsr_average, fsr_average, dsr_average, m2sr_average){
  peaks <- round(as.numeric(peaks), 1)
  gt <- round(as.numeric(gt), 1)
  allele_sr <- round(as.numeric(allele_sr), 1)

  rph_al_1 <- rph_bs_1 <- rph_fs_1 <- rph_ds_1 <- rph_m2s_1 <- rep(0, length(peaks))
  rph_al_2 <- rph_bs_2 <- rph_fs_2 <- rph_ds_2 <- rph_m2s_2 <- rep(0, length(peaks))

  # Heterozygote
  if(gt[1] != gt[2]){
    for(i in 1:2){
      al <- gt[i]
      index_sr_al <- which(is.element(allele_sr, al))
      if(length(index_sr_al) == 1){
        bsr <- bsr_average[index_sr_al]
        fsr <- fsr_average[index_sr_al]
        dsr <- dsr_average[index_sr_al]
        m2sr <- m2sr_average[index_sr_al]
      }else{
        bsr <- 0
        fsr <- 0
        dsr <- 0
        m2sr <- 0
      }

      index_al <- which(is.element(peaks, al))
      bs <- calc_bs_one(al)
      index_bs <- which(is.element(peaks, bs))
      fs <- calc_fs_one(al)
      index_fs <- which(is.element(peaks, fs))
      ds <- calc_ds_one(al)
      index_ds <- which(is.element(peaks, ds))
      m2s <- calc_m2s_one(al)
      index_m2s <- which(is.element(peaks, m2s))

      if(i == 1){
        rph_al_1[index_al] <- ae / (1 + bsr + fsr + dsr + m2sr)
        if(length(index_bs) == 1){
          rph_bs_1[index_bs] <- ae * bsr / (1 + bsr + fsr + dsr + m2sr)
        }
        if(length(index_fs) == 1){
          rph_fs_1[index_fs] <- ae * fsr / (1 + bsr + fsr + dsr + m2sr)
        }
        if(length(index_ds) == 1){
          rph_ds_1[index_ds] <- ae * dsr / (1 + bsr + fsr + dsr + m2sr)
        }
        if(length(index_m2s) == 1){
          rph_m2s_1[index_m2s] <- ae * m2sr / (1 + bsr + fsr + dsr + m2sr)
        }
      }else{
        rph_al_2[index_al] <- ae / (1 + bsr + fsr + dsr + m2sr)
        if(length(index_bs) == 1){
          rph_bs_2[index_bs] <- ae * bsr / (1 + bsr + fsr + dsr + m2sr)
        }
        if(length(index_fs) == 1){
          rph_fs_2[index_fs] <- ae * fsr / (1 + bsr + fsr + dsr + m2sr)
        }
        if(length(index_ds) == 1){
          rph_ds_2[index_ds] <- ae * dsr / (1 + bsr + fsr + dsr + m2sr)
        }
        if(length(index_m2s) == 1){
          rph_m2s_2[index_m2s] <- ae * m2sr / (1 + bsr + fsr + dsr + m2sr)
        }
      }
    }
  }

  return(cbind(rph_al_1, rph_al_2, rph_bs_1, rph_bs_2, rph_fs_1, rph_fs_2, rph_ds_1, rph_ds_2, rph_m2s_1, rph_m2s_2))
}

#' prepare_gamma_data
#'
#' @description The function to prepare data for estimating the gamma distribution
#' @param markers_auto A vector of the autosomal markers
#' @param peak_use A list of the peaks used for estimating the gamma distribution
#' @param gt_use A list of the genotypes used for estimating the gamma distribution
#' @param ae_data_use A matrix of the actual locus-specific amplification efficiencies in each sample
#' @param ave_het_ph_use A vector of the average heterozygote allele peak heights used for estimating the gamma distribution
#' @param list_dt_average_sr A list of the average stutter ratios
prepare_gamma_data <- function(markers_auto, peak_use, gt_use, ae_data_use, ave_het_ph_use, list_dt_average_sr){
  n_gt <- length(gt_use)
  n_mk <- length(markers_auto)

  expect_ae_all_mk <- n_mk * (est_alpha + 1) / (sum(est_alpha) + n_mk)

  sum_peaks_all <- length(unlist(peak_use))
  rphs_gamma <- matrix(-1, sum_peaks_all, 10)
  tphs_gamma <- rep(-1, sum_peaks_all)

  index_start <- 1
  for(i in 1:n_gt){
    peak_s1 <- peak_use[[i]]
    gt_s1 <- gt_use[[i]]
    ave_het_ph_s1 <- ave_het_ph_use[i]

    ae_all_mk <- ae_data_use[, i]

    for(j in 1:n_mk){
      peaks <- peak_s1[[j]]
      n_peak <- length(peaks)
      if(n_peak != 0){
        ae <- ae_all_mk[j]

        mk <- markers_auto[j]
        dt_average_sr <- list_dt_average_sr[[mk]]
        allele_sr <- dt_average_sr$Allele
        bsr_average <- dt_average_sr[["Average BSR"]]
        fsr_average <- dt_average_sr[["Average FSR"]]
        dsr_average <- dt_average_sr[["Average DSR"]]
        m2sr_average <- dt_average_sr[["Average M2SR"]]

        rphs_gamma[index_start:(index_start + n_peak - 1), ] <- calc_relative_ph(peaks, gt_s1[j, ], ae,
                                                                                 allele_sr, bsr_average, fsr_average, dsr_average, m2sr_average)
        tphs_gamma[index_start:(index_start + n_peak - 1)] <- ave_het_ph_s1
        index_start <- index_start + n_peak
      }
    }
  }

  return(list("rphs_gamma" = rphs_gamma,
              "tphs_gamma" = tphs_gamma))
}

#' perform_mcmc_gamma
#'
#' @description The function to perform MCMC for the gamma distribution
#' @param rphs_gamma The relative peak heights for the gamma distribution
#' @param ophs_gamma The observed peak heights for the gamma distribution
#' @param tphs_gamma The total peak heights in one heterozygote allele for the gamma distribution
#' @param model_gamma The stan model of the gamma distribution
#' @param warmup The number of warm-up iterations for each chain
#' @param sampling The number of sampling iterations for each chain
#' @param chains The number of chains
#' @param n_parallel_chains The number of parallel chains
perform_mcmc_gamma <- function(rphs_gamma, ophs_gamma, tphs_gamma,
                               model_gamma,
                               warmup, sampling, chains, n_parallel_chains){
  data_list_gamma <- list(N = length(rphs_gamma), OPH = ophs_gamma, RPH = rphs_gamma, TPH = tphs_gamma)

  fit_gamma <- model_gamma$sample(
    data_list_gamma,
    iter_warmup = warmup,
    iter_sampling = sampling,
    chains = chains,
    parallel_chains = n_parallel_chains
  )

  # Summary data
  fit_gamma_summary <- fit_gamma$summary()

  # Median of the parameter 'cv'
  median_lp <- fit_gamma_summary$median[which(fit_gamma_summary$variable == "lp__")]
  median_cv <- fit_gamma_summary$median[which(fit_gamma_summary$variable == "cv")]

  # rhat of the parameter 'cv' and the likelihood
  rhat_lp <- fit_gamma_summary$rhat[which(fit_gamma_summary$variable == "lp__")]
  rhat_cv <- fit_gamma_summary$rhat[which(fit_gamma_summary$variable == "cv")]

  # Trace plot
  fit_gamma$draws() %>% mcmc_trace(pars = c("cv")) -> traceplot_cv

  return(list("median_lp" = median_lp,
              "median_cv" = median_cv,
              "rhat_lp" = rhat_lp,
              "rhat_cv" = rhat_cv,
              "traceplot_cv" = traceplot_cv))
}

#' estimate_parameters
#'
#' @description The function to estimate parameters of the locus-specific amplification efficiencies and the gamma distribution
estimate_parameters <- function(dt_exp, dt_gt, list_dt_average_sr, dt_setting_est_par,
                                model_ae, model_gamma,
                                warmup, sampling, chains, n_parallel_chains = detectCores()){
  # Extract information on setting
  markers_auto <- dt_setting_est_par$Marker
  n_mk_auto <- length(markers_auto)
  at_auto <- dt_setting_est_par$AT
  names(at_auto) <- markers_auto

  # Information on considering stutter ratios
  sr_consider_auto <- matrix(TRUE, n_mk_auto, 4)
  rownames(sr_consider_auto) <- markers_auto
  colnames(sr_consider_auto) <- c("BS", "FS", "DS", "M2S")
  for(i in 1:n_mk_auto){
    mk <- markers_auto[i]
    dt_average_sr <- list_dt_average_sr[[mk]]
    sr_consider_auto[i, "BS"] <- any(dt_average_sr[["Average BSR"]] != 0)
    sr_consider_auto[i, "FS"] <- any(dt_average_sr[["Average FSR"]] != 0)
    sr_consider_auto[i, "DS"] <- any(dt_average_sr[["Average DSR"]] != 0)
    sr_consider_auto[i, "M2S"] <- any(dt_average_sr[["Average M2SR"]] != 0)
  }

  # Arrange experimental data
  tmp <- arrange_exp_data(dt_exp, dt_gt,
                          markers_auto, sr_consider_auto, at_auto)
  peak_list <- tmp$peak_list
  height_list <- tmp$height_list
  gt_list <- tmp$gt_list

  # Prepare data for AE
  dt_ae_data <- create_dt_ae_data(peak_list, height_list, gt_list)
  prop_ae_data <- matrix(dt_ae_data$PropAE, nrow = n_mk_auto)
  index_prop_ae_use <- which(apply(prop_ae_data, 2, prod) > 0)
  prop_ae_data_use <- prop_ae_data[, index_prop_ae_use, drop = FALSE]

  # Prepare data for average heterozygote allele peak heights
  ave_het_ph_all_data <- calc_ave_het_ph_all_data(height_list)
  ave_het_ph_use <- ave_het_ph_all_data[index_prop_ae_use]

  # Perform MCMC for AE
  tmp <- perform_mcmc_ae_all_mk(prop_ae_data_use, markers_auto,
                                model_ae,
                                warmup, sampling, chains, n_parallel_chains)
  median_alpha <- tmp$median_alpha
  rhat_alpha <- tmp$rhat_alpha
  median_lp_ae <- tmp$median_lp
  rhat_lp_ae <- tmp$rhat_lp
  traceplot_ae <- tmp$traceplot_ae

  # Prepare peak height data for estimating gamma distribution
  tmp <- prepare_gamma_data(markers_auto, peak_list[index_prop_ae_use], gt_list[index_prop_ae_use], prop_ae_data_use * n_mk_auto, ave_het_ph_use, list_dt_average_sr)
  rphs_gamma <- tmp$rphs_gamma
  rphs_gamma <- apply(rphs_gamma, 1, sum)
  tphs_gamma <- tmp$tphs_gamma
  ophs_gamma <- unlist(height_list[index_prop_ae_use])

  index_rm <- which(rphs_gamma == 0)
  if(length(index_rm) != 0){
    rphs_gamma <- rphs_gamma[-index_rm]
    tphs_gamma <- tphs_gamma[-index_rm]
    ophs_gamma <- ophs_gamma[-index_rm]
  }

  # Remove data with oph > 30000 RFU
  index_rm <- which(ophs_gamma > 30000)
  if(length(index_rm) != 0){
    rphs_gamma <- rphs_gamma[-index_rm]
    tphs_gamma <- tphs_gamma[-index_rm]
    ophs_gamma <- ophs_gamma[-index_rm]
  }

  # Perform MCMC for gamma distribution
  tmp <- perform_mcmc_gamma(rphs_gamma, ophs_gamma, tphs_gamma,
                            model_gamma,
                            warmup, sampling, chains, n_parallel_chains)
  median_lp_gamma <- tmp$median_lp
  median_cv <- tmp$median_cv
  rhat_lp_gamma <- tmp$rhat_lp
  rhat_cv <- tmp$rhat_cv
  traceplot_cv <- tmp$traceplot

  return(list("median_alpha" = median_alpha,
              "rhat_alpha" = rhat_alpha,
              "median_lp_ae" = median_lp_ae,
              "rhat_lp_ae" = rhat_lp_ae,
              "traceplot_ae" = traceplot_ae,
              "median_lp_gamma" = median_lp_gamma,
              "median_cv" = median_cv,
              "rhat_lp_gamma" = rhat_lp_gamma,
              "rhat_cv" = rhat_cv,
              "traceplot_cv" = traceplot_cv))
}
