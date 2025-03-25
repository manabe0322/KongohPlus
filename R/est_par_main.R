#' create_dt_ae_data
#'
#' @description The function to create a data.table of the locus-specific amplification efficiencies
#' @param peak_list A list of peaks
#' @param height_list A list of heights
create_dt_ae_data <- function(peak_list, height_list){
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

  current_index <- 1
  for(i in 1:n_sample){
    peak_list_one <- peak_list[[i]]
    height_list_one <- height_list[[i]]

    index_assign <- current_index:(current_index + n_mk - 1)
    sample_name_all[index_assign] <- sample_names[i]
    marker_all[index_assign] <- markers

    temp_amount <- sapply(height_list_one, sum)
    sum_temp_amount <- sum(temp_amount)
    temp_amount_all[index_assign] <- temp_amount
    prop_ae_all[index_assign] <- temp_amount / sum_temp_amount
    ae_all[index_assign] <- n_mk * temp_amount / sum_temp_amount

    current_index <- current_index + n_mk
  }

  dt_ae_data <- data.table(SampleName = sample_name_all,
                           Marker = marker_all,
                           TempAmount = temp_amount_all,
                           PropAE = prop_ae_all,
                           AE = ae_all)

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
#' @param cand_peaks A vector of the candidate detectable peaks
#' @param gt A vector of the genotype in one locus
#' @param ae The amplification efficiency in one locus
#' @param allele_sr A vector of the alleles for stutter ratios
#' @param bsr_average A vector of the average back stutter ratios
#' @param fsr_average A vector of the average forward stutter ratios
#' @param dsr_average A vector of the average double-back stutter ratios
#' @param m2sr_average A vector of the average minus 2-nt stutter ratios
calc_relative_ph <- function(cand_peaks, gt, ae, allele_sr, bsr_average, fsr_average, dsr_average, m2sr_average){
  allele_sr <- round(as.numeric(allele_sr), 1)

  rph_al_1 <- rph_bs_1 <- rph_fs_1 <- rph_ds_1 <- rph_m2s_1 <- rep(0, length(cand_peaks))
  rph_al_2 <- rph_bs_2 <- rph_fs_2 <- rph_ds_2 <- rph_m2s_2 <- rep(0, length(cand_peaks))

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

      index_al <- which(is.element(cand_peaks, al))
      bs <- calc_bs_one(al)
      index_bs <- which(is.element(cand_peaks, bs))
      fs <- calc_fs_one(al)
      index_fs <- which(is.element(cand_peaks, fs))
      ds <- calc_ds_one(al)
      index_ds <- which(is.element(cand_peaks, ds))
      m2s <- calc_m2s_one(al)
      index_m2s <- which(is.element(cand_peaks, m2s))

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
#' @param peak_list A list of the peaks
#' @param height_list A list of the heights
#' @param gt_list A list of the genotypes
#' @param median_alpha A vector of the estimated alpha of the Dirichlet distribution for the locus-specific amplification efficiencies
#' @param ave_het_ph_all_data A vector of the average heterozygote allele peak heights
#' @param list_dt_average_sr A list of the average stutter ratios
#' @param sr_consider_auto Information on consideration of back, forward, double-back, and minus 2-nt stutters in autosomal markers
prepare_gamma_data <- function(markers_auto, peak_list, height_list, gt_list, median_alpha, ave_het_ph_all_data, list_dt_average_sr, sr_consider_auto){
  n_sample <- length(gt_list)
  n_mk <- length(markers_auto)

  expect_ae_all_mk <- n_mk * (median_alpha + 1) / (sum(median_alpha) + n_mk)

  sum_peaks_all <- length(unlist(peak_list))
  rphs_gamma <- matrix(-1, n_sample, 2 * n_mk)
  tphs_gamma <- matrix(-1, n_sample, 2 * n_mk)
  ophs_gamma <- matrix(-1, n_sample, 2 * n_mk)

  index_start <- 1
  for(i in 1:n_sample){
    peak_s1 <- peak_list[[i]]
    height_s1 <- height_list[[i]]
    gt_s1 <- gt_list[[i]]
    ave_het_ph_s1 <- ave_het_ph_all_data[i]

    for(j in 1:n_mk){
      peaks <- peak_s1[[j]]
      n_peak <- length(peaks)
      if(n_peak != 0){
        gt <- gt_s1[j, ]
        gt <- round(as.numeric(gt), 1)

        # Only heterozygote
        if(gt[1] != gt[2]){
          heights <- height_s1[[j]]
          ae <- expect_ae_all_mk[j]
          mk <- markers_auto[j]
          dt_average_sr <- list_dt_average_sr[[mk]]
          allele_sr <- dt_average_sr$Allele
          bsr_average <- dt_average_sr[["Average BSR"]]
          fsr_average <- dt_average_sr[["Average FSR"]]
          dsr_average <- dt_average_sr[["Average DSR"]]
          m2sr_average <- dt_average_sr[["Average M2SR"]]

          cand_peaks <- make_cand_peak(gt, sr_consider_auto[j, ])

          rph_per_component <- calc_relative_ph(cand_peaks, gt, ae, allele_sr, bsr_average, fsr_average, dsr_average, m2sr_average)
          rphs <- apply(rph_per_component, 1, sum)

          index_al <- match(gt, cand_peaks)
          index_al <- index_al[!is.na(index_al)]

          rphs_gamma[i, c(2 * j - 1, 2 * j)] <- rphs[index_al]
          tphs_gamma[i, c(2 * j - 1, 2 * j)] <- ave_het_ph_s1

          index_oph <- match(gt, round(as.numeric(peaks), 1))
          for(k in 1:2){
            if(is.na(index_oph[k])){
              ophs_gamma[i, 2 * j - 2 + k] <- 0
            }else{
              ophs_gamma[i, 2 * j - 2 + k] <- heights[index_oph[k]]
            }
          }
        }
      }
    }
  }

  return(list("rphs_gamma" = rphs_gamma,
              "tphs_gamma" = tphs_gamma,
              "ophs_gamma" = ophs_gamma))
}

#' perform_mcmc_gamma
#'
#' @description The function to perform MCMC for the gamma distribution
#' @param rphs_gamma The relative peak heights for the gamma distribution
#' @param tphs_gamma The total peak heights in one heterozygote allele for the gamma distribution
#' @param ophs_gamma The observed peak heights for the gamma distribution
#' @param ats_gamma The analytical thresholds for the gamma distribution
#' @param model_gamma The stan model of the gamma distribution
#' @param warmup The number of warm-up iterations for each chain
#' @param sampling The number of sampling iterations for each chain
#' @param chains The number of chains
#' @param n_parallel_chains The number of parallel chains
perform_mcmc_gamma <- function(rphs_gamma, tphs_gamma, ophs_gamma, ats_gamma,
                               model_gamma,
                               warmup, sampling, chains, n_parallel_chains){
  data_list_gamma <- list(N = length(rphs_gamma), OPH = ophs_gamma, RPH = rphs_gamma, TPH = tphs_gamma, AT = ats_gamma)

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
  dt_ae_data <- create_dt_ae_data(peak_list, height_list)
  prop_ae_data <- matrix(dt_ae_data$PropAE, nrow = n_mk_auto)
  prop_ae_data_use <- prop_ae_data[, which(apply(prop_ae_data, 2, prod) > 0), drop = FALSE]

  # Prepare data for average heterozygote allele peak heights
  ave_het_ph_all_data <- calc_ave_het_ph_all_data(height_list)

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
  tmp <- prepare_gamma_data(markers_auto, peak_list, height_list, gt_list, median_alpha, ave_het_ph_all_data, list_dt_average_sr, sr_consider_auto)
  rphs_gamma <- tmp$rphs_gamma
  tphs_gamma <- tmp$tphs_gamma
  ophs_gamma <- tmp$ophs_gamma

  ats_gamma <- rep(as.vector(sapply(at_auto, rep, 2)), length(peak_list))

  index_rm <- which(rphs_gamma == -1)
  if(length(index_rm) != 0){
    rphs_gamma <- rphs_gamma[-index_rm]
    tphs_gamma <- tphs_gamma[-index_rm]
    ophs_gamma <- ophs_gamma[-index_rm]
    ats_gamma <- ats_gamma[-index_rm]
  }

  # Perform MCMC for gamma distribution
  tmp <- perform_mcmc_gamma(rphs_gamma, tphs_gamma, ophs_gamma, ats_gamma,
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
