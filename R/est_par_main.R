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

  # Trace plot
  fit_ae$draws() %>% mcmc_trace(pars = paste0("alpha[", 1:n_mk, "]")) -> traceplot_ae

  # Extract parameters
  tmp <- as_draws_df(fit_ae$draws())
  lp_values <- tmp$lp__
  density_lp <- density(lp_values)
  lp_mode_estimated <- density_lp$x[which.max(density_lp$y)]
  index_hpd <- which.min(abs(lp_values - lp_mode_estimated))
  lp <- lp_values[index_hpd]

  alpha <- rep(NA, n_mk)
  for(i in 1:n_mk){
    alpha_mk_values <- tmp[[paste0("alpha[", i, "]")]]
    alpha[i] <- alpha_mk_values[index_hpd]
  }

  # Summary data
  fit_ae_summary <- fit_ae$summary()
  variables <- fit_ae_summary$variable
  rhats <- fit_ae_summary$rhat
  rhat_lp <- rhats[which(variables == "lp__")]
  rhat_alpha <- rhats[which(is.element(variables, paste0("alpha[", 1:n_mk, "]")))]

  return(list("lp" = lp,
              "alpha" = alpha,
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
  lp <- tmp$lp
  alpha <- tmp$alpha
  rhat_lp <- tmp$rhat_lp
  rhat_alpha <- tmp$rhat_alpha
  traceplot_ae <- tmp$traceplot_ae

  names(alpha) <- markers_auto
  names(rhat_alpha) <- markers_auto

  return(list("lp" = lp,
              "alpha" = alpha,
              "rhat_lp" = rhat_lp,
              "rhat_alpha" = rhat_alpha,
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
  rph_al_1 <- rph_bs_1 <- rph_fs_1 <- rph_ds_1 <- rph_m2s_1 <- rep(0, length(cand_peaks))
  rph_al_2 <- rph_bs_2 <- rph_fs_2 <- rph_ds_2 <- rph_m2s_2 <- rep(0, length(cand_peaks))

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

  return(rbind(rph_al_1, rph_al_2, rph_bs_1, rph_bs_2, rph_fs_1, rph_fs_2, rph_ds_1, rph_ds_2, rph_m2s_1, rph_m2s_2))
}

#' prepare_gamma_data
#'
#' @description The function to prepare data for estimating the gamma distribution
#' @param markers_auto A vector of the autosomal markers
#' @param peak_list A list of the peaks
#' @param height_list A list of the heights
#' @param gt_list A list of the genotypes
#' @param alpha A vector of the estimated alpha of the Dirichlet distribution for the locus-specific amplification efficiencies
#' @param ave_het_ph_all_data A vector of the average heterozygote allele peak heights
#' @param list_dt_average_sr A list of the average stutter ratios
#' @param sr_consider_auto Information on consideration of back, forward, double-back, and minus 2-nt stutters in autosomal markers
#' @param at_auto A vector of ATs in each autosomal marker
prepare_gamma_data <- function(markers_auto, peak_list, height_list, gt_list, alpha, ave_het_ph_all_data, list_dt_average_sr, sr_consider_auto, at_auto){
  n_sample <- length(gt_list)
  n_mk <- length(markers_auto)

  expect_ae_all_mk <- n_mk * (alpha + 1) / (sum(alpha) + n_mk)

  sum_peaks_all <- length(unlist(peak_list))
  rphs_gamma_list <- vector(mode = "list", length = n_sample)
  tphs_gamma_list <- vector(mode = "list", length = n_sample)
  ophs_gamma_list <- vector(mode = "list", length = n_sample)
  ats_gamma_list <- vector(mode = "list", length = n_sample)
  n_adopt_all <- rep(0, n_sample * n_mk)

  count <- 0
  for(i in 1:n_sample){
    peak_s1 <- peak_list[[i]]
    height_s1 <- height_list[[i]]
    gt_s1 <- gt_list[[i]]
    ave_het_ph_s1 <- ave_het_ph_all_data[i]

    rphs_gamma_list_s1 <- vector(mode = "list", length = n_mk)
    names(rphs_gamma_list_s1) <- markers_auto
    tphs_gamma_list_s1 <- vector(mode = "list", length = n_mk)
    names(tphs_gamma_list_s1) <- markers_auto
    ophs_gamma_list_s1 <- vector(mode = "list", length = n_mk)
    names(ophs_gamma_list_s1) <- markers_auto
    ats_gamma_list_s1 <- vector(mode = "list", length = n_mk)
    names(ats_gamma_list_s1) <- markers_auto

    for(j in 1:n_mk){
      count <- count + 1

      peaks <- peak_s1[[j]]
      peaks <- round(as.numeric(peaks), 1)
      n_peak <- length(peaks)
      if(n_peak != 0){
        gt <- gt_s1[j, ]
        gt <- round(as.numeric(gt), 1)

        heights <- height_s1[[j]]
        ae <- expect_ae_all_mk[j]
        mk <- markers_auto[j]
        dt_average_sr <- list_dt_average_sr[[mk]]
        allele_sr <- round(as.numeric(dt_average_sr$Allele), 1)
        bsr_average <- dt_average_sr[["Average BSR"]]
        fsr_average <- dt_average_sr[["Average FSR"]]
        dsr_average <- dt_average_sr[["Average DSR"]]
        m2sr_average <- dt_average_sr[["Average M2SR"]]
        at <- at_auto[mk]

        cand_peaks <- make_cand_peak(gt, sr_consider_auto[j, ])

        rphs_per_component <- calc_relative_ph(cand_peaks, gt, ae, allele_sr, bsr_average, fsr_average, dsr_average, m2sr_average)
        rphs <- apply(rphs_per_component, 2, sum)

        index_rph <- which(rphs != 0) # If the stutter ratio = 0, rph for the stutter product is zero.
        if(length(index_rph) == 0){
          rphs_gamma_list_s1[[j]] <- NULL
          tphs_gamma_list_s1[[j]] <- numeric(0)
          ophs_gamma_list_s1[[j]] <- numeric(0)
          ats_gamma_list_s1[[j]] <- numeric(0)
        }else{
          rphs_per_component <- rphs_per_component[, index_rph, drop = FALSE]
          rphs <- rphs[index_rph]
          cand_peaks_tmp <- cand_peaks[index_rph]

          index_gt <- which(is.element(cand_peaks_tmp, gt))
          index_height <- which(is.element(cand_peaks_tmp, peaks))
          index_adopt <- sort(unique(c(index_gt, index_height)))
          n_adopt <- length(index_adopt)

          if(n_adopt == 0){
            rphs_gamma_list_s1[[j]] <- NULL
            tphs_gamma_list_s1[[j]] <- numeric(0)
            ophs_gamma_list_s1[[j]] <- numeric(0)
            ats_gamma_list_s1[[j]] <- numeric(0)
          }else{
            if(all(rphs_per_component[, index_adopt, drop = FALSE] == 0)){
              print(paste0("Sample: ", i, ", Marker: ", j))
            }

            rphs_gamma_list_s1[[j]] <- rphs_per_component[, index_adopt, drop = FALSE]
            tphs_gamma_list_s1[[j]] <- rep(ave_het_ph_s1, n_adopt)
            cand_peaks_adopt <- cand_peaks_tmp[index_adopt]

            ophs <- rep(0, n_adopt)
            index_oph <- which(is.element(peaks, cand_peaks_adopt))
            if(length(index_oph) != 0){
              peaks_adopt <- peaks[index_oph]
              heights_adopt <- heights[index_oph]
              ophs[match(peaks_adopt, cand_peaks_adopt)] <- heights_adopt
            }
            ophs_gamma_list_s1[[j]] <- ophs

            ats_gamma_list_s1[[j]] <- rep(at, n_adopt)

            n_adopt_all[count] <- n_adopt
          }
        }
      }else{
        rphs_gamma_list_s1[[j]] <- NULL
        tphs_gamma_list_s1[[j]] <- numeric(0)
        ophs_gamma_list_s1[[j]] <- numeric(0)
        ats_gamma_list_s1[[j]] <- numeric(0)
      }
    }
    rphs_gamma_list[[i]] <- rphs_gamma_list_s1
    tphs_gamma_list[[i]] <- tphs_gamma_list_s1
    ophs_gamma_list[[i]] <- ophs_gamma_list_s1
    ats_gamma_list[[i]] <- ats_gamma_list_s1
  }

  # Bind rphs
  rphs_gamma <- matrix(NA, 10, sum(n_adopt_all))
  index_start <- 1
  for(i in 1:n_sample){
    for(j in 1:n_mk){
      n_adopt <- n_adopt_all[n_mk * (i - 1) + j]
      if(n_adopt != 0){
        rphs_gamma[, index_start:(index_start + n_adopt - 1)] <- rphs_gamma_list[[i]][[j]]
        index_start <- index_start + n_adopt
      }
    }
  }

  tphs_gamma <- unlist(tphs_gamma_list)
  ophs_gamma <- unlist(ophs_gamma_list)
  ats_gamma <- unlist(ats_gamma_list)

  return(list("rphs_gamma" = rphs_gamma,
              "tphs_gamma" = tphs_gamma,
              "ophs_gamma" = ophs_gamma,
              "ats_gamma" = ats_gamma))
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
  data_list_gamma <- list(N = length(tphs_gamma), RPH = rphs_gamma, TPH = tphs_gamma, OPH = ophs_gamma, AT = ats_gamma)

  fit_gamma <- model_gamma$sample(
    data_list_gamma,
    iter_warmup = warmup,
    iter_sampling = sampling,
    chains = chains,
    parallel_chains = n_parallel_chains
  )

  # Trace plot
  fit_gamma$draws() %>% mcmc_trace(pars = c("cv_allele")) -> traceplot_cv_allele
  fit_gamma$draws() %>% mcmc_trace(pars = c("cv_stutter")) -> traceplot_cv_stutter

  # Extract parameters
  tmp <- as_draws_df(fit_gamma$draws())
  lp_values <- tmp$lp__
  density_lp <- density(lp_values)
  lp_mode_estimated <- density_lp$x[which.max(density_lp$y)]
  index_hpd <- which.min(abs(lp_values - lp_mode_estimated))
  lp <- lp_values[index_hpd]

  cv_allele_values <- tmp$cv_allele
  cv_allele <- cv_allele_values[index_hpd]

  cv_stutter_values <- tmp$cv_stutter
  cv_stutter <- cv_stutter_values[index_hpd]

  fit_gamma_summary <- fit_gamma$summary()
  variables <- fit_gamma_summary$variable
  rhats <- fit_gamma_summary$rhat
  rhat_lp <- rhats[which(variables == "lp__")]
  rhat_cv_allele <- rhats[which(variables == "cv_allele")]
  rhat_cv_stutter <- rhats[which(variables == "cv_stutter")]

  return(list("lp" = lp,
              "cv_allele" = cv_allele,
              "cv_stutter" = cv_stutter,
              "rhat_lp" = rhat_lp,
              "rhat_cv_allele" = rhat_cv_allele,
              "rhat_cv_stutter" = rhat_cv_stutter,
              "traceplot_cv_allele" = traceplot_cv_allele,
              "traceplot_cv_stutter" = traceplot_cv_stutter))
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

  # Perform MCMC for AE
  tmp <- perform_mcmc_ae_all_mk(prop_ae_data_use, markers_auto,
                                model_ae,
                                warmup, sampling, chains, n_parallel_chains)
  lp_ae <- tmp$lp
  alpha <- tmp$alpha
  rhat_lp_ae <- tmp$rhat_lp
  rhat_alpha <- tmp$rhat_alpha
  traceplot_ae <- tmp$traceplot_ae

  # Prepare data for average heterozygote allele peak heights
  ave_het_ph_all_data <- calc_ave_het_ph_all_data(height_list)

  # Prepare peak height data for estimating gamma distribution
  tmp <- prepare_gamma_data(markers_auto, peak_list, height_list, gt_list, alpha, ave_het_ph_all_data, list_dt_average_sr, sr_consider_auto, at_auto)
  rphs_gamma <- tmp$rphs_gamma
  tphs_gamma <- tmp$tphs_gamma
  ophs_gamma <- tmp$ophs_gamma
  ats_gamma <- tmp$ats_gamma

  # Perform MCMC for gamma distribution
  tmp <- perform_mcmc_gamma(rphs_gamma, tphs_gamma, ophs_gamma, ats_gamma,
                            model_gamma,
                            warmup, sampling, chains, n_parallel_chains)
  lp_gamma <- tmp$lp
  cv_allele <- tmp$cv_allele
  cv_stutter <- tmp$cv_stutter
  rhat_lp_gamma <- tmp$rhat_lp
  rhat_cv_allele <- tmp$rhat_cv_allele
  rhat_cv_stutter <- tmp$rhat_cv_stutter
  traceplot_cv_allele <- tmp$traceplot_cv_allele
  traceplot_cv_stutter <- tmp$traceplot_cv_stutter

  return(list("lp_ae" = lp_ae,
              "alpha" = alpha,
              "rhat_lp_ae" = rhat_lp_ae,
              "rhat_alpha" = rhat_alpha,
              "traceplot_ae" = traceplot_ae,
              "lp_gamma" = lp_gamma,
              "cv_allele" = cv_allele,
              "cv_stutter" = cv_stutter,
              "rhat_lp_gamma" = rhat_lp_gamma,
              "rhat_cv_allele" = rhat_cv_allele,
              "rhat_cv_stutter" = rhat_cv_stutter,
              "traceplot_cv_allele" = traceplot_cv_allele,
              "traceplot_cv_stutter" = traceplot_cv_stutter))
}
