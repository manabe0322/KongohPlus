#' create_dt_ae_data
#'
#' @description The function to create a data.table of the locus-specific amplification efficiencies
#' @param height_list A list of heights
create_dt_ae_data <- function(height_list){
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
    height_list_one <- height_list[[i]]
    temp_amount <- sapply(height_list_one, sum)

    index_assign <- current_index:(current_index + n_mk - 1)
    sample_name_all[index_assign] <- sample_names[i]
    marker_all[index_assign] <- markers
    temp_amount_all[index_assign] <- temp_amount
    prop_ae_all[index_assign] <- temp_amount / sum(temp_amount)
    ae_all[index_assign] <- n_mk * temp_amount / sum(temp_amount)
    if(!is.element(0, temp_amount)){
      bool_ae_use_all[index_assign] <- TRUE
    }
  }

  dt_ae_data <- data.table(SampleName = sample_name_all,
                           Marker = marker_all,
                           TempAmount = temp_amount_all,
                           PropAE = prop_ae_all,
                           AE = ae_all,
                           Bool_Use = bool_ae_use_all)

  return(dt_ae_data)
}

#' run_mcmc_ae
#'
#' @description The function to perform the MCMC for the locus-specific amplification efficiency
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

perform_mcmc_ae_all_mk <- function(pae_data_use, markers_auto,
                                   model_ae,
                                   warmup, sampling, chains, n_parallel_chains){
  n_mk <- length(markers_auto)
  data_list_ae <- list(N = ncol(pae_data_use), N_Mk = nrow(pae_data_use), PAE = pae_data_use)

  tmp <- run_mcmc_ae(model_ae, data_list_ae, n_mk, warmup, sampling, chains, n_parallel_chains)
  median_lp <- tmp$median_lp
  median_alpha <- tmp$median_alpha
  rhat_alpha <- tmp$rhat_alpha
  rhat_lp <- tmp$rhat_lp
  traceplot_ae <- tmp$traceplot_ae

  dt_result_est_ae <- data.table(Marker = markers_auto,
                                 Alpha = median_alpha,
                                 Rhat_Alpha = rhat_alpha)

  return(list("dt_result_est_ae" = dt_result_est_ae,
              "traceplot_ae" = traceplot_ae,
              "median_lp" = median_lp,
              "rhat_lp" = rhat_lp))
}

tmp_func2 <- function(peaks, gt){
  peaks <- as.numeric(peaks)
  gt <- as.numeric(gt)

  rph_base <- ae_actual

  rph_al_1 <- rph_bs_1 <- rph_fs_1 <- rph_ds_1 <- rph_m2s_1 <- rep(0, length(peaks))
  rph_al_2 <- rph_bs_2 <- rph_fs_2 <- rph_ds_2 <- rph_m2s_2 <- rep(0, length(peaks))

  # Heterozygote
  if(gt[1] != gt[2]){
    for(i in 1:length(gt)){
      al <- gt[i]

      index_al <- which(is.element(peaks, al))

      bs <- calc_bs_one(al)
      index_bs <- which(is.element(peaks, bs))

      fs <- calc_fs_one(al)
      index_fs <- which(is.element(peaks, fs))

      ds <- calc_ds_one(al)
      index_ds <- which(is.element(peaks, ds))

      m2s <- calc_m2s_one(al)
      index_m2s <- which(is.element(peaks, m2s))


      index_sr_al <- which(is.element(allele_sr, al))

      if(length(index_al) == 1){
        if(sr_consider["BS"]){
          bs <- calc_bs_one(al)
          index_bs <- which(is.element(peaks, bs))
          if(length(index_sr_al) == 1){
            bsr <- bsr_average[index_sr_al]
          }else{
            bsr <- 0
          }
        }else{
          bsr <- 0
        }

        if(sr_consider["FS"]){
          fs <- calc_fs_one(al)
          index_fs <- which(is.element(peaks, fs))
          if(length(index_fsr) == 1){
            fsr <- fsr_average[index_sr_al]
          }else{
            fsr <- 0
          }
        }else{
          fsr <- 0
        }

        if(sr_consider["DS"]){
          ds <- calc_ds_one(al)
          index_ds <- which(is.element(peaks, ds))
          if(length(index_sr_al) == 1){
            dsr <- dsr_average[index_sr_al]
          }else{
            dsr <- 0
          }
        }else{
          dsr <- 0
        }
      }
    }
  }

}

tmp_func <- function(peak_use, gt_use, prop_ae_data_use, list_dt_average_sr){
  n_gt <- length(gt_use)
  n_mk <- length(markers_auto)

  for(i in 1:n_gt){
    peak_s1 <- peak_use[[i]]
    gt_s1 <- gt_use[[i]]
    pae_s1 <- prop_ae_data_use[, i]

    for(j in 1:n_mk){
      peaks <- peak_s1[[j]]
      n_peak <- length(peaks)
      if(n_peak != 0){
        ae_actual <- n_mk * pae_s1[j]

        mk <- markers_auto[j]
        dt_average_sr <- list_dt_average_sr[[mk]]
        allele_sr <- as.numeric(dt_average_sr$Allele)
        bsr_average <- dt_average_sr[["Average BSR"]]
        fsr_average <- dt_average_sr[["Average FSR"]]
        dsr_average <- dt_average_sr[["Average DSR"]]
        m2sr_average <- dt_average_sr[["Average M2SR"]]
      }
    }
  }
}

estimate_parameters <- function(dt_exp, dt_gt, list_dt_average_sr, dt_kit, dt_setting_est_par, kit_name,
                                model_ae, model_gamma,
                                warmup, sampling, chains, pu_remove, pu_dist, at_sex,
                                n_parallel_chains = detectCores()){
  # Extract information on the kit
  markers_all <- unique(dt_kit$Marker)
  n_mk_all <- length(markers_all)
  dye_all_mk <- extract_dye(dt_kit, markers_all)

  # Extract information on setting
  markers_auto <- dt_setting_est_par$Marker
  n_mk_auto <- length(markers_auto)
  at_auto <- dt_setting_est_par$AT
  names(at_auto) <- markers_auto

  # Make an object for ATs of all markers
  at_all_mk <- rep(NA, n_mk_all)
  names(at_all_mk) <- markers_all
  at_all_mk[match(markers_auto, markers_all)] <- at_auto
  at_all_mk[is.na(at_all_mk)] <- at_sex

  # Information on considering stutter ratios
  sr_consider_auto <- matrix(TRUE, n_mk_auto, 4)
  rownames(sr_consider_auto) <- markers_auto
  colnames(sr_consider_auto) <- c("BS", "FS", "DS", "M2S")
  for(i in 1:n_mk_auto){
    mk <- markers_auto[i]
    dt_average_sr <- list_dt_average_sr[[mk]]
    sr_consider_auto[i, "BS"] <- any(!is.na(dt_average_sr[["Average BSR"]]))
    sr_consider_auto[i, "FS"] <- any(!is.na(dt_average_sr[["Average FSR"]]))
    sr_consider_auto[i, "DS"] <- any(!is.na(dt_average_sr[["Average DSR"]]))
    sr_consider_auto[i, "M2S"] <- any(!is.na(dt_average_sr[["Average M2SR"]]))
  }

  # Arrange experimental data
  tmp <- arrange_exp_data(dt_exp, dt_gt,
                          markers_auto, sr_consider_auto,
                          markers_all, at_all_mk, dye_all_mk,
                          pu_remove, pu_dist)
  peak_list <- tmp$peak_list
  size_list <- tmp$size_list
  height_list <- tmp$height_list
  gt_list <- tmp$gt_list

  # Prepare data for AE
  dt_ae_data <- create_dt_ae_data(height_list)
  prop_ae_data <- matrix(dt_ae_data$PropAE, nrow = n_mk_auto)
  index_prop_ae_use <- which(apply(prop_ae_data, 2, prod) > 0)
  prop_ae_data_use <- prop_ae_data[, index_prop_ae_use, drop = FALSE]

  #tph_all_data <- calc_tph_all_data(height_list)
  #tph_data_use <- tph_all_data[index_pae_use]

  # Perform MCMC for AE
  tmp <- perform_mcmc_ae_all_mk(dt_ae_data, markers_auto,
                                model_ae,
                                warmup, sampling, chains, n_parallel_chains)
  dt_result_est_ae <- tmp$dt_result_est_ae
  traceplot_ae <- tmp$traceplot_ae
  median_lp <- tmp$median_lp
  rhat_lp <- tmp$rhat_lp
}
