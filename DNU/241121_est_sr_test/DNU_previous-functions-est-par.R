#' calc_aic
#'
#' @description The function to calculate the AIC
#' @param log_like A logarithmic likelihood
#' @param n_par The number of parameters
calc_aic <- function(log_like, n_par){
  return(-2 * log_like + 2 * n_par)
}

#' calc_tph_all_data
#'
#' @description The function to calculate the average heterozygote peak heights of all marker
#' @param height_list The list of peak heights
calc_tph_all_data <- function(height_list){
  n_data <- length(height_list)
  tph_all_data <- rep(0, n_data)
  for(i in 1:n_data){
    height_list_one_data <- height_list[[i]]
    tph_per_mk <- sapply(height_list_one_data, sum)
    tph_all_data[i] <- mean(tph_per_mk) / 2
  }
  return(tph_all_data)
}

#' make_dt_par
#'
#' @description The function to make a data.table for estimated parameters
make_dt_par <- function(dt_result_ae, dt_result_bsr, dt_result_fsr, dt_result_dsr, dt_result_m2sr, dt_result_gamma, kit_par){
  markers <- dt_result_ae[, Marker]
  n_mk <- length(markers)
  dt_par <- data.table(Kit = rep(kit_par, 5 * n_mk + 1),
                       Marker = c(rep(markers, 5), NA),
                       Parameter = c(rep("AE", n_mk), rep("BSR", n_mk), rep("FSR", n_mk), rep("DSR", n_mk), rep("M2SR", n_mk), "Gamma"),
                       Model = c(rep("Dirichlet", n_mk), dt_result_bsr[, Model], dt_result_fsr[, Model], dt_result_dsr[, Model], dt_result_m2sr[, Model], "Gamma"),
                       Alpha = c(dt_result_ae[, Alpha], rep(NA, 4 * n_mk + 1)),
                       Myu = c(rep(NA, n_mk), dt_result_bsr[, Myu], dt_result_fsr[, Myu], dt_result_dsr[, Myu], dt_result_m2sr[, Myu], NA),
                       Slope = c(rep(NA, n_mk), dt_result_bsr[, Slope], dt_result_fsr[, Slope], dt_result_dsr[, Slope], rep(NA, n_mk), NA),
                       Intercept = c(rep(NA, n_mk), dt_result_bsr[, Intercept], dt_result_fsr[, Intercept], dt_result_dsr[, Intercept], rep(NA, n_mk), NA),
                       Sigma = c(rep(NA, n_mk), dt_result_bsr[, Sigma], dt_result_fsr[, Sigma], dt_result_dsr[, Sigma], dt_result_m2sr[, Sigma], NA),
                       MinUS = c(rep(NA, n_mk), dt_result_bsr[, MinUS], dt_result_fsr[, MinUS], dt_result_dsr[, MinUS], rep(NA, n_mk), NA),
                       Cv = c(rep(NA, 5 * n_mk), dt_result_gamma[, Cv]),
                       Rhat_Alpha = c(dt_result_ae[, Rhat_Alpha], rep(NA, 4 * n_mk + 1)),
                       Rhat_Myu = c(rep(NA, n_mk), dt_result_bsr[, Rhat_Myu], dt_result_fsr[, Rhat_Myu], dt_result_dsr[, Rhat_Myu], dt_result_m2sr[, Rhat_Myu], NA),
                       Rhat_Slope = c(rep(NA, n_mk), dt_result_bsr[, Rhat_Slope], dt_result_fsr[, Rhat_Slope], dt_result_dsr[, Rhat_Slope], rep(NA, n_mk), NA),
                       Rhat_Intercept = c(rep(NA, n_mk), dt_result_bsr[, Rhat_Intercept], dt_result_fsr[, Rhat_Intercept], dt_result_dsr[, Rhat_Intercept], rep(NA, n_mk), NA),
                       Rhat_Sigma = c(rep(NA, n_mk), dt_result_bsr[, Rhat_Sigma], dt_result_fsr[, Rhat_Sigma], dt_result_dsr[, Rhat_Sigma], dt_result_m2sr[, Rhat_Sigma], NA),
                       Rhat_Cv = c(rep(NA, 5 * n_mk), dt_result_gamma[, Rhat_Cv]),
                       Rhat_Lp = c(dt_result_ae[, Rhat_Lp], dt_result_bsr[, Rhat_Lp], dt_result_fsr[, Rhat_Lp], dt_result_dsr[, Rhat_Lp], dt_result_m2sr[, Rhat_Lp], dt_result_gamma[, Rhat_Lp]))
  return(dt_par)
}

#' estimate_parameters
#'
#' @description The function to estimate parameters
#' @param dt_exp The data.table of the experimental data
#' @param dt_gt The data.table of the genotype data
#' @param dt_seq The data.table of the sequence data
#' @param dt_kit The data.table of information on the kit
#' @param dt_kit_set The data.table of the kit-setting
#' @param kit_par Kit name
#' @param model_ae The stan model for the locus-specific amplification efficiency
#' @param model_sr_gen The stan model for the general stutter ratio
#' @param model_sr_uni The stan model for the uniform stutter ratio
#' @param model_gamma The stan model for gamma distribution
#' @param warmup The number of warm-up iterations for each chain
#' @param sampling The number of sampling iterations for each chain
#' @param chains The number of chains
#' @param pu_remove Whether pull-up peaks are automatically removed or not
#' @param pu_dist Maximum pull-up distance (bp)
#' @param at_sex The analytical threshold of the sex-chromosomal markers
#' @param min_n_motif_lus The minimum number of repeat motifs for considering the LUS model
#' @param min_n_motif_multi The minimum number of repeat motifs for considering the multi-sequence model
#' @param n_parallel_chains The number of parallel chains
estimate_parameters <- function(dt_exp, dt_gt, dt_seq, dt_kit, dt_kit_set, kit_par,
                                model_ae, model_sr_gen, model_sr_uni, model_gamma,
                                warmup, sampling, chains,
                                pu_remove, pu_dist, at_sex, min_n_motif_lus, min_n_motif_multi, max_n_motif_multi = 10,
                                n_parallel_chains = detectCores()){
  # Extract information on the kit
  markers_all <- unique(dt_kit[, Marker])
  n_mk_all <- length(markers_all)
  dye_all_mk <- extract_dye(dt_kit, markers_all)

  # Extract information on setting
  markers_auto <- dt_kit_set[, Marker]
  n_mk_auto <- length(markers_auto)
  at_auto <- dt_kit_set[, AT]
  names(at_auto) <- markers_auto
  meth_model_bsr_auto <- dt_kit_set[, Model_BSR]
  meth_model_fsr_auto <- dt_kit_set[, Model_FSR]
  meth_model_dsr_auto <- dt_kit_set[, Model_DSR]
  meth_model_m2sr_auto <- dt_kit_set[, Model_M2SR]

  # Make an object for ATs of all markers
  at_all_mk <- rep(NA, n_mk_all)
  names(at_all_mk) <- markers_all
  at_all_mk[match(markers_auto, markers_all)] <- at_auto
  at_all_mk[is.na(at_all_mk)] <- at_sex

  # Information on considering stutter ratios
  sr_consider_auto <- matrix(TRUE, n_mk_auto, 4)
  rownames(sr_consider_auto) <- markers_auto
  colnames(sr_consider_auto) <- c("BS", "FS", "DS", "M2S")
  sr_consider_auto[, 1] <- meth_model_bsr_auto != "Not consider"
  sr_consider_auto[, 2] <- meth_model_fsr_auto != "Not consider"
  sr_consider_auto[, 3] <- meth_model_dsr_auto != "Not consider"
  sr_consider_auto[, 4] <- meth_model_m2sr_auto != "Not consider"

  # Arrange experimental data
  tmp <- arrange_exp_data(dt_exp, dt_gt,
                          markers_auto, sr_consider_auto,
                          markers_all, at_all_mk, dye_all_mk,
                          pu_remove, pu_dist)
  peak_list <- tmp$peak_list
  size_list <- tmp$size_list
  height_list <- tmp$height_list
  gt_list <- tmp$gt_list

  # Prepare data for the locus-specific amplification efficiency
  pae_data <- get_pae_data(peak_list, height_list, gt_list)
  index_pae_use <- which(apply(pae_data, 2, prod) > 0)
  pae_data_use <- pae_data[, index_pae_use, drop = FALSE]
  tph_all_data <- calc_tph_all_data(height_list)
  tph_data_use <- tph_all_data[index_pae_use]

  # Prepare data for the stutter ratios
  tmp <- get_sr_data(peak_list, height_list, gt_list, sr_consider_auto)
  bsr_data <- tmp$bsr_data
  fsr_data <- tmp$fsr_data
  dsr_data <- tmp$dsr_data
  m2sr_data <- tmp$m2sr_data
  parent_al_data <- tmp$parent_al_data
  parent_ph_data <- tmp$parent_ph_data
  total_al_prod_data <- tmp$total_al_prod_data

  # Prepare data related to the sequence data
  tmp <- analyze_seq(dt_seq, markers_auto)
  lus_auto <- tmp$lus_auto
  reps_auto <- tmp$reps_auto
  seq_al_auto <- tmp$seq_al_auto
  seq_count_auto <- tmp$seq_count_auto

  # Perform the MCMC for the locus-specific amplification efficiency
  ae_result_list <- analyze_ae_main(markers_auto, model_ae, pae_data_use,
                                    warmup, sampling, chains, n_parallel_chains)

  # Perform the MCMC for the back stutter ratio
  bsr_result_list <- analyze_bfdstr_main("BSR", markers_auto, meth_model_bsr_auto, model_sr_gen, model_sr_uni,
                                         bsr_data, parent_al_data, total_al_prod_data,
                                         lus_auto, reps_auto, seq_al_auto, seq_count_auto,
                                         warmup, sampling, chains, n_parallel_chains, min_n_motif_lus, min_n_motif_multi, max_n_motif_multi)

  # Perform the MCMC for the forward stutter ratio
  fsr_result_list <- analyze_bfdstr_main("FSR", markers_auto, meth_model_fsr_auto, model_sr_gen, model_sr_uni,
                                         fsr_data, parent_al_data, total_al_prod_data,
                                         lus_auto, reps_auto, seq_al_auto, seq_count_auto,
                                         warmup, sampling, chains, n_parallel_chains, min_n_motif_lus, min_n_motif_multi, max_n_motif_multi)

  # Perform the MCMC for the double-back stutter ratio
  dsr_result_list <- analyze_bfdstr_main("DSR", markers_auto, meth_model_dsr_auto, model_sr_gen, model_sr_uni,
                                         dsr_data, parent_al_data, total_al_prod_data,
                                         lus_auto, reps_auto, seq_al_auto, seq_count_auto,
                                         warmup, sampling, chains, n_parallel_chains, min_n_motif_lus, min_n_motif_multi, max_n_motif_multi)

  # Perform the MCMC for the minus 2-nt stutter ratio
  m2sr_result_list <- analyze_m2sr_main(markers_auto, meth_model_m2sr_auto, model_sr_uni,
                                        m2sr_data, parent_al_data, total_al_prod_data,
                                        lus_auto,
                                        warmup, sampling, chains, n_parallel_chains)

  dt_result_ae <- ae_result_list$dt_result_ae
  dt_result_bsr <- bsr_result_list$dt_result
  dt_result_fsr <- fsr_result_list$dt_result
  dt_result_dsr <- dsr_result_list$dt_result
  dt_result_m2sr <- m2sr_result_list$dt_result
  dt_alcor_bsr <- bsr_result_list$dt_alcor
  dt_alcor_fsr <- fsr_result_list$dt_alcor
  dt_alcor_dsr <- dsr_result_list$dt_alcor

  tmp <- make_data_gamma(markers_auto, gt_list[index_pae_use], peak_list[index_pae_use], pae_data_use, tph_data_use, sr_consider_auto, dt_alcor_bsr, dt_alcor_fsr, dt_alcor_dsr, dt_result_m2sr)
  rphs_gamma <- unlist(tmp$rphs_gamma)
  tphs_gamma <- unlist(tmp$tphs_gamma)
  ophs_gamma <- unlist(height_list[index_pae_use])

  #remove NA
  index_not_na <- which(!is.na(rphs_gamma))
  rphs_gamma <- rphs_gamma[index_not_na]
  tphs_gamma <- tphs_gamma[index_not_na]
  ophs_gamma <- ophs_gamma[index_not_na]

  data_list_gamma <- list(N = length(ophs_gamma), OPH = ophs_gamma, RPH = rphs_gamma, TPH = tphs_gamma)
  tmp <- run_mcmc_gamma(model_gamma, data_list_gamma, warmup, sampling, chains, n_parallel_chains)
  dt_result_gamma <- data.table(Cv = tmp$median_cv,
                                Rhat_Cv = tmp$rhat_cv,
                                Rhat_Lp = tmp$rhat_lp)
  gamma_result_list <- list(dt_result = dt_result_gamma,
                            traceplot = tmp$traceplot)

  #### Data.table for stutter ratios ####

  dt_alcor_combined <- make_dt_alcor_combined(dt_alcor_bsr, dt_alcor_fsr, dt_alcor_dsr)
  dt_sr <- make_dt_sr(dt_alcor_combined, dt_result_m2sr)

  #### Results of estimated parameters ####

  dt_par <- make_dt_par(dt_result_ae, dt_result_bsr, dt_result_fsr, dt_result_dsr, dt_result_m2sr, dt_result_gamma, kit_par)

  return(list("markers_auto" = markers_auto,
              "ae_result_list" = ae_result_list,
              "bsr_result_list" = bsr_result_list,
              "fsr_result_list" = fsr_result_list,
              "dsr_result_list" = dsr_result_list,
              "m2sr_result_list" = m2sr_result_list,
              "gamma_result_list" = gamma_result_list,
              "dt_sr" = dt_sr,
              "dt_par" = dt_par))
}

#' analyze_ae_main
#'
#' @description The function to analyze AE data
#' @param markers_auto The vector of autosomal markers
#' @param model_ae The stan model for AE
#' @param pae_data_use The data of the proportion of the AE which is used for MCMC
#' @param warmup The number of warm-up iterations for each chain
#' @param sampling The number of sampling iterations for each chain
#' @param chains The number of chains
#' @param n_parallel_chains The number of parallel chains
analyze_ae_main <- function(markers_auto, model_ae, pae_data_use,
                            warmup, sampling, chains, n_parallel_chains){

  cat("MCMC for the locus-specific amplification efficiency\n")
  n_mk <- length(markers_auto)
  data_list_ae <- list(N = ncol(pae_data_use), N_Mk = nrow(pae_data_use), PAE = pae_data_use)
  tmp <- run_mcmc_ae(model_ae, data_list_ae, n_mk, warmup, sampling, chains, n_parallel_chains)
  median_alpha <- tmp$median_alpha
  rhat_alpha <- tmp$rhat_alpha
  rhat_lp <- tmp$rhat_lp
  traceplot_ae <- tmp$traceplot_ae

  dt_result_ae <- data.table(Marker = markers_auto,
                             Alpha = median_alpha,
                             Rhat_Alpha = rhat_alpha,
                             Rhat_Lp = c(rhat_lp, rep(NA, n_mk - 1)))

  return(list("dt_result_ae" = dt_result_ae,
              "traceplot_ae" = traceplot_ae))
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

  return(list("median_alpha" = median_alpha,
              "rhat_lp" = rhat_lp,
              "rhat_alpha" = rhat_alpha,
              "traceplot_ae" = traceplot_ae))
}

#' calc_sr
#'
#' @description The function to calculate stutter ratios
#' @param peaks A character vector of peaks
#' @param heights A numeric vector of heights
#' @param gt A character vector of the genotype in one marker
calc_sr <- function(peaks, heights, gt, sr_consider){

  # bsr: back stutter ratio
  # fsr: forward stutter ratio
  # dsr: double-back stutter ratio
  # m2sr: minus 2-nt stutter ratio
  # parent_al: parent allele data for stutter ratios
  # parent_ph: peak height data for parent alleles
  # total_al_prod: total allelic product data for stutter ratios
  bsr <- fsr <- dsr <- m2sr <- parent_al <- parent_ph <- total_al_prod <- rep(0, 2)

  peaks <- as.numeric(peaks)

  gt <- as.numeric(gt)
  gt <- sort(unique(gt[!is.na(gt)]))

  ##############
  # Homozygote #
  ##############

  if(length(gt) == 1){
    al_1 <- gt

    # No allelic drop-out
    if(is.element(al_1, peaks)){
      parent_al[1] <- al_1
      parent_ph[1] <- h_al <- heights[peaks == al_1]
      total_al_prod[1] <- sum(heights)

      # Back stutter
      if(sr_consider["BS"]){
        bs <- calc_bs_one(al_1)
        if(is.element(bs, peaks)){
          bsr[1] <- heights[peaks == bs] / h_al
        }
      }

      # Forward stutter
      if(sr_consider["FS"]){
        fs <- calc_fs_one(al_1)
        if(is.element(fs, peaks)){
          fsr[1] <- heights[peaks == fs] / h_al
        }
      }

      # Double-back stutter
      if(sr_consider["DS"]){
        ds <- calc_ds_one(al_1)
        if(is.element(ds, peaks)){
          dsr[1] <- heights[peaks == ds] / h_al
        }
      }

      # Minus 2-nt stutter
      if(sr_consider["M2S"]){
        m2s <- calc_m2s_one(al_1)
        if(is.element(m2s, peaks)){
          m2sr[1] <- heights[peaks == m2s] / h_al
        }
      }
    }

    ################
    # Heterozygote #
    ################

  }else{

    # Other than one repeat difference
    if(abs(round(gt[1] - gt[2], 1)) != 1){

      for(i in 1:2){
        al_1 <- gt[i]

        # No allelic drop-out
        if(is.element(al_1, peaks)){
          parent_al[i] <- al_1
          h_sum <- h_al <- heights[peaks == al_1]
          parent_ph[i] <- h_al

          # Candidates of shared peaks
          al_2 <- gt[3 - i]
          cand_share <- al_2
          if(sr_consider["BS"]){
            cand_share <- c(cand_share, calc_bs_one(al_2))
          }
          if(sr_consider["FS"]){
            cand_share <- c(cand_share, calc_fs_one(al_2))
          }
          if(sr_consider["DS"]){
            cand_share <- c(cand_share, calc_ds_one(al_2))
          }
          if(sr_consider["M2S"]){
            cand_share <- c(cand_share, calc_m2s_one(al_2))
          }
          cand_share <- round(sort(unique(cand_share)), 1)

          # Back stutter
          if(sr_consider["BS"]){
            bs <- calc_bs_one(al_1)
            if(is.element(bs, peaks) && !is.element(bs, cand_share)){
              h_bs <- heights[peaks == bs]
              bsr[i] <- h_bs / h_al
              h_sum <- h_sum + h_bs
            }
          }

          # Forward stutter
          if(sr_consider["FS"]){
            fs <- calc_fs_one(al_1)
            if(is.element(fs, peaks) && !is.element(fs, cand_share)){
              h_fs <- heights[peaks == fs]
              fsr[i] <- h_fs / h_al
              h_sum <- h_sum + h_fs
            }
          }

          # Double-back stutter
          if(sr_consider["DS"]){
            ds <- calc_ds_one(al_1)
            if(is.element(ds, peaks) && !is.element(ds, cand_share)){
              h_ds <- heights[peaks == ds]
              dsr[i] <- h_ds / h_al
              h_sum <- h_sum + h_ds
            }
          }

          # Minus 2-nt stutter
          if(sr_consider["M2S"]){
            m2s <- calc_m2s_one(al_1)
            if(is.element(m2s, peaks) && !is.element(m2s, cand_share)){
              h_m2s <- heights[peaks == m2s]
              m2sr[i] <- h_m2s / h_al
              h_sum <- h_sum + h_m2s
            }
          }
          total_al_prod[i] <- h_sum
        }
      }
    }
  }
  return(list(bsr, fsr, dsr, m2sr, parent_al, parent_ph, total_al_prod))
}

#' get_sr_data
#'
#' @description The function to get data of stutter ratios
#' @param peak_list A list of peaks
#' @param height_list A list of heights
#' @param sr_consider_all_mk Information on consideration of back, forward, double-back, and minus 2-nt stutters in all markers
get_sr_data <- function(peak_list, height_list, gt_list, sr_consider_all_mk){
  n_sample <- length(peak_list)
  markers <- names(peak_list[[1]])
  n_mk <- length(markers)

  bsr_data <- fsr_data <- dsr_data <- m2sr_data <- parent_al_data <- parent_ph_data <- total_al_prod_data <- matrix(0, 2 * n_sample, n_mk)

  # Define row and column names
  tmp <- sapply(names(height_list), rep, 2)
  tmp[1, ] <- paste0(tmp[1, ], "_data1")
  tmp[2, ] <- paste0(tmp[2, ], "_data2")
  rownames(bsr_data) <- rownames(fsr_data) <- rownames(dsr_data) <- rownames(m2sr_data) <- rownames(parent_al_data) <- rownames(parent_ph_data) <- rownames(total_al_prod_data) <- as.vector(tmp)
  colnames(bsr_data) <- colnames(fsr_data) <- colnames(dsr_data) <- colnames(m2sr_data) <- colnames(parent_al_data) <- colnames(parent_ph_data) <- colnames(total_al_prod_data) <- markers

  for(i in 1:n_sample){
    peak_list_one <- peak_list[[i]]
    height_list_one <- height_list[[i]]
    gt_list_one <- gt_list[[i]]

    for(j in 1:n_mk){
      peaks <- peak_list_one[[j]]
      heights <- height_list_one[[j]]
      gt <- gt_list_one[j, ]

      if(length(peaks) > 0){
        sr_consider <- sr_consider_all_mk[j, ]

        sr_one_mk <- calc_sr(peaks, heights, gt, sr_consider)
        bsr_data[c(2 * i - 1, 2 * i), j] <- sr_one_mk[[1]]
        fsr_data[c(2 * i - 1, 2 * i), j] <- sr_one_mk[[2]]
        dsr_data[c(2 * i - 1, 2 * i), j] <- sr_one_mk[[3]]
        m2sr_data[c(2 * i - 1, 2 * i), j] <- sr_one_mk[[4]]
        parent_al_data[c(2 * i - 1, 2 * i), j] <- sr_one_mk[[5]]
        parent_ph_data[c(2 * i - 1, 2 * i), j] <- sr_one_mk[[6]]
        total_al_prod_data[c(2 * i - 1, 2 * i), j] <- sr_one_mk[[7]]
      }
    }
  }

  return(list("bsr_data" = bsr_data,
              "fsr_data" = fsr_data,
              "dsr_data" = dsr_data,
              "m2sr_data" = m2sr_data,
              "parent_al_data" = parent_al_data,
              "parent_ph_data" = parent_ph_data,
              "total_al_prod_data" = total_al_prod_data))
}

#' analyze_seq
#'
#' @description The function to analyze sequence data
#' @param dt_seq The data.table of sequence data
#' @param markers_auto The vector of autosomal markers
analyze_seq <- function(dt_seq, markers_auto){
  lus_auto <- list()
  reps_auto <- list()
  seq_al_auto <- list()
  seq_count_auto <- list()

  for(i in 1:length(markers_auto)){
    # Extract data of one marker
    dt_seq_one_mk <- dt_seq[Marker == markers_auto[i]]
    alleles <- as.numeric(dt_seq_one_mk[, Allele])
    al_names <- sort(unique(alleles))
    counts <- dt_seq_one_mk[, Count]
    sequences <- dt_seq_one_mk[, RepeatRegion]

    # Investigate the number of repeats in each motif
    reps_one_mk <- lapply(sequences, pick_rep)

    # Calculate weighted LUS values
    lus_all_reps <- sapply(reps_one_mk, max)
    lus_one_mk <- rep(0, length(al_names))
    for(j in 1:length(al_names)){
      index_al <- which(alleles == al_names[j])
      lus_one_mk[j] <- calc_weighted_lus(lus_all_reps[index_al], counts[index_al])
    }
    names(lus_one_mk) <- al_names

    # Assign data of one marker
    lus_auto[[i]] <- lus_one_mk
    reps_auto[[i]] <- reps_one_mk
    seq_al_auto[[i]] <- alleles
    seq_count_auto[[i]] <- counts
  }

  names(lus_auto) <- markers_auto
  names(reps_auto) <- markers_auto
  names(seq_al_auto) <- markers_auto
  names(seq_count_auto) <- markers_auto

  return(list("lus_auto" = lus_auto,
              "reps_auto" = reps_auto,
              "seq_al_auto" = seq_al_auto,
              "seq_count_auto" = seq_count_auto))
}

#' make_sr_ca_one_mk
#'
#' @description The function to make data of corrected allele numbers in one marker
#' @param min_us Minimum uninterrupted stretch
#' @param reps_one_mk Repeat motifs of one marker
#' @param seq_al_one_mk Alleles for sequence data in one marker
#' @param seq_count_one_mk Counts for sequence data in one marker
#' @param sr_al_one_mk A vector of the parent alleles of each stutter ratio
make_sr_ca_one_mk <- function(min_us, reps_one_mk, seq_al_one_mk, seq_count_one_mk, sr_al_one_mk){
  sr_ca_one_mk <- rep(0, length(sr_al_one_mk))
  al_names <- sort(unique(seq_al_one_mk))
  for(i in 1:length(al_names)){
    index_seq_al <- which(seq_al_one_mk == al_names[i])
    sr_ca_one_mk[which(sr_al_one_mk == al_names[i])] <- calc_weighted_ca(min_us, reps_one_mk[index_seq_al], seq_count_one_mk[index_seq_al])
  }
  return(sr_ca_one_mk)
}

#' extract_sr_use
#'
#' @description The function to extract SR data used
#' @param sr_one_mk A vector of the stutter ratios in one marker
#' @param sr_al_one_mk A vector of the parent alleles of each stutter ratio
#' @param sr_h_one_mk A vector of the total allelic products of each stutter ratio
#' @param lus_one_mk A vector of the LUS of each allele in one marker
extract_sr_use <- function(sr_one_mk, sr_al_one_mk, sr_h_one_mk, lus_one_mk){
  # Remove data without stutter ratio
  index_non_zero <- which(sr_one_mk != 0)
  sr_one_mk <- as.numeric(sr_one_mk[index_non_zero])
  sr_al_one_mk <- as.numeric(sr_al_one_mk[index_non_zero])
  sr_h_one_mk <- as.numeric(sr_h_one_mk[index_non_zero])

  # Assign the LUS value to each allele
  al_names <- as.numeric(names(lus_one_mk))
  sr_lus_one_mk <- rep(-1, length(sr_one_mk))
  for(i in 1:length(lus_one_mk)){
    index_al <- which(sr_al_one_mk == al_names[i])
    if(length(index_al) > 0){
      sr_lus_one_mk[index_al] <- lus_one_mk[i]
    }
  }

  # Remove data without sequence information
  index_no_seq <- which(sr_lus_one_mk == -1)
  if(length(index_no_seq) > 0){
    sr_one_mk <- sr_one_mk[- index_no_seq]
    sr_al_one_mk <- sr_al_one_mk[- index_no_seq]
    sr_h_one_mk <- sr_h_one_mk[- index_no_seq]
    sr_lus_one_mk <- sr_lus_one_mk[- index_no_seq]
  }

  # Make an output object
  sr_use_list <- list(sr_one_mk, sr_al_one_mk, sr_h_one_mk, sr_lus_one_mk)
  names(sr_use_list) <- c("stutter ratio", "parent allele", "total allelic product", "LUS of parent allele")
  return(sr_use_list)
}
