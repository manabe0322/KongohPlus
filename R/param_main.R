#' calc_aic
#'
#' @description The function to calculate the AIC
#' @param log_like A logarithmic likelihood
#' @param n_par The number of parameters
calc_aic <- function(log_like, n_par){
  return(-2 * log_like + 2 * n_par)
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
