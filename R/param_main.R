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
make_dt_par <- function(dt_result_ae, dt_result_bsr, dt_result_fsr, dt_result_dsr, dt_result_m2sr, dt_result_gamma){
  markers <- dt_result_ae[, Marker]
  n_mk <- length(markers)
  dt_par <- data.table(Marker = c(rep(markers, 5), NA),
                       Parameter = c(rep("AE", n_mk), rep("BSR", n_mk), rep("FSR", n_mk), rep("DSR", n_mk), rep("M2SR", n_mk), "Gamma"),
                       Model = c(rep("Log-normal", n_mk), dt_result_bsr[, Model], dt_result_fsr[, Model], dt_result_dsr[, Model], dt_result_m2sr[, Model], "Gamma"),
                       Myu = c(dt_result_ae[, Myu], dt_result_bsr[, Myu], dt_result_fsr[, Myu], dt_result_dsr[, Myu], dt_result_m2sr[, Myu], NA),
                       Slope = c(rep(NA, n_mk), dt_result_bsr[, Slope], dt_result_fsr[, Slope], dt_result_dsr[, Slope], rep(NA, n_mk), NA),
                       Intercept = c(rep(NA, n_mk), dt_result_bsr[, Intercept], dt_result_fsr[, Intercept], dt_result_dsr[, Intercept], rep(NA, n_mk), NA),
                       Sigma = c(rep(NA, n_mk), dt_result_bsr[, Sigma], dt_result_fsr[, Sigma], dt_result_dsr[, Sigma], rep(NA, n_mk), NA),
                       MinUS = c(rep(NA, n_mk), dt_result_bsr[, MinUS], dt_result_fsr[, MinUS], dt_result_dsr[, MinUS], rep(NA, n_mk), NA),
                       Beta = c(rep(NA, 5 * n_mk), dt_result_gamma[, Beta]),
                       Rhat_Myu = c(dt_result_ae[, Rhat_Myu], dt_result_bsr[, Rhat_Myu], dt_result_fsr[, Rhat_Myu], dt_result_dsr[, Rhat_Myu], dt_result_m2sr[, Rhat_Myu], NA),
                       Rhat_Slope = c(rep(NA, n_mk), dt_result_bsr[, Rhat_Slope], dt_result_fsr[, Rhat_Slope], dt_result_dsr[, Rhat_Slope], rep(NA, n_mk), NA),
                       Rhat_Intercept = c(rep(NA, n_mk), dt_result_bsr[, Rhat_Intercept], dt_result_fsr[, Rhat_Intercept], dt_result_dsr[, Rhat_Intercept], rep(NA, n_mk), NA),
                       Rhat_Sigma = c(rep(NA, n_mk), dt_result_bsr[, Rhat_Sigma], dt_result_fsr[, Rhat_Sigma], dt_result_dsr[, Rhat_Sigma], rep(NA, n_mk), NA),
                       Rhat_Beta = c(rep(NA, 5 * n_mk), dt_result_gamma[, Rhat_Beta]),
                       Rhat_Likelihood = c(dt_result_ae[, Rhat_Likelihood], dt_result_bsr[, Rhat_Likelihood], dt_result_fsr[, Rhat_Likelihood], dt_result_dsr[, Rhat_Likelihood], dt_result_m2sr[, Rhat_Likelihood], dt_result_gamma[, Rhat_Likelihood]))
  return(dt_par)
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
#' @param min_n_motif_lus The minimum number of repeat motifs for considering the LUS model
#' @param min_n_motif_multi The minimum number of repeat motifs for considering the multi-sequence model
estimate_parameters <- function(dt_exp, dt_gt, dt_seq, dt_kit_set,
                                model_ae, model_sr_gen, model_sr_together, model_sr_uni, model_gamma,
                                n_warmup, n_sampling, n_chains,
                                min_n_motif_lus, min_n_motif_multi, max_n_motif_multi = 10,
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
  ae_result_list <- analyze_ae_main(markers_par, model_ae,
                                    ae_data, height_ae_data,
                                    n_warmup, n_sampling, n_chains, n_parallel_chains)

  # Perform the MCMC for the back stutter ratio
  bsr_result_list <- analyze_bfdstr_main("BSR", markers_par, model_bsr_all_mk, model_sr_gen, model_sr_together,
                                         bsr_data, allele_sr_data, height_sr_data,
                                         lus_all_mk, reps_all_mk, seq_al_all_mk, seq_count_all_mk,
                                         n_warmup, n_sampling, n_chains, n_parallel_chains, min_n_motif_lus, min_n_motif_multi, max_n_motif_multi)

  # Perform the MCMC for the forward stutter ratio
  fsr_result_list <- analyze_bfdstr_main("FSR", markers_par, model_fsr_all_mk, model_sr_gen, model_sr_together,
                                         fsr_data, allele_sr_data, height_sr_data,
                                         lus_all_mk, reps_all_mk, seq_al_all_mk, seq_count_all_mk,
                                         n_warmup, n_sampling, n_chains, n_parallel_chains, min_n_motif_lus, min_n_motif_multi, max_n_motif_multi)

  # Perform the MCMC for the double-back stutter ratio
  dsr_result_list <- analyze_bfdstr_main("DSR", markers_par, model_dsr_all_mk, model_sr_gen, model_sr_together,
                                         dsr_data, allele_sr_data, height_sr_data,
                                         lus_all_mk, reps_all_mk, seq_al_all_mk, seq_count_all_mk,
                                         n_warmup, n_sampling, n_chains, n_parallel_chains, min_n_motif_lus, min_n_motif_multi, max_n_motif_multi)

  # Perform the MCMC for the minus 2-nt stutter ratio
  m2sr_result_list <- analyze_m2sr_main(markers_par, model_m2sr_all_mk, model_sr_uni, model_sr_together,
                                        m2sr_data, allele_sr_data, height_sr_data,
                                        lus_all_mk,
                                        n_warmup, n_sampling, n_chains, n_parallel_chains)

  #### Gamma parameter ####

  cat("MCMC for the gamma parameter\n")

  # Calculate expected peak heights
  dt_result_ae <- ae_result_list$dt_result
  dt_result_bsr <- bsr_result_list$dt_result
  dt_result_fsr <- fsr_result_list$dt_result
  dt_result_dsr <- dsr_result_list$dt_result
  dt_result_m2sr <- m2sr_result_list$dt_result
  dt_alcor_bsr <- bsr_result_list$dt_alcor
  dt_alcor_fsr <- fsr_result_list$dt_alcor
  dt_alcor_dsr <- dsr_result_list$dt_alcor
  eph_all_data <- calc_eph_all_data(markers_par, n_oph, tph_all_data, gt_list, peak_list,
                                    dt_result_ae, dt_result_bsr, dt_result_fsr, dt_result_dsr, dt_result_m2sr,
                                    dt_alcor_bsr, dt_alcor_fsr, dt_alcor_dsr)

  # Remove data without sequence information (i.e., NA data)
  index_adopt <- which(!is.na(eph_all_data))
  oph_all_data <- oph_all_data[index_adopt]
  eph_all_data <- eph_all_data[index_adopt]

  # Perform MCMC
  data_list_gamma <- list(N = length(oph_all_data), OPH = oph_all_data, EPH = eph_all_data)
  tmp <- run_mcmc_gamma(model_gamma, data_list_gamma, n_warmup, n_sampling, n_chains, n_parallel_chains)
  dt_result_gamma <- data.table(Beta = tmp$median_beta,
                                Rhat_Beta = tmp$rhat_beta,
                                Rhat_Likelihood = tmp$rhat_lp)
  gamma_result_list <- list(dt_result = dt_result_gamma,
                            traceplot = tmp$traceplot)

  #### Results of estimated parameters ####

  dt_par <- make_dt_par(dt_result_ae, dt_result_bsr, dt_result_fsr, dt_result_dsr, dt_result_m2sr, dt_result_gamma)

  return(list("markers_par" = markers_par,
              "ae_result_list" = ae_result_list,
              "bsr_result_list" = bsr_result_list,
              "fsr_result_list" = fsr_result_list,
              "dsr_result_list" = dsr_result_list,
              "m2sr_result_list" = m2sr_result_list,
              "gamma_result_list" = gamma_result_list,
              "dt_par" = dt_par))
}
