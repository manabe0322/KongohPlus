#' estimate_parameters
#'
#' @description The function to estimate parameters
#' @param dt_exp The data.table of the experimental data
#' @param dt_gt The data.table of the genotype data
#' @param dt_seq The data.table of the sequence data
#' @param dt_kit_set The data.table of the kit-setting
estimate_parameters <- function(dt_exp, dt_gt, dt_seq, dt_kit_set){
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

  # Load the MCMC models
  model_ae <- readRDS(paste0(path_pack, "/model/model_ae.rds"))
  model_sr_gen <- readRDS(paste0(path_pack, "/model/model_sr_gen.rds"))
  model_sr_together <- readRDS(paste0(path_pack, "/model/model_sr_together.rds"))
  model_sr_uni <- readRDS(paste0(path_pack, "/model/model_sr_uni.rds"))

  # Perform the MCMC for the locus-specific amplification efficiency
  n_ae_data <- nrow(ae_data)
  ae_median_myu_all_mk <- ae_median_sigma_all_mk <- ae_rhat_myu_all_mk <- ae_rhat_sigma_all_mk <- rep(NA, n_mk)
  ae_traceplot_all_mk <- list()
  for(i in 1:n_mk){
    data_list_ae <- list(N = n_ae_data, AE = ae_data[, i], Height = height_ae_data[, i])
    tmp <- run_mcmc_ae(model_ae, data_list_ae, n_warmup, n_sampling, n_chains, n_parallel_chains)
    ae_median_myu_all_mk[i] <- tmp$median_myu
    ae_median_sigma_all_mk[i] <- tmp$median_sigma
    ae_rhat_myu_all_mk[i] <- tmp$rhat_myu
    ae_rhat_sigma_all_mk[i] <- tmp$rhat_sigma
    ae_traceplot_all_mk[[i]] <- tmp$traceplot
  }

  # Perform the MCMC for the back stutter ratio
  model_bsr_all_mk

  # Perform the MCMC for the forward stutter ratio

  # Perform the MCMC for the double-back stutter ratio

  # Perform the MCMC for the minus 2-nt stutter ratio
}
