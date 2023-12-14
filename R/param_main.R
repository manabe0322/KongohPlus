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

  # Perform the MCMC for the locus-specific amplification efficiency
  model_ae <- readRDS("data/model_ae.rds")
  mcmc_data_ae <- run_mcmc_ae(model_ae, ae_data, height_ae_data,
                              n_warmup, n_sampling, n_chains, n_parallel_chains)

  # Perform the MCMC for the back stutter ratio

  # model_sr_gen <- cmdstan_model("inst/stan/model_sr_gen.stan")
  # saveRDS(model_sr_gen, "data/model_sr_gen.rds")
  # model_sr_gen <- readRDS("data/model_sr_gen.rds")
  for(i in 1:n_mk){
    # Extract data of one marker
    tmp <- extract_sr_use(bsr_data[, i], allele_sr_data[, i], height_sr_data[, i], lus_all_mk[[i]])
    bsr_one_mk <- tmp[[1]]
    bsr_al_one_mk <- tmp[[2]]
    bsr_h_one_mk <- tmp[[3]]
    bsr_lus_one_mk <- tmp[[4]]
    reps_one_mk <- reps_all_mk[[i]]
    seq_al_one_mk <- seq_al_all_mk[[i]]
    seq_count_one_mk <- seq_count_all_mk[[i]]

    # Allele model
    data_list_bsr_allele <- list(N = length(bsr_one_mk), SR = bsr_one_mk, XVal = bsr_al_one_mk, Height = bsr_h_one_mk)
    fit_bsr_allele <- model_sr_gen$sample(
      data_list_bsr_allele,
      iter_warmup = 1000,
      iter_sampling = 1000,
      chains = 4,
      parallel_chains = 4
    )

    # LUS model
    data_list_bsr_lus <- list(N = length(bsr_one_mk), SR = bsr_one_mk, XVal = bsr_lus_one_mk, Height = bsr_h_one_mk)
    fit_bsr_lus <- model_sr_gen$sample(
      data_list_bsr_lus,
      iter_warmup = 10000,
      iter_sampling = 10000,
      chains = 4,
      parallel_chains = 4
    )

    # Multi-seq model
    max_cand_us <- max(unlist(reps_one_mk))
    for(j in 1:max_cand_us){
      min_us <- j
      bsr_ca_one_mk <- make_sr_ca_one_mk(min_us, reps_one_mk, seq_al_one_mk, seq_count_one_mk, bsr_al_one_mk)
      data_list_bsr_multi <- list(N = length(bsr_one_mk), SR = bsr_one_mk, XVal = bsr_ca_one_mk, Height = bsr_h_one_mk)
      fit_bsr_multi <- model_sr_gen$sample(
        data_list_bsr_multi,
        iter_warmup = 10000,
        iter_sampling = 10000,
        chains = 4,
        parallel_chains = 4,
        show_messages = FALSE,
        show_exceptions = FALSE
      )

      print(fit_bsr_multi)
    }
  }
}
