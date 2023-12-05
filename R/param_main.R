test_parameter <- function(){

  ##############
  # Load files #
  ##############

  # Load experimental data
  #dt_exp_init <- fread("C:/Users/manab/Documents/RStudio/KongohPlus/inst/extdata/examples/expData_kakeiken.csv")
  dt_exp_init <- fread("D:/RStudio_GitHub/KongohPlus/inst/extdata/examples/expData_kakeiken.csv")

  # Change data type of dt_exp_init
  dt_exp_init <- change_data_type_exp(dt_exp_init)

  # Load genotype data
  #dt_gt_init <- fread("C:/Users/manab/Documents/RStudio/KongohPlus/inst/extdata/examples/gtData_kakeiken.csv")
  dt_gt_init <- fread("D:/RStudio_GitHub/KongohPlus/inst/extdata/examples/gtData_kakeiken.csv")

  # Change data type of dt_gt_init
  dt_gt_init <- change_data_type_gt(dt_gt_init)

  # Load sequence data
  # Already registerd or user prepared
  #dt_seq_init <- fread("C:/Users/manab/Documents/RStudio/KongohPlus/inst/extdata/examples/Sequence-data_GF_Example.csv")
  dt_seq_init <- fread("D:/RStudio_GitHub/KongohPlus/inst/extdata/examples/Sequence-data_GF_Example.csv")

  # Change data type of dt_seq_init
  dt_seq_init <- change_data_type_seq(dt_seq_init)

  # Check files
  # sample names
  # markers


  ############
  # Settings #
  ############

  # Separate autosomal and sex chromosomal markers
  markers <- unique(dt_exp_init[, Marker])
  sex_markers <- c("Yindel", "AMEL", "DYS391")
  par_markers <- setdiff(markers, sex_markers)

  # Analytical thresholds
  at_all_mk <- c(rep(12, 5), rep(20, 6), rep(13, 4), rep(16, 5), rep(19, 4))
  names(at_all_mk) <- markers

  # Dyes
  dye_all_mk <- c(rep("B", 5), rep("G", 6), rep("Y", 4), rep("R", 5), rep("P", 4))
  names(dye_all_mk) <- markers

  # Consideration of stutters
  sr_consider_all_mk <- matrix(TRUE, length(par_markers), 4)
  rownames(sr_consider_all_mk) <- par_markers
  colnames(sr_consider_all_mk) <- c("BS", "FS", "DS", "M2S")
  sr_consider_all_mk[, 4] <- FALSE
  sr_consider_all_mk[17, 4] <- TRUE
  sr_consider_all_mk[19, 4] <- TRUE

  ############
  # Analysis #
  ############

  # Arrange experimental data
  tmp <- arrange_exp_data(dt_exp_init, dt_gt_init, par_markers, at_all_mk, dye_all_mk, sr_consider_all_mk, max_dist_pu = 1)
  peak_list <- tmp[[1]]
  size_list <- tmp[[2]]
  height_list <- tmp[[3]]
  gt_list <- tmp[[4]]

  # Locus-specific amplification efficiency
  tmp <- get_ae_data(height_list)
  ae_data <- tmp[[1]]
  height_ae <- tmp[[2]]

  # Stutter ratios
  tmp <- get_sr_data(peak_list, height_list, gt_list, sr_consider_all_mk)
  bsr_data <- tmp[[1]]
  fsr_data <- tmp[[2]]
  dsr_data <- tmp[[3]]
  m2sr_data <- tmp[[4]]
  allele_sr_data <- tmp[[5]]
  height_sr_data <- tmp[[6]]

  # Sequence data
  tmp <- analyze_seq(dt_seq_init, par_markers)
  lus_all_mk <- tmp[[1]]
  reps_all_mk <- tmp[[2]]
  seq_al_all_mk <- tmp[[3]]
  seq_count_all_mk <- tmp[[4]]

  ########
  # Stan #
  ########

  # model_ae <- cmdstan_model(paste0(path_pack, "/stan/model_ae.stan"))
  # model_ae <- cmdstan_model("inst/stan/model_ae.stan") # cpp_options = list(stan_threads = TRUE)
  # saveRDS(model_ae, "data/model_ae.rds")
  # readRDS("data/model_ae.rds")
  n_ae_data <- nrow(ae_data)
  n_mk <- ncol(ae_data)
  for(i in 1:n_mk){
    data_list_ae <- list(N = n_ae_data, AE = ae_data[, i], Height = height_ae[, i])

    fit_ae <- model_ae$sample(
      data_list_ae,
      iter_warmup = 10000,
      iter_sampling = 10000,
      chains = 4,
      parallel_chains = 4
    )
  }

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
