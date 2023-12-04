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
  n_ae_data <- nrow(ae_data)
  n_mk <- ncol(ae_data)
  for(i in 1:n_mk){
    data_list_ae <- list(N = n_ae_data, AE = ae_data[, i], Height = height_ae[, i])

    fit_ae <- model_ae$sample(
      data_list_ae,
      iter_warmup = 10000,
      iter_sampling = 10000,
      chains = 1,
      parallel_chains = 4
    )
  }

  # model_sr_allele <- cmdstan_model("inst/stan/model_sr_allele.stan")
  for(i in 1:n_mk){
    index_bsr <- which(bsr_data[, i] != 0)
    bsr <- bsr_data[index_bsr, i]
    bsr_al <- allele_sr_data[index_bsr, i]
    bsr_h <- height_sr_data[index_bsr, i]

    data_list_bsr <- list(N = length(bsr), SR = bsr, Allele = bsr_al, Height = bsr_h)

    fit_bsr <- model_sr_allele$sample(
      data_list_bsr,
      iter_warmup = 10000,
      iter_sampling = 10000,
      chains = 1,
      parallel_chains = 4
    )
  }

}
