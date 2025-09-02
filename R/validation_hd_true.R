#' check_gt_non_contrib
#'
#' @description The function to check the loaded genotypes for non-contributors
#' @param dt_gt_non_contrib The data.table of the genotypes for non-contributors
#' @param markers_use A vector of the autosomal markers used
check_gt_non_contrib <- function(dt_gt_non_contrib, markers_use){
  error_message <- ""

  setindex(dt_gt_non_contrib, SampleName)
  sn_all_non_contrib <- unique(dt_gt_non_contrib[, SampleName])
  for(i in 1:length(sn_all_non_contrib)){
    sn1 <- sn_all_non_contrib[i]
    dt_gt_non_contrib_sn1 <- dt_gt_non_contrib[SampleName == sn1]
    if(!all(is.element(markers_use, dt_gt_non_contrib_sn1[, Marker]))){
      if(error_message == ""){
        error_message <- paste0("The following non-contributors lacks genotype information for the markers used in the deconvolution: <br>- ", sn1, "<br>")
      }else{
        error_message <- paste0(error_message, "- ", sn1, "<br>")
      }
    }
  }

  return(error_message)
}

#' generate_rand_gt
#'
#' @description The function to generate random genotypes
#' @param n_simu The number of persons generated
#' @param pop_al_list A list of the alleles in population database
#' @param pop_prob_list A list of the allele probabilities in population database
generate_rand_gt <- function(n_simu, pop_al_list, pop_prob_list){
  mapply(sample, sapply(pop_al_list, as.numeric), 2 * n_simu, prob = pop_prob_list, replace = TRUE)
}

#' create_dt_rand_gt
#'
#' @description The function to create the data.table for random genotypes
#' @param rand_gt The output object of the function 'generate_rand_gt'
#' @param markers_use A vector of the autosomal markers used
create_dt_rand_gt <- function(rand_gt, markers_use){
  n_simu <- ncol(rand_gt) / 2

  alleles1 <- as.vector(rand_gt[, seq(1, ncol(rand_gt), 2), drop = FALSE])
  alleles2 <- as.vector(rand_gt[, seq(2, ncol(rand_gt), 2), drop = FALSE])
  mat_alleles <- cbind(alleles1, alleles2)
  mat_alleles <- t(apply(mat_alleles, 1, sort))

  dt_rand_gt <- data.table(SampleName = as.vector(sapply(paste0("Sample", 1:n_simu), rep, length(markers_use))),
                           Marker = rep(markers_use, n_simu),
                           Allele1 = mat_alleles[, 1],
                           Allele2 = mat_alleles[, 2])
  return(dt_rand_gt)
}

#' generate_importance_gt
#'
#' @description The function to generate genotypes using importance sampling
#' @param n_simu The number of persons generated
#' @param dt_pg_list_h1 A list of data.tables for probabilistic genotyping results in each marker
#' @param pop_al_list A list of the alleles in population database
#' @param pop_prob_list A list of the allele probabilities in population database
generate_importance_gt <- function(n_simu, dt_pg_list_h1, pop_al_list, pop_prob_list){
  n_mk <- length(dt_pg_list_h1)
  hnc <- length(grep("Allele", names(dt_pg_list_h1[[1]]))) / 2
  importance_gt <- matrix(NA, n_mk, 2 * n_simu)
  bias_all_mk <- matrix(NA, n_mk, n_simu)

  for(i in 1:n_mk){

    # --- Change the format of the pg table from per combination to a minor contributor --- #

    dt_pg <- dt_pg_list_h1[[i]]
    dt_pg_per_contrib <- aggregate_pg_contrib(dt_pg, paste0("U", hnc))

    # --- Calculate genotype probabilities --- #

    tmp_gt <- gsub("Q", "99", dt_pg_per_contrib[[paste0("U", hnc)]])
    tmp_gt <- strsplit(tmp_gt, ", ")
    gt_mk <- lapply(tmp_gt, as.numeric)
    pop_al <- pop_al_list[[i]]
    pop_prob <- pop_prob_list[[i]]
    n_obs <- make_n_obs(pop_al, prof_ref = numeric(0))
    obs_al <- as.numeric(names(n_obs))
    gt_prob <- calc_gts_prob_mk(gt_mk, pop_al, pop_prob, theta = 0, prof_ref_99 = numeric(0), hyp_id = 0, n_obs, obs_al, index_dropout_cpp = integer(0))

    # --- Randomly select genotype combinations based on 'Probability' --- #

    prob_sampling <- dt_pg_per_contrib$Probability
    index_extract <- sample(1:length(prob_sampling), n_simu, prob = prob_sampling, replace = TRUE)
    importance_gt[i, ] <- unlist(gt_mk[index_extract])

    # --- Calculate bias --- #

    bias_all_mk[i, ] <- gt_prob[index_extract] / prob_sampling[index_extract]
  }

  return(list(importance_gt = importance_gt, bias = apply(bias_all_mk, 2, prod)))
}

#' perform_hd_true
#'
#' @description The function to perform Hd-true test
#' @param hnc The targeted number of contributors
#' @param data_deconvo A list of the deconvolution results
#' @param n_simu The number of simulations
#' @param dt_gt_non_contrib The data.table of the genotypes for non-contributors
#' @param importance_sampling Whether importance sampling will be performed or not
perform_hd_true <- function(hnc, data_deconvo, n_simu, dt_gt_non_contrib = NULL, importance_sampling = FALSE){

  # --- Extract data from data_deconvo --- #

  hnc_all <- data_deconvo$hnc_all
  markers_use <- data_deconvo$markers_use
  peaks_list <- data_deconvo$peaks_list
  pop_al_list <- data_deconvo$pop_al_list
  pop_prob_list <- data_deconvo$pop_prob_list
  theta <- data_deconvo$theta
  data_hnc_list <- data_deconvo$data_hnc_list

  data_h1 <- data_hnc_list[[which(hnc_all == hnc)]]
  error_message_h1 <- data_h1$error_message

  if(is.null(error_message_h1)){

    lr_rand <- matrix(NA, length(markers_use) + 1, n_simu)
    bias <- rep(NA, n_simu)

    # --- Generate random genotypes --- #

    dt_pg_list_h1 <- data_h1$dt_pg_list

    if(importance_sampling){
      tmp <- generate_importance_gt(n_simu, dt_pg_list_h1, pop_al_list, pop_prob_list)
      dt_gt_non_contrib <- create_dt_rand_gt(tmp$importance_gt, markers_use)
      bias <- tmp$bias
    }else{
      if(is.null(dt_gt_non_contrib)){
        tmp <- set_al_prob_randgt(markers_use, data_deconvo$dt_af)
        rand_gt <- t(generate_rand_gt(n_simu, tmp$pop_al_list, tmp$pop_prob_list))
        dt_gt_non_contrib <- create_dt_rand_gt(rand_gt, markers_use)
      }else{
        dt_gt_non_contrib <- dt_gt_non_contrib[is.element(dt_gt_non_contrib$Marker, markers_use)]
      }
    }

    # --- Create order of contributors in Hp and Hd --- #

    # Hp
    hyp_id_all_hp <- make_hyp_id_all(hnc, "Ref", "Ref")

    # Hd
    hyp_id_all_hd <- make_hyp_id_all(hnc, NULL, "Ref")

    # References' IDs included in any hypotheses
    ref_id_lr <- 1

    setindex(dt_gt_non_contrib, SampleName)
    sn_all_non_contrib <- unique(dt_gt_non_contrib[, SampleName])
    for(i in 1:n_simu){

      # --- Extract a reference profile --- #

      dt_ref_rand <- dt_gt_non_contrib[SampleName == sn_all_non_contrib[i]]

      # --- Calculate log10 likelihood for Hp --- #

      loglike_all_hyp_hp <- calc_loglike_all_hyp(markers_use, peaks_list, pop_al_list, pop_prob_list, theta, dt_pg_list_h1, dt_ref_rand, hyp_id_all_hp, ref_id_lr)
      tmp <- calc_subsub_log10like_all_mk(loglike_all_hyp_hp)
      subsub_log10like_all_mk_hp <- tmp$subsub_log10like_all_mk
      subsub_log10like_hp <- c(subsub_log10like_all_mk_hp, sum(subsub_log10like_all_mk_hp))

      # --- Calculate log10 likelihood for Hd --- #

      loglike_all_hyp_hd <- calc_loglike_all_hyp(markers_use, peaks_list, pop_al_list, pop_prob_list, theta, dt_pg_list_h1, dt_ref_rand, hyp_id_all_hd, ref_id_lr)
      tmp <- calc_subsub_log10like_all_mk(loglike_all_hyp_hd)
      subsub_log10like_all_mk_hd <- tmp$subsub_log10like_all_mk
      subsub_log10like_hd <- c(subsub_log10like_all_mk_hd, sum(subsub_log10like_all_mk_hd))

      # --- Calculate LR --- #

      lr_rand[, i] <- 10^(subsub_log10like_hp - subsub_log10like_hd)
    }

    lr_rand <- cbind(c(markers_use, "LR"), lr_rand)
    dt_lr_rand <- as.data.table(lr_rand)
    names(dt_lr_rand) <- c("Marker", sn_all_non_contrib)
  }else{
    dt_lr_rand <- NULL
    bias <- NULL
    dt_gt_non_contrib <- NULL
  }

  return(list(dt_lr_rand = dt_lr_rand, bias = bias, dt_gt_non_contrib = dt_gt_non_contrib))
}
