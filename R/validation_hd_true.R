#' generate_rand_gt
#'
#' @description The function to generate random genotypes
#' @param np The number of persons generated
#' @param pop_al_list A list of the alleles in population database
#' @param pop_prob_list A list of the allele probabilities in population database
generate_rand_gt <- function(np, pop_al_list, pop_prob_list){
  mapply(sample, sapply(pop_al_list, as.numeric), 2 * np, prob = pop_prob_list, replace = TRUE)
}

#' generate_importance_gt
#'
#' @description The function to generate genotypes using importance sampling
#' @param np The number of persons generated
#' @param dt_pg_list_h1 A list of data.tables for probabilistic genotyping results in each marker
#' @param pop_al_list A list of the alleles in population database
#' @param pop_prob_list A list of the allele probabilities in population database
generate_importance_gt <- function(np, dt_pg_list_h1, pop_al_list, pop_prob_list){
  n_mk <- length(dt_pg_list_h1)
  hnc <- length(grep("Allele", names(dt_pg_list_h1[[1]]))) / 2
  importance_gt <- matrix(NA, n_mk, 2 * np)
  bias_all_mk <- matrix(NA, n_mk, np)

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
    index_extract <- sample(1:length(prob_sampling), np, prob = prob_sampling, replace = TRUE)
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
#' @param n_simu The number of simulations
#' @param data_deconvo A list of the deconvolution results
#' @param importance_sampling Whether importance sampling will be performed or not
perform_hd_true <- function(hnc, n_simu, data_deconvo, importance_sampling = FALSE){

  lr_rand <- rep(NA, n_simu)
  bias <- rep(NA, n_simu)

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
    dt_pg_list_h1 <- data_h1$dt_pg_list

    # --- Generate random genotypes --- #

    if(importance_sampling){
      tmp <- generate_importance_gt(n_simu, dt_pg_list_h1, pop_al_list, pop_prob_list)
      randgt <- tmp$importance_gt
      bias <- tmp$bias
    }else{
      randgt <- t(generate_rand_gt(n_simu, pop_al_list, pop_prob_list))
    }

    # --- Create order of contributors in Hp and Hd --- #

    # Hd
    hyp_id_all_hd <- make_hyp_id_all(hnc, NULL, "Ref")

    # Hp
    hyp_id_all_hp <- hyp_id_all_hd
    hyp_id_all_hp[length(hyp_id_all_hp)] <- 1 # Assuming 'Ref' as the minor contributor

    # References' IDs included in any hypotheses
    ref_id_lr <- 1

    for(i in 1:n_simu){

      # --- Extract a reference profile --- #

      ref1 <- randgt[, c(2 * i - 1, 2 * i)]
      ref1 <- t(apply(ref1, 1, sort))
      refs <- paste0("Ref", i)
      dt_ref_rand <- data.table(SampleName = rep(refs, length(markers_use)),
                                Marker = markers_use,
                                Allele1 = ref1[, 1],
                                Allele2 = ref1[, 2])

      # --- Calculate log10 likelihood for Hp --- #

      loglike_all_hyp_hp <- calc_loglike_all_hyp(markers_use, peaks_list, pop_al_list, pop_prob_list, theta, dt_pg_list_h1, dt_ref_rand, hyp_id_all_hp, ref_id_lr)
      tmp <- calc_subsub_log10like_all_mk(loglike_all_hyp_hp)
      subsub_log10like_hp <- sum(tmp$subsub_log10like_all_mk)

      # --- Calculate log10 likelihood for Hd --- #

      loglike_all_hyp_hd <- calc_loglike_all_hyp(markers_use, peaks_list, pop_al_list, pop_prob_list, theta, dt_pg_list_h1, dt_ref_rand, hyp_id_all_hd, ref_id_lr)
      tmp <- calc_subsub_log10like_all_mk(loglike_all_hyp_hd)
      subsub_log10like_hd <- sum(tmp$subsub_log10like_all_mk)

      # --- Calculate LR --- #

      lr_rand[i] <- 10^(subsub_log10like_hp - subsub_log10like_hd)
    }
  }

  return(list(lr_rand = lr_rand, bias = bias))
}
