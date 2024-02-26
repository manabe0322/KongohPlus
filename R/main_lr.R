assign_log10like <- function(data_deconvo, dt_ref, hyp, refs){
  # Extract data from data_deconvo
  hnc <- data_deconvo$hnc
  markers_auto <- data_deconvo$markers_auto
  pop_al_list <- data_deconvo$pop_al_list
  pop_prob_list <- data_deconvo$pop_prob_list
  theta <- data_deconvo$theta
  gtc_adopt <- data_deconvo$gtc_adopt
  lpd_gtc_adopt <- data_deconvo$lpd_gtc_adopt

  n_mk <- length(markers_auto)

  prof_ref_all <- make_prof_ref_all(dt_ref, markers_auto)

  hyp_id_all <- make_hyp_id_all(hnc, hyp, refs)
  n_hyp <- nrow(hyp_id_all)

  gtc_logprob_all <- vector(mode = "list", length = n_hyp)
  for(i in 1:n_hyp){
    gtc_logprob_all[[i]] <- vector(mode = "list", length = n_mk)
  }

  for(i in 1:n_mk){
    # Extract data of one marker
    mk <- markers_auto[i]
    prof_ref <- prof_ref_all[mk, ]
    pop_al <- pop_al_list[[mk]]
    pop_prob <- pop_prob_list[[mk]]
    gtc_mk <- gtc_adopt[[mk]]
    gtc_mk <- asplit(gtc_mk, 1)

    # Calculate genotype probabilities
    n_obs <- make_n_obs(pop_al, prof_ref)
    for(j in 1:n_hyp){
      gtc_prob_mk <- calc_gtc_prob_mk(gtc_mk, pop_al, pop_prob, theta, prof_ref, hyp_id_all[j, ], n_obs)
      gtc_logprob_all[[j]][[i]] <- log(gtc_prob_mk)
    }
  }

  log10like_tmp <- matrix(0, n_mk + 1, n_hyp)
  for(i in 1:n_hyp){
    log10like <- calc_log10like_all_mk(lpd_gtc_adopt, gtc_logprob_all[[i]])
    combined_log10like <- calc_combined_log10like(log10like)
    log10like_tmp[, i] <- c(log10like, combined_log10like)

    print(gtc_logprob_all[[i]])
  }

  return(log10like_tmp[, which.max(log10like_tmp[n_mk + 1, ])])
}
