assign_log10like <- function(data_deconvo, dt_ref, hyp, refs){
  # Extract data from data_deconvo
  hnc <- data_deconvo$hnc
  markers_auto <- data_deconvo$markers_auto
  pop_al_list <- data_deconvo$pop_al_list
  pop_prob_list <- data_deconvo$pop_prob_list
  theta <- data_deconvo$theta
  gtc_adopt <- data_deconvo$gtc_adopt
  lpd_gtc_adopt <- data_deconvo$lpd_gtc_adopt
  dt_mpdeg <- data_deconvo$dt_mpdeg

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

    # Calculate genotype probabilities
    n_obs <- make_n_obs(pop_al, prof_ref)
    for(j in 1:n_hyp){
      gtc_prob_mk <- calc_gtc_prob_mk(gtc_mk, pop_al, pop_prob, theta, prof_ref, hyp_id_all[j, ], n_obs)
      gtc_logprob_all[[j]][[i]] <- log(gtc_prob_mk)
    }
  }

  log10like_tmp <- matrix(0, n_mk, n_hyp)
  combined_log10like_tmp <- rep(0, n_hyp)
  for(i in 1:n_hyp){
    log10like_tmp[, i] <- log10like <- calc_log10like_all_mk(lpd_gtc_adopt, gtc_logprob_all[[i]])
    combined_log10like <- calc_combined_log10like(log10like)
    combined_log10like_tmp[i] <- combined_log10like
  }

  index_max_like <- which.max(combined_log10like_tmp)
  log10like_adopt <- log10like_tmp[, index_max_like]
  combined_log10like_adopt <- combined_log10like_tmp[index_max_like]

  hyp_id_adopt <- hyp_id_all[index_max_like, ]

  label_contrib <- make_label_contrib(refs, hyp_id_adopt)
  dt_mpdeg_hyp <- data.table(Contributor = label_contrib,
                             MP = dt_mpdeg[, MP],
                             Deg = dt_mpdeg[, Deg])

  dt_pg_list_hyp <- make_dt_pg_list(markers_auto, gtc_adopt, lpd_gtc_adopt, gtc_logprob_all[[index_max_like]], label_contrib)

  return(list("log10like_adopt" = log10like_adopt,
              "combined_log10like_adopt" = combined_log10like_adopt,
              "dt_mpdeg_hyp" = dt_mpdeg_hyp,
              "dt_pg_list_hyp" = dt_pg_list_hyp))
}

#' make_label_contrib
#'
#' @description The function to make labels of each contributor for assigning them to MP and Deg
#' @param refs The vector of all references
#' @param hyp_id A vector of IDs in a hypothesis
make_label_contrib <- function(refs, hyp_id){
  hnc <- length(hyp_id)

  # Reverse sequence of 'hyp_id' (minor to major -> major to minor)
  hyp_id <- hyp_id[hnc:1]

  label_contrib <- rep("", hnc)
  index_unk <- which(hyp_id == 0)
  n_unk <- length(index_unk)
  if(n_unk != 0){
    label_contrib[index_unk] <- paste0("U", 1:n_unk)
  }

  index_k <- which(hyp_id != 0)
  n_k <- length(index_k)
  if(n_k != 0){
    label_contrib[index_k] <- refs[hyp_id[index_k]]
  }

  return(label_contrib)
}

#' make_dt_lr
#'
#' @description The function to make a data.table of the LR
#' @param markers_auto A vector of the autosomal markers
#' @param log10like_hp A vector of the log10 likelihoods of all markers in Hp
#' @param combined_log10like_hp The combined log10 likelihood in Hp
#' @param log10like_hd A vector of the log10 likelihoods of all markers in Hd
#' @param combined_log10like_hd The combined log10 likelihood in Hd
make_dt_lr <- function(markers_auto, log10like_hp, combined_log10like_hp, log10like_hd, combined_log10like_hd){
  dt_lr <- data.table(Marker = c(markers_auto, "Total"),
                      Log10_Likelihood_Hp = c(log10like_hp, combined_log10like_hp),
                      Log10_Likelihood_Hd = c(log10like_hd, combined_log10like_hd),
                      LR = c(10^(log10like_hp - log10like_hd), 10^(combined_log10like_hp - combined_log10like_hd)))
  return(dt_lr)
}
