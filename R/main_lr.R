#' make_prof_ref_all
#'
#' @description The function to make a matrix for reference profiles of all autosomal markers
#' @param dt_ref The data.table of the reference profile
#' @param markers_auto A vector of the autosomal markers
make_prof_ref_all <- function(dt_ref, markers_auto){
  n_mk <- length(markers_auto)
  refs <- unique(dt_ref[, SampleName])
  n_ref <- length(refs)
  setindex(dt_ref, Marker)

  prof_ref_all <- matrix(0, n_mk, 2 * n_ref)
  index_al1 <- 2 * 1:n_ref - 1
  index_al2 <- 2 * 1:n_ref
  for(i in 1:n_mk){
    mk <- markers_auto[i]
    dt_ref_one_mk <- dt_ref[Marker == mk]
    prof_ref_all[i, index_al1] <- as.numeric(dt_ref_one_mk[, Allele1])
    prof_ref_all[i, index_al2] <- as.numeric(dt_ref_one_mk[, Allele2])
  }
  rownames(prof_ref_all) <- markers_auto
  colnames(prof_ref_all) <- as.vector(sapply(refs, rep, 2))
  return(prof_ref_all)
}

#' make_hyp_id_all
#'
#' @description The function to make ids of all possible hypotheses
#' @param hnc The hypothesized number of contributors
#' @param hyp Known contributors in a hypothesis
#' @param refs The names of all references
make_hyp_id_all <- function(hnc, hyp, refs){
  known_id <- which(is.element(refs, hyp))
  n_known <- length(known_id)
  if(n_known == 0){
    hyp_id_all <- matrix(0, 1, hnc)
  }else{
    index_known <- permutations(n = hnc, r = n_known, repeats.allowed = FALSE)
    n_hyp <- nrow(index_known)
    hyp_id_all <- matrix(0, n_hyp, hnc)
    for(i in 1:n_hyp){
      hyp_id_all[i, index_known[i, ]] <- known_id
    }
  }
  return(hyp_id_all)
}

#' assign_log10like
#'
#' @description The function to assign the log10 likelihood
#' @param hnc The hypothesized number of contributors
#' @param markers_use A vector of the autosomal markers used
#' @param pop_al_list A list of the alleles in population database
#' @param pop_prob_list A list of the allele probabilities in population database
#' @param theta The theta correction to consider sub-population
#' @param gtc_use A list of the genotype combinations in used markers
#' @param dt_pg_list A list of the data.tables for the results of probabilistic genotyping in used markers
#' @param dt_mpdeg The data.table for mixture proportions and degradation parameters of each contributor
#' @param dt_ref The data.table for the reference profiles
#' @param hyp Known contributors in a hypothesis
#' @param refs The names of all references
assign_log10like <- function(hnc, markers_use, pop_al_list, pop_prob_list, theta, gtc_use, dt_pg_list, dt_mpdeg, dt_ref, hyp, refs){
  n_mk <- length(markers_use)

  prof_ref_all <- make_prof_ref_all(dt_ref, markers_use)

  hyp_id_all <- make_hyp_id_all(hnc, hyp, refs)
  n_hyp <- nrow(hyp_id_all)

  gtc_logprob_all <- vector(mode = "list", length = n_hyp)
  for(i in 1:n_hyp){
    gtc_logprob_all[[i]] <- vector(mode = "list", length = n_mk)
  }

  for(i in 1:n_mk){
    # Extract data of one marker
    mk <- markers_use[i]
    prof_ref <- prof_ref_all[mk, ]
    pop_al <- pop_al_list[[mk]]
    pop_prob <- pop_prob_list[[mk]]
    gtc_mk <- gtc_use[[mk]]

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
    log10like_tmp[, i] <- log10like <- calc_log10like_all_mk(lpd_gtc_use, gtc_logprob_all[[i]])
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

  dt_pg_list_hyp <- make_dt_pg_list(markers_use, gtc_use, lpd_gtc_use, gtc_logprob_all[[index_max_like]], label_contrib)

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
#' @param markers_use A vector of the autosomal markers
#' @param log10like_hp A vector of the log10 likelihoods of all markers in Hp
#' @param combined_log10like_hp The combined log10 likelihood in Hp
#' @param log10like_hd A vector of the log10 likelihoods of all markers in Hd
#' @param combined_log10like_hd The combined log10 likelihood in Hd
make_dt_lr <- function(markers_use, log10like_hp, combined_log10like_hp, log10like_hd, combined_log10like_hd){
  dt_lr <- data.table(Marker = c(markers_use, "Total"),
                      Log10_Likelihood_Hp = c(log10like_hp, combined_log10like_hp),
                      Log10_Likelihood_Hd = c(log10like_hd, combined_log10like_hd),
                      LR = c(10^(log10like_hp - log10like_hd), 10^(combined_log10like_hp - combined_log10like_hd)))
  return(dt_lr)
}
