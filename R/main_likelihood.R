#' extract_prof_ref
#'
#' @description The function to extract reference profiles of one marker
#' @param dt_ref The data.table of the reference profile
#' @param mk The target marker
#' @param refs The reference names
extract_prof_ref <- function(dt_ref, mk, refs){
  dt_ref_one_mk <- dt_ref[dt_ref$Marker == mk]
  n_ref <- length(refs)
  prof_ref <- rep(-1, 2 * n_ref)
  for(i in 1:n_ref){
    index_ref <- which(dt_ref_one_mk$SampleName == refs[i])
    if(length(index_ref) != 0){
      tmp <- dt_ref_one_mk[index_ref, ]
      a1 <- as.numeric(tmp$Allele1)
      if(!is.na(a1)){
        prof_ref[2 * i - 1] <- a1
      }

      a2 <- as.numeric(tmp$Allele2)
      if(!is.na(a2)){
        prof_ref[2 * i] <- a2
      }
    }
  }
  return(prof_ref)
}

#' replace_prof_ref_99
#'
#' @description The function to replace drop-out alleles with 99 for reference profiles in one marker
#' @param prof_ref Reference profiles of one marker
#' @param peaks A vector of peaks (observed peaks + 99) in one marker
replace_prof_ref_99 <- function(prof_ref, peaks){
  prof_ref_99 <- prof_ref
  prof_ref_99[!is.element(prof_ref, peaks)] <- 99

  n_ref <- length(prof_ref) / 2
  for(i in 1:n_ref){
    prof_ref_99[c(2 * i - 1, 2 * i)] <- sort(prof_ref_99[c(2 * i - 1, 2 * i)])
  }

  return(prof_ref_99)
}

#' create_prof_ref_all
#'
#' @description The function to create two matrices for reference profiles of all autosomal markers
#' @param dt_ref The data.table of the reference profile
#' @param markers_auto A vector of the autosomal markers
#' @param peaks_list A list of peaks
create_prof_ref_all <- function(dt_ref, markers_auto, peaks_list){
  n_mk <- length(markers_auto)
  refs <- unique(dt_ref$SampleName)
  n_ref <- length(refs)
  prof_ref_all <- matrix(-1, n_mk, 2 * n_ref)
  rownames(prof_ref_all) <- markers_auto
  colnames(prof_ref_all) <- as.vector(sapply(refs, rep, 2))
  prof_ref_99_all <- prof_ref_all
  for(i in 1:n_mk){
    mk <- markers_auto[i]
    prof_ref_all[i, ] <- extract_prof_ref(dt_ref, mk, refs)
    prof_ref_99_all[i, ] <- replace_prof_ref_99(prof_ref_all[i, ], peaks_list[[mk]])
  }
  return(list("prof_ref_all" = prof_ref_all,
              "prof_ref_99_all" = prof_ref_99_all))
}

#' make_hyp_id_all
#'
#' @description The function to make ids of all possible hypotheses
#' @param hnc The assumed number of contributors
#' @param hyp The contributor names in a hypothesis
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

#' select_hyp_id_assumed_contrib
#'
#' @description The function to select hyp_id fixing assumed contributors' position
#' @param hyp_id_all Contributors' IDs of all possible hypotheses
#' @param hyp_id_best_assumed_contrib The best hyp_id when considering assumed contributors
select_hyp_id_assumed_contrib <- function(hyp_id_all, hyp_id_best_assumed_contrib){
  index_assumed_contrib <- which(hyp_id_best_assumed_contrib != 0)
  n_assumed_contrib <- length(index_assumed_contrib)
  bool_mat_hyp <- matrix(TRUE, nrow = nrow(hyp_id_all), ncol = n_assumed_contrib)
  for(i in 1:n_assumed_contrib){
    index_tmp <- index_assumed_contrib[i]
    bool_mat_hyp[, i] <- hyp_id_all[, index_tmp] == hyp_id_best_assumed_contrib[index_tmp]
  }
  hyp_id_all_new <- hyp_id_all[apply(bool_mat_hyp, 1, all), , drop = FALSE]
  return(hyp_id_all_new)
}

#' make_n_obs
#'
#' @description The function to make an object for the number of observed alleles in reference profiles
#' @param pop_al A vector of the alleles in the population database
#' @param prof_ref_lr A vector of the reference alleles in one marker. The references are included in any hypotheses.
make_n_obs <- function(pop_al, prof_ref_lr){
  prob_al <- sort(unique(c(pop_al, prof_ref_lr)))
  n_obs <- rep(0, length(prob_al))
  names(n_obs) <- prob_al
  if(length(prof_ref_lr) != 0){
    for(i in 1:length(prof_ref_lr)){
      index_al <- which(prob_al == prof_ref_lr[i])
      n_obs[index_al] <- n_obs[index_al] + 1
    }
  }
  return(n_obs)
}

#' calc_loglike_r
#'
#' @description The function to calculate the log likelihood of all markers
#' @param dt_pg_list A list of the PG data.table for all markers
#' @param gtc_logprob_all_mk A list of the lob probabilities of each genotype combination for all markers when assuming POI or known contributors
#' @param index_hyp The index of the targeted hypothesis
calc_loglike_r <- function(dt_pg_list, gtc_logprob_all_mk, index_hyp){
  n_mk <- length(dt_pg_list)
  loglike_all_mk <- rep(0, n_mk)
  names(loglike_all_mk) <- names(dt_pg_list)

  # log10 likelihood in each marker
  n_rep <- length(grep("Loglike_ph_rep", names(dt_pg_list[[1]])))
  for(i in 1:n_mk){
    dt_pg <- dt_pg_list[[i]]
    loglike_gtc <- gtc_logprob_all_mk[[i]][[index_hyp]]

    # With reference genotype(s)
    if(!is.null(loglike_gtc)){
      for(j in 1:n_rep){
        loglike_gtc <- loglike_gtc + dt_pg[[paste0("Loglike_ph_rep", j)]]
      }
      loglike_mk <- log_sum_exp(loglike_gtc)
      if(is.nan(loglike_mk)){
        loglike_all_mk[i] <- -Inf
      }else{
        loglike_all_mk[i] <- loglike_mk
      }
    }
  }

  return(loglike_all_mk)
}

#' calc_gtc_logprob_all_hyp
#'
#' @description The function to calculate log-probability of each genotype combination for all hypotheses
#' @param mk The marker name
#' @param hyp_id_all Contributors' IDs of all possible hypotheses
#' @param prof_ref_all The object is created by the function 'create_prof_ref_all'
#' @param prof_ref_99_all The object is created by the function 'create_prof_ref_all'
#' @param peaks_list A list of peaks
#' @param pop_al_list A list of the alleles in population database
#' @param pop_prob_list A list of the allele probabilities in population database
#' @param dt_pg_list A list of the PG data.table for all markers
#' @param theta The theta value to consider sub-population effect
#' @param ref_id_lr References' IDs included in any hypotheses
calc_gtc_logprob_all_hyp <- function(mk, hyp_id_all, prof_ref_all, prof_ref_99_all, peaks_list, pop_al_list, pop_prob_list, dt_pg_list, theta, ref_id_lr){
  n_hyp <- nrow(hyp_id_all)
  gtc_logprob_all_hyp <- vector(mode = "list", length = n_hyp)

  # Check reference profiles of assumed contributors are obtained or not
  prof_ref <- prof_ref_all[mk, ]
  index_ref_contrib <- sort(c(2 * ref_id_lr - 1, 2 * ref_id_lr))
  if(all(prof_ref[index_ref_contrib] != -1)){
    # Extract data of one marker
    prof_ref_99 <- prof_ref_99_all[mk, ]
    peaks <- peaks_list[[mk]]
    pop_al <- pop_al_list[[mk]]
    pop_prob <- pop_prob_list[[mk]]
    dt_pg <- dt_pg_list[[mk]]

    # Create observed allele lists
    n_obs <- make_n_obs(pop_al, prof_ref[index_ref_contrib])
    obs_al <- as.numeric(names(n_obs))

    # Indices of drop-out alleles
    dropout_alleles <- setdiff(prof_ref[index_ref_contrib], peaks)
    index_dropout_cpp <- which(is.element(obs_al, dropout_alleles)) - 1

    # Make a list of genotype combinations for C++ format
    gtc_mk <- as.matrix(dt_pg[, grep("Allele", names(dt_pg), value = TRUE), with = FALSE])
    gtc_mk <- split(gtc_mk, row(gtc_mk))

    # Calculate genotype probabilities
    for(i in 1:n_hyp){
      gtc_prob_mk <- calc_gts_prob_mk(gtc_mk, pop_al, pop_prob, theta, prof_ref_99, hyp_id_all[i, ], n_obs, obs_al, index_dropout_cpp)
      gtc_logprob_all_hyp[[i]] <- log(gtc_prob_mk)
    }
  }

  return(gtc_logprob_all_hyp)
}

#' calc_loglike_all_hyp
#'
#' @description The function to assign the log10 likelihood
#' @param markers_use A vector of the autosomal markers used
#' @param peaks_list A list of peaks
#' @param pop_al_list A list of the alleles in population database
#' @param pop_prob_list A list of the allele probabilities in population database
#' @param theta The theta correction to consider sub-population
#' @param dt_pg_list A list of the PG data.table for all markers
#' @param dt_ref The data.table for the reference profiles
#' @param hyp_id_all Contributors' IDs of all possible hypotheses
#' @param ref_id_lr References' IDs included in any hypotheses
calc_loglike_all_hyp <- function(markers_use, peaks_list, pop_al_list, pop_prob_list, theta, dt_pg_list, dt_ref, hyp_id_all, ref_id_lr){
  n_mk <- length(markers_use)
  n_hyp <- nrow(hyp_id_all)

  tmp <- create_prof_ref_all(dt_ref, markers_use, peaks_list)
  prof_ref_all <- tmp$prof_ref_all
  prof_ref_99_all <- tmp$prof_ref_99_all

  gtc_logprob_all_mk <- vector(mode = "list", length = n_mk)
  for(i in 1:n_mk){
    gtc_logprob_all_mk[[i]] <- calc_gtc_logprob_all_hyp(markers_use[i], hyp_id_all, prof_ref_all, prof_ref_99_all, peaks_list, pop_al_list, pop_prob_list, dt_pg_list, theta, ref_id_lr)
  }

  loglike_all_hyp <- matrix(0, n_mk, n_hyp)
  for(i in 1:n_hyp){
    loglike_all_hyp[, i] <- calc_loglike_r(dt_pg_list, gtc_logprob_all_mk, i)
  }

  return(loglike_all_hyp)
}

#' calc_subsub_log10like_all_mk
#'
#' @description The function to calculate sub-sub-source log10 likelihood
#' @param loglike_all_hyp A matrix of log10 likelihoods of all markers in all possible hypotheses
calc_subsub_log10like_all_mk <- function(loglike_all_hyp){
  overall_loglike_all_hyp <- apply(loglike_all_hyp, 2, sum)
  index_max_loglike <- which.max(overall_loglike_all_hyp)
  subsub_log10like_all_mk <- loglike_all_hyp[, index_max_loglike] / log(10)

  return(list(index_max_loglike = index_max_loglike,
              subsub_log10like_all_mk = subsub_log10like_all_mk))
}

#' create_dt_lr
#'
#' @description The function to create the data.table for LR
#' @param subsub_log10like_all_mk_hp A vector of the sub-sub-source log10 likelihoods of all markers in Hp
#' @param subsub_log10like_all_mk_hd A vector of the sub-sub-source log10 likelihoods of all markers in Hd
#' @param markers_use A vector of the autosomal markers used
create_dt_lr <- function(subsub_log10like_all_mk_hp, subsub_log10like_all_mk_hd, markers_use){

  n_mk <- length(markers_use)
  index_adopt_hp <- which(subsub_log10like_all_mk_hp != 0)
  index_adopt_hd <- which(subsub_log10like_all_mk_hd != 0)
  index_adopt_lr <- sort(intersect(index_adopt_hp, index_adopt_hd))

  log10like_hp_adopt <- log10like_hd_adopt <- rep(NA, n_mk + 1)
  log10like_hp_adopt[index_adopt_lr] <- subsub_log10like_all_mk_hp[index_adopt_lr]
  log10like_hp_adopt[n_mk + 1] <- sum(subsub_log10like_all_mk_hp[index_adopt_lr])
  log10like_hd_adopt[index_adopt_lr] <- subsub_log10like_all_mk_hd[index_adopt_lr]
  log10like_hd_adopt[n_mk + 1] <- sum(subsub_log10like_all_mk_hd[index_adopt_lr])

  dt_lr <- data.table(Marker = c(markers_use, "LR"),
                      Log10_Likelihood_Hp = log10like_hp_adopt,
                      Log10_Likelihood_Hd = log10like_hd_adopt,
                      LR = 10^(log10like_hp_adopt - log10like_hd_adopt))

  return(dt_lr)
}

#' calculate_lr
#'
#' @description The function to calculate a likelihood ratio
#' @param noc_hp The assumed number of contributors in Hp
#' @param hp The contributor names in Hp
#' @param noc_hd The assumed number of contributors in Hd
#' @param hd The contributor names in Hd
#' @param refs The names of all references
#' @param markers_use A vector of the autosomal markers used
#' @param peaks_list A list of peaks
#' @param pop_al_list A list of the alleles in population database
#' @param pop_prob_list A list of the allele probabilities in population database
#' @param theta The theta correction to consider sub-population
#' @param dt_pg_list_hp A list of the PG data.table for all markers in Hp
#' @param dt_pg_list_hd A list of the PG data.table for all markers in Hd
#' @param dt_ref The data.table for the reference profiles
#' @param waiter_ui Information on UI for waiter
calculate_lr <- function(noc_hp, hp, noc_hd, hd, refs,
                         markers_use, peaks_list, pop_al_list, pop_prob_list, theta, dt_pg_list_hp, dt_pg_list_hd, dt_ref,
                         waiter_ui = NULL){

  # --- Contributor IDs for all possible hypotheses --- #

  hyp_id_all_hp <- make_hyp_id_all(noc_hp, hp, refs)
  hyp_id_all_hd <- make_hyp_id_all(noc_hd, hd, refs)
  ref_id_lr <- sort(setdiff(c(as.vector(hyp_id_all_hp), as.vector(hyp_id_all_hd)), 0))

  if(noc_hp == noc_hd && any(is.element(hp, hd))){

    # --- Calculate the likelihood considering assumed contributors --- #

    assumed_contrib <- hp[is.element(hp, hd)]
    hyp_id_all_assumed_contrib <- make_hyp_id_all(noc_hp, assumed_contrib, refs)
    loglike_all_hyp_assumed_contrib <- calc_loglike_all_hyp(markers_use, peaks_list, pop_al_list, pop_prob_list, theta, dt_pg_list_hp, dt_ref, hyp_id_all_assumed_contrib, ref_id_lr)
    tmp <- calc_subsub_log10like_all_mk(loglike_all_hyp_assumed_contrib)
    index_max_loglike_assumed_contrib <- tmp$index_max_loglike
    hyp_id_best_assumed_contrib <- hyp_id_all_assumed_contrib[index_max_loglike_assumed_contrib, ]

    # --- Select possible hypotheses from hyp_id_all_hp and hyp_id_all_hd fixing assumed contributors' position --- #

    hyp_id_all_hp <- select_hyp_id_assumed_contrib(hyp_id_all_hp, hyp_id_best_assumed_contrib)
    hyp_id_all_hd <- select_hyp_id_assumed_contrib(hyp_id_all_hd, hyp_id_best_assumed_contrib)
  }

  # --- Calculate log10 likelihood of Hp --- #

  if(!is.null(waiter_ui)){
    waiter_ui$update(
      html = tagList(
        spin_3k(),
        h3("Calculating log10 likelihood of Hp...")
      )
    )
  }

  loglike_all_hyp_hp <- calc_loglike_all_hyp(markers_use, peaks_list, pop_al_list, pop_prob_list, theta, dt_pg_list_hp, dt_ref, hyp_id_all_hp, ref_id_lr)
  tmp <- calc_subsub_log10like_all_mk(loglike_all_hyp_hp)
  index_max_loglike_hp <- tmp$index_max_loglike
  subsub_log10like_all_mk_hp <- tmp$subsub_log10like_all_mk

  # --- Calculate log10 likelihood of Hd --- #

  if(!is.null(waiter_ui)){
    waiter_ui$update(
      html = tagList(
        spin_3k(),
        h3("Calculating log10 likelihood of Hd...")
      )
    )
  }

  loglike_all_hyp_hd <- calc_loglike_all_hyp(markers_use, peaks_list, pop_al_list, pop_prob_list, theta, dt_pg_list_hd, dt_ref, hyp_id_all_hd, ref_id_lr)
  tmp <- calc_subsub_log10like_all_mk(loglike_all_hyp_hd)
  index_max_loglike_hd <- tmp$index_max_loglike
  subsub_log10like_all_mk_hd <- tmp$subsub_log10like_all_mk

  # --- Calculate LR --- #

  if(!is.null(waiter_ui)){
    waiter_ui$update(
      html = tagList(
        spin_3k(),
        h3("Calculating LR...")
      )
    )
  }

  dt_lr <- create_dt_lr(subsub_log10like_all_mk_hp, subsub_log10like_all_mk_hd, markers_use)
  names(dt_lr) <- c("Marker", "log10 (f(E | Hp))", "log10 (f(E | Hd))", "LR")

  # --- Create output data --- #

  contrib_order_hp <- create_contrib_order(hyp_id_all_hp[index_max_loglike_hp, ], refs)
  contrib_order_hd <- create_contrib_order(hyp_id_all_hd[index_max_loglike_hd, ], refs)
  dt_lr_calc_record_1 <- data.table("Hp" = paste0(contrib_order_hp, collapse = ", "),
                                    "Hd" = paste0(contrib_order_hd, collapse = ", "),
                                    "LR" = dt_lr$LR[length(markers_use) + 1])

  lr_data <- list("hp" = paste0(contrib_order_hp, collapse = ", "),
                  "hd" = paste0(contrib_order_hd, collapse = ", "),
                  "dt_lr" = dt_lr)

  return(list("lr_data" = lr_data, "dt_lr_calc_record_1" = dt_lr_calc_record_1))
}

#' create_contrib_order
#'
#' @description The function to create a vector of the ordered contributor
#' @param hyp_id_best A vector of the contributors' IDs for the best hypothesis
#' @param refs A vector of the references' names
create_contrib_order <- function(hyp_id_best, refs){
  hnc <- length(hyp_id_best)
  contrib_order <- rep("U", hnc)
  for(i in 1:hnc){
    id1 <- hyp_id_best[i]
    if(id1 != 0){
      contrib_order[i] <- refs[id1]
    }
  }
  return(contrib_order)
}
