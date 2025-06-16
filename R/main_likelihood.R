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
      prof_ref[2 * i - 1] <- as.numeric(tmp$Allele1)
      prof_ref[2 * i] <- as.numeric(tmp$Allele2)
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

#' make_n_obs
#'
#' @description The function to make an object for the number of observed alleles in reference profiles
#' @param pop_al A vector of the alleles in the population database
#' @param prof_ref A vector of the reference alleles in one marker
make_n_obs <- function(pop_al, prof_ref){
  prob_al <- sort(unique(c(pop_al, prof_ref)))
  n_obs <- rep(0, length(prob_al))
  names(n_obs) <- prob_al
  if(length(prof_ref) != 0){
    for(i in 1:length(prof_ref)){
      index_al <- which(prob_al == prof_ref[i])
      n_obs[index_al] <- n_obs[index_al] + 1
    }
  }
  return(n_obs)
}

#' calc_log10like_use
#'
#' @description The function to calculate the log10 likelihood of all markers
#' @param dt_pg_list A list of the PG data.table for all markers
#' @param gtc_logprob_all_mk A list of the lob probabilities of each genotype combination for all markers when assuming POI or known contributors
#' @param index_hyp The index of the targeted hypothesis
calc_log10like_use <- function(dt_pg_list, gtc_logprob_all_mk = NULL, index_hyp = NULL){
  n_mk <- length(dt_pg_list)
  log10like_use <- rep(0, n_mk)
  names(log10like_use) <- names(dt_pg_list)

  # log10 likelihood in each marker
  if(is.null(gtc_logprob_all_mk)){
    for(i in 1:n_mk){
      dt_pg <- dt_pg_list[[i]]
      loglike_gtc <- dt_pg$Loglike
      log10like_use[i] <- log10(sum(exp(loglike_gtc)))
    }
  }else{
    n_rep <- length(grep("Loglike_ph_rep", names(dt_pg_list[[1]])))
    for(i in 1:n_mk){
      dt_pg <- dt_pg_list[[i]]
      loglike_gtc <- gtc_logprob_all_mk[[i]][[index_hyp]]

      # With reference genotype(s)
      if(!is.null(loglike_gtc)){
        for(j in 1:n_rep){
          loglike_gtc <- loglike_gtc + dt_pg[[paste0("Loglike_ph_rep", j)]]
        }
        log10like_use[i] <- log10(sum(exp(loglike_gtc)))
      }
    }
  }

  return(log10like_use)
}

#' calc_gtc_logprob_all_hyp
#'
#' @description The function to calculate log-probability of each genotype combination for all hypotheses
#' @param mk The marker name
#' @param hyp_id_all IDs of all possible hypotheses
#' @param prof_ref_all The object is created by the function 'create_prof_ref_all'
#' @param prof_ref_99_all The object is created by the function 'create_prof_ref_all'
#' @param pop_al_list A list of the alleles in population database
#' @param pop_prob_list A list of the allele probabilities in population database
#' @param dt_pg_list A list of the PG data.table for all markers
#' @param theta The theta value to consider sub-population effect
calc_gtc_logprob_all_hyp <- function(mk, hyp_id_all, prof_ref_all, prof_ref_99_all, pop_al_list, pop_prob_list, dt_pg_list, theta){
  n_hyp <- nrow(hyp_id_all)
  gtc_logprob_all_hyp <- vector(mode = "list", length = n_hyp)

  # Check reference profiles of assumed contributors are obtained or not
  prof_ref <- prof_ref_all[mk, ]
  tmp <- setdiff(as.vector(hyp_id_all), 0)
  index_ref_contrib <- sort(c(2 * tmp - 1, 2 * tmp))
  if(all(prof_ref[index_ref_contrib] != -1)){
    # Extract data of one marker
    prof_ref_99 <- prof_ref_99_all[mk, ]
    pop_al <- pop_al_list[[mk]]
    pop_prob <- pop_prob_list[[mk]]
    dt_pg <- dt_pg_list[[mk]]

    # Calculate genotype probabilities
    n_obs <- make_n_obs(pop_al, prof_ref[prof_ref != -1])
    obs_al <- as.numeric(names(n_obs))

    # Make a list of genotype combinations for C++ format
    gtc_mk <- as.matrix(dt_pg[, grep("Allele", names(dt_pg), value = TRUE), with = FALSE])
    gtc_mk <- split(gtc_mk, row(gtc_mk))

    # Calculate genotype probabilities
    for(i in 1:n_hyp){
      gtc_prob_mk <- calc_gtc_prob_mk(gtc_mk, pop_al, pop_prob, theta, prof_ref_99, hyp_id_all[i, ], n_obs, obs_al)
      gtc_logprob_all_hyp[[i]] <- log(gtc_prob_mk)
    }
  }

  return(gtc_logprob_all_hyp)
}

#' log_sum_10
#'
#' @description The function to calculate log10 sum 10^log10_vals
#' @param log10_vals A vector of the log10 values
calc_log_sum_10 <- function(log10_vals){
  max_log10 <- max(log10_vals)
  return(max_log10 + log10(sum(10^(log10_vals - max_log10))))
}

#' calc_log10like_all_hyp
#'
#' @description The function to assign the log10 likelihood
#' @param markers_use A vector of the autosomal markers used
#' @param peaks_list A list of peaks
#' @param pop_al_list A list of the alleles in population database
#' @param pop_prob_list A list of the allele probabilities in population database
#' @param theta The theta correction to consider sub-population
#' @param dt_pg_list A list of the PG data.table for all markers
#' @param dt_ref The data.table for the reference profiles
#' @param hyp_id_all IDs of all possible hypotheses
calc_log10like_all_hyp <- function(markers_use, peaks_list, pop_al_list, pop_prob_list, theta, dt_pg_list, dt_ref, hyp_id_all){
  n_mk <- length(markers_use)
  n_hyp <- nrow(hyp_id_all)

  # Only unknowns
  if(all(as.vector(hyp_id_all) == 0)){
    log10like_all_hyp <- matrix(0, n_mk, 1)
    log10like_all_hyp[, 1] <- calc_log10like_use(dt_pg_list)
  }

  # Including contributors with profiles
  else{
    tmp <- create_prof_ref_all(dt_ref, markers_use, peaks_list)
    prof_ref_all <- tmp$prof_ref_all
    prof_ref_99_all <- tmp$prof_ref_99_all

    gtc_logprob_all_mk <- vector(mode = "list", length = n_mk)
    for(i in 1:n_mk){
      gtc_logprob_all_mk[[i]] <- calc_gtc_logprob_all_hyp(markers_use[i], hyp_id_all, prof_ref_all, prof_ref_99_all, pop_al_list, pop_prob_list, dt_pg_list, theta)
    }

    log10like_all_hyp <- matrix(0, n_mk, n_hyp)
    for(i in 1:n_hyp){
      log10like_all_hyp[, i] <- calc_log10like_use(dt_pg_list, gtc_logprob_all_mk, i)
    }
  }

  return(log10like_all_hyp)
}

#' select_alternative_hyp
#'
#' @description The function to select alternative hypothesis
#' @param hyp_id_best The best contributors' IDs
#' @param hyp_id_all_alternative IDs of all possible alternative hypotheses
#' @param overall_log10like_all_hyp_alternative A vector of overall log10 likelihoods in all possible alternative hypotheses
select_alternative_hyp <- function(hyp_id_best, hyp_id_all_alternative, overall_log10like_all_hyp_alternative){
  known_ids_best <- setdiff(hyp_id_best, 0)
  ids_share <- known_ids_best[is.element(known_ids_best, setdiff(hyp_id_all_alternative[1, ], 0))]
  n_ids_share <- length(ids_share)

  if(n_ids_share != 0){
    bool_mat <- matrix(FALSE, nrow(hyp_id_all_alternative), n_ids_share)

    for(i in 1:n_ids_share){
      id1 <- ids_share[i]
      index_id1 <- which(hyp_id_best == id1)
      bool_mat[, i] <- hyp_id_all_alternative[, index_id1] == id1
    }

    best_log10like_alternative <- max(overall_log10like_all_hyp_alternative[apply(bool_mat, 1, all)])
    index_alternative <- which(overall_log10like_all_hyp_alternative == best_log10like_alternative)

  }else{
    index_alternative <- integer(0)
  }

  return(index_alternative)
}

#' create_dt_lr
#'
#' @description The function to create the data.table for LR
#' @param log10like_all_hyp_hp A matrix of log10 likelihoods of all markers in all possible Hp
#' @param log10like_all_hyp_hd A matrix of log10 likelihoods of all markers in all possible Hd
#' @param hyp_id_all_hp IDs of all possible hypotheses in Hp
#' @param hyp_id_all_hd IDs of all possible hypotheses in Hd
#' @param markers_use A vector of the autosomal markers used
create_dt_lr <- function(log10like_all_hyp_hp, log10like_all_hyp_hd, hyp_id_all_hp, hyp_id_all_hd, markers_use){

  # --- Calculate sub-sub-source log10 likelihood --- #

  noc_hp <- ncol(hyp_id_all_hp)
  known_ids_hp <- sort(setdiff(hyp_id_all_hp[1, ], 0))

  noc_hd <- ncol(hyp_id_all_hd)
  known_ids_hd <- sort(setdiff(hyp_id_all_hd[1, ], 0))

  overall_log10like_all_hyp_hp <- apply(log10like_all_hyp_hp, 2, sum)
  overall_log10like_all_hyp_hd <- apply(log10like_all_hyp_hd, 2, sum)

  index_max_log10like_hp <- which.max(overall_log10like_all_hyp_hp)
  max_log10like_hp <- overall_log10like_all_hyp_hp[index_max_log10like_hp]
  index_max_log10like_hd <- which.max(overall_log10like_all_hyp_hd)
  max_log10like_hd <- overall_log10like_all_hyp_hd[index_max_log10like_hd]

  # Select the best hypothesis
  # - NOC in Hp is the same as that in Hd, and
  # - There is at least one assumed contributor in each hypothesis
  if(noc_hp == noc_hd && any(is.element(known_ids_hp, known_ids_hd))){

    if(max_log10like_hp > max_log10like_hd){

      # Select alternative hypothesis
      hyp_id_best <- hyp_id_all_hp[index_max_log10like_hp, ]
      hyp_id_all_alternative <- hyp_id_all_hd
      overall_log10like_all_hyp_alternative <- overall_log10like_all_hyp_hd
      index_alternative <- select_alternative_hyp(hyp_id_best, hyp_id_all_alternative, overall_log10like_all_hyp_alternative)

      # Select Hp based on the maximum likelihood
      subsub_log10like_all_mk_hp <- log10like_all_hyp_hp[, index_max_log10like_hp]
      subsub_log10like_hp <- max_log10like_hp

      # Select Hd based on the selected alternative hypothesis
      subsub_log10like_all_mk_hd <- log10like_all_hyp_hd[, index_alternative]
      subsub_log10like_hd <- overall_log10like_all_hyp_hd[index_alternative]
    }else{

      # Select alternative hypothesis
      hyp_id_best <- hyp_id_all_hd[index_max_log10like_hd, ]
      hyp_id_all_alternative <- hyp_id_all_hp
      overall_log10like_all_hyp_alternative <- overall_log10like_all_hyp_hp
      index_alternative <- select_alternative_hyp(hyp_id_best, hyp_id_all_alternative, overall_log10like_all_hyp_alternative)

      # Select Hp based on the selected alternative hypothesis
      subsub_log10like_all_mk_hp <- log10like_all_hyp_hp[, index_alternative]
      subsub_log10like_hp <- overall_log10like_all_hyp_hp[index_alternative]

      # Select Hd based on the maximum likelihood
      subsub_log10like_all_mk_hd <- log10like_all_hyp_hd[, index_max_log10like_hd]
      subsub_log10like_hd <- max_log10like_hd
    }
  }

  # Select the best hypothesis
  # - NOC in Hp is not the same as that in Hd, or
  # - There is no assumed contributor in each hypothesis
  else{

    # Select Hp based on the maximum likelihood
    subsub_log10like_all_mk_hp <- log10like_all_hyp_hp[, index_max_log10like_hp]
    subsub_log10like_hp <- max_log10like_hp

    # Select Hd based on the maximum likelihood
    subsub_log10like_all_mk_hd <- log10like_all_hyp_hd[, index_max_log10like_hd]
    subsub_log10like_hd <- max_log10like_hd
  }

  # --- Calculate sub-source log10 likelihood --- #

  n_unk_hp <- length(which(is.element(hyp_id_all_hp[1, ], 0)))
  sub_log10like_hp <- log10(factorial(n_unk_hp)) + calc_log_sum_10(overall_log10like_all_hyp_hp)

  n_unk_hd <- length(which(is.element(hyp_id_all_hd[1, ], 0)))
  sub_log10like_hd <- log10(factorial(n_unk_hd)) + calc_log_sum_10(overall_log10like_all_hyp_hd)

  # --- Create a data.table for LR --- #

  n_mk <- length(markers_use)
  index_adopt_hp <- which(subsub_log10like_all_mk_hp != 0)
  index_adopt_hd <- which(subsub_log10like_all_mk_hd != 0)
  index_adopt_lr <- sort(intersect(index_adopt_hp, index_adopt_hd))

  log10like_hp_adopt <- log10like_hd_adopt <- rep(NA, n_mk + 2)
  log10like_hp_adopt[index_adopt_lr] <- subsub_log10like_all_mk_hp[index_adopt_lr]
  log10like_hp_adopt[n_mk + 1] <- subsub_log10like_hp
  log10like_hp_adopt[n_mk + 2] <- sub_log10like_hp
  log10like_hd_adopt[index_adopt_lr] <- subsub_log10like_all_mk_hd[index_adopt_lr]
  log10like_hd_adopt[n_mk + 1] <- subsub_log10like_hd
  log10like_hd_adopt[n_mk + 2] <- sub_log10like_hd

  dt_lr <- data.table(Marker = c(markers_use, "Sub-Sub-Source LR", "Sub-Source LR"),
                      Log10_Likelihood_Hp = log10like_hp_adopt,
                      Log10_Likelihood_Hd = log10like_hd_adopt,
                      LR = 10^(log10like_hp_adopt - log10like_hd_adopt))

  return(dt_lr)
}
