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

#' assign_log10like
#'
#' @description The function to assign the log10 likelihood
#' @param hnc The hypothesized number of contributors
#' @param markers_use A vector of the autosomal markers used
#' @param peaks_list A list of peaks
#' @param pop_al_list A list of the alleles in population database
#' @param pop_prob_list A list of the allele probabilities in population database
#' @param theta The theta correction to consider sub-population
#' @param dt_pg_list A list of the PG data.table for all markers
#' @param dt_ref The data.table for the reference profiles
#' @param hyp Known contributors in a hypothesis
#' @param refs The names of all references
assign_log10like <- function(hnc, markers_use, peaks_list, pop_al_list, pop_prob_list, theta, dt_pg_list, dt_ref, hyp, refs){
  n_mk <- length(markers_use)
  hyp_id_all <- make_hyp_id_all(hnc, hyp, refs)
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

  overall_log10like <- apply(log10like_all_hyp, 2, sum)

  # --- sub-sub-source log10 likelihood --- #

  index_max_like <- which.max(overall_log10like)
  subsub_source_log10like_mk <- log10like_all_hyp[, index_max_like]
  subsub_source_log10like <- overall_log10like[index_max_like]

  # --- sub-source log10 likelihood --- #

  n_unk <- length(which(is.element(hyp_id_all[1, ], 0)))
  sub_source_log10like <- log10(factorial(n_unk)) + calc_log_sum_10(overall_log10like)

  log10like_adopt <- c(subsub_source_log10like_mk, subsub_source_log10like, sub_source_log10like)
  names(log10like_adopt) <- c(markers_use, "sub-sub-source", "sub-source")

  return(log10like_adopt)
}

#' create_dt_lr
#'
#' @description The function to create the data.table for LR
#' @param markers_use A vector of the autosomal markers used
#' @param log10like_hp Log 10 likelihood in each marker, sub-sub-source, and sub-source under Hp
#' @param log10like_hd Log 10 likelihood in each marker, sub-sub-source, and sub-source under Hd
create_dt_lr <- function(markers_use, log10like_hp, log10like_hd){
  n_mk <- length(markers_use)
  index_adopt_hp <- which(log10like_hp[1:n_mk] != 0)
  index_adopt_hd <- which(log10like_hd[1:n_mk] != 0)
  index_adopt_lr <- sort(intersect(index_adopt_hp, index_adopt_hd))

  log10like_hp_adopt <- log10like_hd_adopt <- rep(NA, n_mk + 2)
  log10like_hp_adopt[index_adopt_lr] <- log10like_hp[index_adopt_lr]
  log10like_hp_adopt[n_mk + 1] <- log10like_hp[n_mk + 1]
  log10like_hp_adopt[n_mk + 2] <- log10like_hp[n_mk + 2]
  log10like_hd_adopt[index_adopt_lr] <- log10like_hd[index_adopt_lr]
  log10like_hd_adopt[n_mk + 1] <- log10like_hd[n_mk + 1]
  log10like_hd_adopt[n_mk + 2] <- log10like_hd[n_mk + 2]

  dt_lr <- data.table(Marker = c(markers_use, "Sub-Sub-Source LR", "Sub-Source LR"),
                      Log10_Likelihood_Hp = log10like_hp_adopt,
                      Log10_Likelihood_Hd = log10like_hd_adopt,
                      LR = 10^(log10like_hp_adopt - log10like_hd_adopt))
  return(dt_lr)
}
