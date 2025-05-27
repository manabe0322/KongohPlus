#' make_dt_pg_list
#'
#' @description The function to make a list of data.tables for the result of probabilistic genotyping
#' @param markers_use A character vector of the used markers
#' @param label_contrib A character vector of the contributor names
#' @param gtc_use Genotype combinations in the used markers
#' @param gtc_logprob_use A list of the genotype probabilities in each used marker
#' @param lpd_gtc_use A list of the log likelihood of each genotype in all markers and replicates
make_dt_pg_list <- function(markers_use, label_contrib, gtc_use, gtc_logprob_use, lpd_gtc_use){
  hnc <- length(label_contrib)
  n_mk <- length(gtc_use)
  n_rep <- length(lpd_gtc_use)
  dt_pg_list <- vector(mode = "list", length = n_mk)

  for(i in 1:n_mk){

    # Paste two alleles of each contributor and reverse the contributor order (minor to major -> major to minor)
    gtc_mk <- matrix(unlist(gtc_use[[i]]), ncol = 2 * hnc, byrow = TRUE)
    gtc_mk[which(gtc_mk == 99, arr.ind = TRUE)] <- "Q"
    gtc_mk_new <- matrix("", nrow(gtc_mk), hnc)
    for(j in 1:hnc){
      gtc_mk_new[, j] <- apply(gtc_mk[, c(2 * (hnc - j) + 1, 2 * (hnc - j) + 2)], 1, paste, collapse = ", ")
    }

    # Calculate posterior probabilities of each genotype combination
    gtc_logprob_mk <- gtc_logprob_use[[i]]
    loglike_gtc <- gtc_logprob_mk
    for(j in 1:n_rep){
      loglike_gtc <- loglike_gtc + lpd_gtc_use[[j]][[i]]
    }
    tmp <- exp(loglike_gtc)
    post_prob_gtc <- tmp / sum(tmp)

    # Make a data.table
    mat_pg <- cbind(gtc_mk_new, gtc_logprob_mk)
    for(j in 1:n_rep){
      mat_pg <- cbind(mat_pg, lpd_gtc_use[[j]][[i]])
    }
    mat_pg <- cbind(mat_pg, loglike_gtc, post_prob_gtc)
    mat_pg <- mat_pg[post_prob_gtc != 0, , drop = FALSE]
    dt_pg <- as.data.table(mat_pg)
    names(dt_pg) <- c(label_contrib, "Logprob_gtc", paste0("Loglike_ph_rep", 1:n_rep), "Loglike", "Probability")

    # Change data type
    col_numeric <- c("Logprob_gtc", paste0("Loglike_ph_rep", 1:n_rep), "Loglike", "Probability")
    dt_pg[, (col_numeric) := lapply(.SD, as.numeric), .SDcols = col_numeric]

    # Descending order for posterior probabilities
    dt_pg <- dt_pg[order(dt_pg[, Probability], decreasing = TRUE)]

    dt_pg_list[[i]] <- dt_pg
  }
  names(dt_pg_list) <- markers_use
  return(dt_pg_list)
}

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
#' @param gtc_logprob_list A list of the lob probabilities of each genotype combination for all markers when assuming POI or known contributors
calc_log10like_use <- function(dt_pg_list, gtc_logprob_list = NULL){
  n_mk <- length(dt_pg_list)
  log10like_use <- rep(0, n_mk)
  names(log10like_use) <- names(dt_pg_list)

  # log10 likelihood in each marker
  if(is.null(gtc_logprob_list)){
    for(i in 1:n_mk){
      dt_pg <- dt_pg_list[[i]]
      loglike_gtc <- dt_pg$Loglike
      log10like_use[i] <- log10(sum(exp(loglike_gtc)))
    }
  }else{
    n_rep <- length(grep("Loglike_ph_rep", names(dt_pg_list[[1]])))
    for(i in 1:n_mk){
      dt_pg <- dt_pg_list[[i]]
      loglike_gtc <- gtc_logprob_list[[i]]

      # With reference genotype(s)
      if(!is.null(loglike_gtc)){
        for(j in 1:n_rep){
          loglike_gtc <- loglike_gtc + dt_pg[[paste0("Loglike_ph_rep", j)]]
        }
        log10like_use[i] <- log10(sum(exp(loglike_gtc)))
      }

      # Without reference genotype(s)
      else{
        log10like_use[i] <- NA
      }
    }
  }

  return(log10like_use)
}

#' calc_sub_source_like
#'
#' @description The function to calculate sub-source likelihood
#' @param hyp_id_all IDs of all possible hypotheses
#' @param log10like_all_hyp Log 10 likelihoods in each marker for all hypotheses
calc_sub_source_like <- function(hyp_id_all, log10like_all_hyp){
  n_unk <- length(which(is.element(hyp_id_all[1, ], 0)))
  log10like_adopt <- log10(factorial(n_unk) * apply(10^log10like_all_hyp, 1, sum))
  return(log10like_adopt)
}

#' calc_subsub_source_like
#'
#' @description The function to calculate sub-sub-source likelihood
#' @param log10like_all_hyp Log 10 likelihoods for all hypotheses
calc_subsub_source_like <- function(log10like_all_hyp){
  overall_like <- apply(log10like_all_hyp, 2, sum)
  index_max_like <- which.max(overall_like)
  log10like_adopt <- log10like_all_hyp[, index_max_like]
  return(log10like_adopt)
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
#' @param sub_source Whether the sub-source level likelihood is calculated or not
assign_log10like <- function(hnc, markers_use, peaks_list, pop_al_list, pop_prob_list, theta, dt_pg_list, dt_ref, hyp, refs, sub_source = TRUE){
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

    gtc_logprob_all_hyp <- vector(mode = "list", length = n_hyp)
    for(i in 1:n_hyp){
      gtc_logprob_all_hyp[[i]] <- vector(mode = "list", length = n_mk)
    }

    for(i in 1:n_mk){

      # Check reference profiles of assumed contributors are obtained or not
      mk <- markers_use[i]
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
        gtc_mk <- matrix("", nrow(dt_pg), 2 * hnc)
        for(j in 1:hnc){
          gtc_mk[, c(2 * j - 1, 2 * j)] <- matrix(unlist(strsplit(dt_pg[[paste0("U", j)]], ", ")), ncol = 2, byrow = TRUE)
        }
        gtc_mk[which(gtc_mk == "Q", arr.ind = TRUE)] <- "99"
        gtc_mk <- split(gtc_mk, row(gtc_mk))
        gtc_mk <- lapply(gtc_mk, as.numeric)

        # Calculate genotype probabilities
        for(j in 1:n_hyp){
          gtc_prob_mk <- calc_gtc_prob_mk(gtc_mk, pop_al, pop_prob, theta, prof_ref_99, hyp_id_all[j, ], n_obs, obs_al)
          gtc_logprob_all_hyp[[j]][[i]] <- log(gtc_prob_mk)
        }
      }
    }

    log10like_all_hyp <- matrix(0, n_mk, n_hyp)
    for(i in 1:n_hyp){
      log10like_all_hyp[, i] <- calc_log10like_use(dt_pg_list, gtc_logprob_all_hyp[[i]])
    }
  }

  if(sub_source){
    log10like_adopt <- calc_sub_source_like(hyp_id_all, log10like_all_hyp)
  }else{
    log10like_adopt <- calc_subsub_source_like(log10like_all_hyp)
  }
  names(log10like_adopt) <- markers_use

  return(log10like_adopt)
}

#' create_dt_lr
#'
#' @description The function to create the data.table for LR
#' @param markers_use A vector of the autosomal markers used
#' @param log10like_hp Log 10 likelihood in each marker under Hp
#' @param log10like_hd Log 10 likelihood in each marker under Hd
create_dt_lr <- function(markers_use, log10like_hp, log10like_hd){
  n_mk <- length(markers_use)
  index_adopt_hp <- which(!is.na(log10like_hp))
  index_adopt_hd <- which(!is.na(log10like_hd))
  index_adopt_lr <- sort(intersect(index_adopt_hp, index_adopt_hd))

  log10like_hp_adopt <- log10like_hd_adopt <- rep(NA, n_mk + 1)
  tmp_hp <- log10like_hp[index_adopt_lr]
  log10like_hp_adopt[index_adopt_lr] <- tmp_hp
  log10like_hp_adopt[n_mk + 1] <- sum(tmp_hp)
  tmp_hd <- log10like_hd[index_adopt_lr]
  log10like_hd_adopt[index_adopt_lr] <- tmp_hd
  log10like_hd_adopt[n_mk + 1] <- sum(tmp_hd)

  dt_lr <- data.table(Marker = c(markers_use, "Overall"),
                      Log10_Likelihood_Hp = log10like_hp_adopt,
                      Log10_Likelihood_Hd = log10like_hd_adopt,
                      LR = 10^(log10like_hp_adopt - log10like_hd_adopt))
  return(dt_lr)
}
