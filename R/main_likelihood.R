#' extract_prof_ref
#'
#' @description The function to extract reference profiles of one marker
#' @param dt_ref The data.table of the reference profile
#' @param mk The target marker
extract_prof_ref <- function(dt_ref, mk){
  dt_ref_one_mk <- dt_ref[dt_ref$Marker == mk]
  n_ref <- nrow(dt_ref_one_mk)
  prof_ref <- rep(0, 2 * n_ref)
  for(i in 1:n_ref){
    tmp <- dt_ref_one_mk[i, ]
    prof_ref[2 * i - 1] <- as.numeric(tmp$Allele1)
    prof_ref[2 * i] <- as.numeric(tmp$Allele2)
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
  refs <- unique(dt_ref[, SampleName])
  prof_ref_all <- matrix(0, n_mk, 2 * length(refs))
  rownames(prof_ref_all) <- markers_auto
  colnames(prof_ref_all) <- as.vector(sapply(refs, rep, 2))
  prof_ref_99_all <- prof_ref_all
  for(i in 1:n_mk){
    mk <- markers_auto[i]
    prof_ref_all[i, ] <- extract_prof_ref(dt_ref, mk)
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

#' calc_loglike_gtc
#'
#' @description The function to calculate the log likelihood of each genotype combination in one marker
#' @param gtc_logprob_mk A vector of the genotype probabilities in one marker
#' @param lpd_gtc_mk_1st A vector of the log likelihood of each genotype in one marker (replicate 1)
#' @param lpd_gtc_mk_2nd A vector of the log likelihood of each genotype in one marker (replicate 2)
calc_loglike_gtc <- function(gtc_logprob_mk, lpd_gtc_mk_1st, lpd_gtc_mk_2nd = NULL){
  if(is.null(lpd_gtc_mk_2nd)){
    loglike_gtc <- gtc_logprob_mk + lpd_gtc_mk_1st
  }else{
    loglike_gtc <- gtc_logprob_mk + lpd_gtc_mk_1st + lpd_gtc_mk_2nd
  }
  return(loglike_gtc)
}

#' calc_log10like_use
#'
#' @description The function to calculate the log10 likelihood of the used markers
#' @param gtc_logprob_use A list of the genotype probabilities in each used marker
#' @param lpd_gtc_use_1st A list of the log likelihood of each genotype in each used marker (replicate 1)
#' @param lpd_gtc_use_2nd A list of the log likelihood of each genotype in each used marker (replicate 2)
calc_log10like_use <- function(gtc_logprob_use, lpd_gtc_use_1st, lpd_gtc_use_2nd = NULL){
  n_mk <- length(lpd_gtc_use_1st)
  log10like_use <- rep(0, n_mk)

  if(is.null(lpd_gtc_use_2nd)){
    lpd_gtc_use_2nd <- vector(mode = "list", length = n_mk)
  }

  for(i in 1:n_mk){
    loglike_gtc <- calc_loglike_gtc(gtc_logprob_use[[i]], lpd_gtc_use_1st[[i]], lpd_gtc_use_2nd[[i]])
    log10like_use[i] <- log10(sum(exp(loglike_gtc)))
  }
  return(log10like_use)
}

#' make_dt_pg
#'
#' @description The function to make a data.table for the result of probabilistic genotyping in one marker
#' @param label_contrib A character vector of the contributor names
#' @param gtc_mk Genotype combinations in one marker
#' @param gtc_logprob_mk A vector of the genotype probabilities in one marker
#' @param lpd_gtc_mk_1st A vector of the log likelihood of each genotype in one marker (replicate 1)
#' @param lpd_gtc_mk_2nd A vector of the log likelihood of each genotype in one marker (replicate 2)
make_dt_pg <- function(label_contrib, gtc_mk, gtc_logprob_mk, lpd_gtc_mk_1st, lpd_gtc_mk_2nd = NULL){
  hnc <- length(gtc_mk[[1]]) / 2
  gtc_mk <- matrix(unlist(gtc_mk), ncol = 2 * hnc, byrow = TRUE)
  gtc_mk[which(gtc_mk == 99, arr.ind = TRUE)] <- "Q"
  gtc_mk_new <- matrix("", nrow(gtc_mk), hnc)

  # Paste two alleles reversing sequence (minor to major -> major to minor)
  for(i in 1:hnc){
    gtc_mk_new[, i] <- apply(gtc_mk[, c(2 * (hnc - i) + 1, 2 * (hnc - i) + 2)], 1, paste, collapse = ", ")
  }

  # Assign NA to lpd_gtc_mk_2nd when performing one PCR
  if(is.null(lpd_gtc_mk_2nd)){
    lpd_gtc_mk_2nd <- rep(NA, length(lpd_gtc_mk_1st))
  }

  # Cut data with genotype probability = 0
  index_keep <- which(gtc_logprob_mk != -Inf)
  gtc_mk_new <- gtc_mk_new[index_keep, , drop = FALSE]
  gtc_logprob_mk <- gtc_logprob_mk[index_keep]
  lpd_gtc_mk_1st <- lpd_gtc_mk_1st[index_keep]
  lpd_gtc_mk_2nd <- lpd_gtc_mk_2nd[index_keep]

  # Calculate posterior probabilities of each genotype combination
  if(any(is.na(lpd_gtc_mk_2nd))){
    loglike_gtc <- calc_loglike_gtc(gtc_logprob_mk, lpd_gtc_mk_1st, NULL)
  }else{
    loglike_gtc <- calc_loglike_gtc(gtc_logprob_mk, lpd_gtc_mk_1st, lpd_gtc_mk_2nd)
  }
  tmp <- exp(loglike_gtc)
  post_prob_gtc <- tmp / sum(tmp)

  # Make a data.table
  mat_pg <- cbind(gtc_mk_new, gtc_logprob_mk, lpd_gtc_mk_1st, lpd_gtc_mk_2nd, post_prob_gtc)
  mat_pg <- mat_pg[gtc_logprob_mk != -Inf, , drop = FALSE]
  dt_pg <- as.data.table(mat_pg)
  names(dt_pg) <- c(label_contrib, "Log_gtc_prob", "Loglike_gtc_1st", "Loglike_gtc_2nd", "Probability")

  # Change data type
  col_numeric <- c("Log_gtc_prob", "Loglike_gtc_1st", "Loglike_gtc_2nd", "Probability")
  dt_pg[, (col_numeric) := lapply(.SD, as.numeric), .SDcols = col_numeric]

  # Descending order for weights
  dt_pg <- dt_pg[order(dt_pg[, Probability], decreasing = TRUE)]

  return(dt_pg)
}

#' make_dt_pg_list
#'
#' @description The function to make a list of data.tables for the result of probabilistic genotyping
#' @param markers_use A character vector of the used markers
#' @param label_contrib A character vector of the contributor names
#' @param gtc_use Genotype combinations in the used markers
#' @param gtc_logprob_use A list of the genotype probabilities in each used marker
#' @param lpd_gtc_use_1st A list of the log likelihood of each genotype in each used marker (replicate 1)
#' @param lpd_gtc_use_2nd A list of the log likelihood of each genotype in each used marker (replicate 2)
make_dt_pg_list <- function(markers_use, label_contrib, gtc_use, gtc_logprob_use, lpd_gtc_use_1st, lpd_gtc_use_2nd = NULL){
  n_mk <- length(gtc_use)
  dt_pg_list <- vector(mode = "list", length = n_mk)

  if(is.null(lpd_gtc_use_2nd)){
    lpd_gtc_use_2nd <- vector(mode = "list", length = n_mk)
  }

  for(i in 1:n_mk){
    dt_pg_list[[i]] <- make_dt_pg(label_contrib, gtc_use[[i]], gtc_logprob_use[[i]], lpd_gtc_use_1st[[i]], lpd_gtc_use_2nd[[i]])
  }
  names(dt_pg_list) <- markers_use
  return(dt_pg_list)
}

calc_sub_source_like <- function(hyp_id_all, log10like_tmp){
  n_mk <- nrow(log10like_tmp) - 1
  n_unk <- length(which(is.element(hyp_id_all[1, ], 0)))
  log10like_adopt <- log10(factorial(n_unk) * apply(10^log10like_tmp, 1, sum))
  log10like_adopt[n_mk + 1] <- sum(log10like_adopt[1:n_mk])
  return(log10like_adopt)
}

calc_subsub_source_like <- function(log10like_tmp){
  index_max_like <- which.max(log10like_tmp[nrow(log10like_tmp), ])
  log10like_adopt <- log10like_tmp[, index_max_like]
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
#' @param gtc_use A list of the genotype combinations in used markers
#' @param lpd_gtc_use_1st A list of the log likelihood of each genotype in each used marker (replicate 1)
#' @param lpd_gtc_use_2nd A list of the log likelihood of each genotype in each used marker (replicate 2)
#' @param dt_ref The data.table for the reference profiles
#' @param hyp Known contributors in a hypothesis
#' @param refs The names of all references
#' @param sub_source Whether the sub-source level likelihood is calculated or not
assign_log10like <- function(hnc, markers_use, peaks_list, pop_al_list, pop_prob_list, theta, gtc_use, lpd_gtc_use_1st, lpd_gtc_use_2nd, dt_ref, hyp, refs, sub_source = TRUE){
  n_mk <- length(markers_use)

  tmp <- create_prof_ref_all(dt_ref, markers_use, peaks_list)
  prof_ref_all <- tmp$prof_ref_all
  prof_ref_99_all <- tmp$prof_ref_99_all

  hyp_id_all <- make_hyp_id_all(hnc, hyp, refs)
  n_hyp <- nrow(hyp_id_all)

  gtc_logprob_all_hyp <- vector(mode = "list", length = n_hyp)
  for(i in 1:n_hyp){
    gtc_logprob_all_hyp[[i]] <- vector(mode = "list", length = n_mk)
  }

  for(i in 1:n_mk){
    # Extract data of one marker
    mk <- markers_use[i]
    prof_ref <- prof_ref_all[mk, ]
    prof_ref_99 <- prof_ref_99_all[mk, ]

    pop_al <- pop_al_list[[mk]]
    pop_prob <- pop_prob_list[[mk]]
    gtc_mk <- gtc_use[[mk]]

    # Calculate genotype probabilities
    n_obs <- make_n_obs(pop_al, prof_ref)
    obs_al <- as.numeric(names(n_obs))
    for(j in 1:n_hyp){
      gtc_prob_mk <- calc_gtc_prob_mk(gtc_mk, pop_al, pop_prob, theta, prof_ref_99, hyp_id_all[j, ], n_obs, obs_al)
      gtc_logprob_all_hyp[[j]][[i]] <- log(gtc_prob_mk)
    }
  }

  log10like_tmp <- matrix(0, n_mk + 1, n_hyp)
  for(i in 1:n_hyp){
    log10like_tmp[1:n_mk, i] <- log10like_use <- calc_log10like_use(gtc_logprob_all_hyp[[i]], lpd_gtc_use_1st, lpd_gtc_use_2nd)
    log10like_tmp[n_mk + 1, i] <- sum(log10like_use)
  }

  if(sub_source){
    log10like_adopt <- calc_sub_source_like(hyp_id_all, log10like_tmp)
  }else{
    log10like_adopt <- calc_subsub_source_like(log10like_tmp)
  }
  names(log10like_adopt) <- c(markers_use, "Total")

  return(log10like_adopt)
}
