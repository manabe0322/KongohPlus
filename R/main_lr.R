#' make_prof_ref_all
#'
#' @description The function to make a matrix for reference profiles of all autosomal markers
#' @param dt_ref The data.table of the reference profile
#' @param markers_auto A vector of the autosomal markers
#' @param peaks_list A list of peaks
make_prof_ref_all <- function(dt_ref, markers_auto, peaks_list){
  n_mk <- length(markers_auto)
  refs <- unique(dt_ref[, SampleName])
  n_ref <- length(refs)
  setindex(dt_ref, Marker)

  prof_ref_all <- matrix(0, n_mk, 2 * n_ref)
  for(i in 1:n_mk){
    mk <- markers_auto[i]
    dt_ref_one_mk <- dt_ref[Marker == mk]
    peaks <- peaks_list[[mk]]
    for(j in 1:n_ref){
      tmp <- dt_ref_one_mk[j, ]
      al_1 <- as.numeric(tmp[, Allele1])
      if(!is.element(al_1, peaks)){
        al_1 <- 99
      }
      al_2 <- as.numeric(tmp[, Allele2])
      if(!is.element(al_2, peaks)){
        al_2 <- 99
      }
      prof_ref_all[i, c(2 * j - 1, 2 * j)] <- sort(c(al_1, al_2))
    }
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

  prof_ref_all <- make_prof_ref_all(dt_ref, markers_use, peaks_list)

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
    pop_al <- pop_al_list[[mk]]
    pop_prob <- pop_prob_list[[mk]]
    gtc_mk <- gtc_use[[mk]]

    # Calculate genotype probabilities
    n_obs <- make_n_obs(pop_al, prof_ref)
    for(j in 1:n_hyp){
      gtc_prob_mk <- calc_gtc_prob_mk(gtc_mk, pop_al, pop_prob, theta, prof_ref, hyp_id_all[j, ], n_obs)
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
