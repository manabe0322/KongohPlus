#' generate_randgt
#'
#' @description The function to generate random genotypes
#' @param np The number of persons generated
#' @param pop_al_list A list of the alleles in population database
#' @param pop_prob_list A list of the allele probabilities in population database
generate_randgt <- function(np, pop_al_list, pop_prob_list){
  mapply(sample, sapply(pop_al_list, as.numeric), 2 * np, prob = pop_prob_list, replace = TRUE)
}

#' make_hyp_id_randgt
#'
#' @description The function to make ids of the hypothesis for the random genotype
#' @param hnc The hypothesized number of contributors
#' @param hyp Known contributors in a hypothesis
#' @param refs The names of all references
make_hyp_id_randgt <- function(hnc, hyp, refs){
  known_id <- which(is.element(refs, hyp))
  n_known <- length(known_id)
  hyp_id <- rep(0, hnc)
  if(n_known != 0){
    hyp_id[1] <- 1
  }
  return(hyp_id)
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
#' @param dt_mpdeg The data.table for mixture proportions and degradation parameters of each contributor
#' @param dt_ref The data.table for the reference profiles
#' @param hyp Known contributors in a hypothesis
#' @param refs The names of all references
assign_log10like_randgt <- function(hnc, markers_use, peaks_list, pop_al_list, pop_prob_list, theta, gtc_use, lpd_gtc_use_1st, lpd_gtc_use_2nd, dt_ref, hyp, refs){
  n_mk <- length(markers_use)
  prof_ref_all <- make_prof_ref_all(dt_ref, markers_use, peaks_list)
  hyp_id <- make_hyp_id_randgt(hnc, hyp, refs)

  gtc_logprob_use <- vector(mode = "list", length = n_mk)
  for(i in 1:n_mk){
    # Extract data of one marker
    mk <- markers_use[i]
    prof_ref <- prof_ref_all[mk, ]
    pop_al <- pop_al_list[[mk]]
    pop_prob <- pop_prob_list[[mk]]
    gtc_mk <- gtc_use[[mk]]

    # Calculate genotype probabilities
    n_obs <- make_n_obs(pop_al, prof_ref)
    gtc_prob_mk <- calc_gtc_prob_mk(gtc_mk, pop_al, pop_prob, theta, prof_ref, hyp_id, n_obs)
    gtc_logprob_use[[i]] <- log(gtc_prob_mk)
  }

  log10like_adopt <- calc_log10like_use(gtc_logprob_use, lpd_gtc_use_1st, lpd_gtc_use_2nd)
  combined_log10like_adopt <- sum(log10like_adopt)

  return(list("log10like_adopt" = log10like_adopt,
              "combined_log10like_adopt" = combined_log10like_adopt))
}
