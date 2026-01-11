#' correct_af_dirichlet
#'
#' @description The function to correct allele frequencies based on the Dirichlet distribution
#' @param pop_al A vector of the alleles in the population database
#' @param pop_freq A vector of the allele frequencies in the population database
correct_af_dirichlet <- function(pop_al, pop_freq){

  # Add a dimension for an unobserved allele named allele '0'
  pop_al <- c(pop_al, 0)
  pop_freq <- c(pop_freq, 0)

  pop_prob <- (pop_freq + 1) / (sum(pop_freq) + length(pop_al))

  order_al <- order(pop_al)
  pop_al <- pop_al[order_al]
  pop_prob <- pop_prob[order_al]

  return(list("pop_al" = pop_al,
              "pop_prob" = pop_prob))
}

#' calc_al_prob_drop
#'
#' @description The function to calculate an allele probability for drop-out
#' @param al_csp A vector of alleles detected in crime stain profiles
#' @param pop_al A vector of the alleles in the population database
#' @param pop_prob A vector of the allele probabilities
calc_al_prob_drop <- function(al_csp, pop_al, pop_prob){
  index_drop <- which(!is.element(pop_al, al_csp))
  if(length(index_drop) != 0){
    pop_prob_drop <- sum(pop_prob[index_drop])
  }else{
    pop_prob_drop <- 0
  }
  return(pop_prob_drop)
}

#' set_al_prob
#'
#' @description The function to set allele probabilities
#' @param markers_use A vector of the autosomal markers used
#' @param dt_af The data.table of the allele frequencies
#' @param peaks_list The peak list of all markers
set_al_prob <- function(markers_use, dt_af, peaks_list){
  n_mk <- length(markers_use)
  pop_al_all <- dt_af[, Allele]
  pop_al_list <- pop_prob_list <- vector(mode = "list", length = n_mk)

  for(i in 1:n_mk){
    mk <- markers_use[i]
    pop_freq <- as.numeric(unlist(dt_af[, mk, with = FALSE]))
    pop_al <- pop_al_all[!is.na(pop_freq)]
    pop_freq <- pop_freq[!is.na(pop_freq)]
    peaks <- peaks_list[[mk]]
    al_csp <- setdiff(peaks, 99)

    unobs_al <- setdiff(al_csp, pop_al)

    # --- Calculate the expected allele probabilities based on the Dirichlet distribution --- #

    tmp <- correct_af_dirichlet(pop_al, pop_freq, unobs_al)
    pop_al <- tmp$pop_al
    pop_prob <- tmp$pop_prob

    # --- Calculate the allele probability for drop-out alleles --- #

    pop_prob_drop <- calc_al_prob_drop(al_csp, pop_al, pop_prob)
    pop_al <- c(pop_al, 99)
    pop_prob <- c(pop_prob, pop_prob_drop)

    index_target <- match(peaks, pop_al)
    pop_al_list[[i]] <- pop_al[index_target]
    pop_prob_list[[i]] <- pop_prob[index_target]
  }

  names(pop_al_list) <- markers_use
  names(pop_prob_list) <- markers_use
  return(list("pop_al_list" = pop_al_list,
              "pop_prob_list" = pop_prob_list))
}

#' set_al_prob_randgt
#'
#' @description The function to set allele probabilities for generating random genotypes
#' @param markers_use A vector of the autosomal markers used
#' @param dt_af The data.table of the allele frequencies
set_al_prob_randgt <- function(markers_use, dt_af){
  n_mk <- length(markers_use)
  pop_al_all <- dt_af[, Allele]
  pop_al_list <- pop_prob_list <- vector(mode = "list", length = n_mk)
  for(i in 1:n_mk){
    mk <- markers_use[i]
    pop_freq <- as.numeric(unlist(dt_af[, mk, with = FALSE]))
    pop_al <- pop_al_all[!is.na(pop_freq)]
    pop_freq <- pop_freq[!is.na(pop_freq)]

    # --- Calculate the expected allele probabilities based on the Dirichlet distribution --- #

    tmp <- correct_af_dirichlet(pop_al, pop_freq, unobs_al = numeric(0), add_unobs = FALSE)
    pop_al_list[[i]] <- tmp$pop_al
    pop_prob_list[[i]] <- tmp$pop_prob
  }

  names(pop_al_list) <- markers_use
  names(pop_prob_list) <- markers_use
  return(list("pop_al_list" = pop_al_list,
              "pop_prob_list" = pop_prob_list))
}
