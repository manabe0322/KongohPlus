#' correct_af_dirichlet
#'
#' @description The function to correct allele frequencies based on the Dirichlet distribution
#' @param pop_al A vector of the alleles in the population database
#' @param pop_freq A vector of the allele frequencies in the population database
#' @param unobs_al A vector of the unobserved alleles in the population database
correct_af_dirichlet <- function(pop_al, pop_freq, unobs_al = numeric(0)){
  n_unobs <- length(unobs_al)
  if(n_unobs > 0){
    pop_al <- c(pop_al, unobs_al)
    pop_freq <- c(pop_freq, rep(0, n_unobs))
  }
  pop_prob <- (pop_freq + 1) / (sum(pop_freq) + length(pop_al))

  order_al <- order(pop_al)
  pop_al <- pop_al[order_al]
  pop_prob <- pop_prob[order_al]

  return(list("pop_al" = pop_al,
              "pop_prob" = pop_prob))
}

#' add_pop_prob_drop
#'
#' @description The function to add an allele probability for drop-out
#' @param al_csp A vector of alleles detected in crime stain profiles
#' @param pop_al A vector of the alleles in the population database
#' @param pop_prob A vector of the allele probabilities
add_pop_prob_drop <- function(al_csp, pop_al, pop_prob){
  index_drop <- which(!is.element(pop_al, al_csp))
  if(length(index_drop) != 0){
    pop_prob_drop <- sum(pop_prob[index_drop])
  }else{
    pop_prob_drop <- 0
  }
  pop_al <- c(pop_al, 99)
  pop_prob <- c(pop_prob, pop_prob_drop)
  return(list("pop_al" = pop_al,
              "pop_prob" = pop_prob))
}

#' set_pop_prob
#'
#' @description The function to set allele probabilities
#' @param markers_auto A vector of the autosomal markers
#' @param dt_af The data.table of the allele frequencies
#' @param peaks_list The peak list of all samples
set_pop_prob <- function(markers_auto, dt_af, peaks_list){
  n_mk <- length(markers_auto)
  n_sn_csp <- length(peaks_list)
  pop_al_all <- dt_af[, Allele]
  pop_al_list <- pop_prob_list <- vector(mode = "list", length = n_mk)
  for(i in 1:n_mk){
    mk <- markers_auto[i]
    pop_freq <- as.numeric(unlist(dt_af[, mk, with = FALSE]))
    pop_al <- pop_al_all[!is.na(pop_freq)]
    pop_freq <- pop_freq[!is.na(pop_freq)]

    al_csp <- numeric(0)
    for(j in 1:n_sn_csp){
      peaks_rm99 <- setdiff(peaks_list[[j]][[mk]], 99)
      al_csp <- c(al_csp, peaks_rm99)
    }

    unobs_al <- setdiff(al_csp, pop_al)

    tmp <- correct_af_dirichlet(pop_al, pop_freq, unobs_al)
    tmp <- add_pop_prob_drop(al_csp, tmp$pop_al, tmp$pop_prob)
    pop_al_list[[i]] <- tmp$pop_al
    pop_prob_list[[i]] <- tmp$pop_prob
  }

  names(pop_al_list) <- markers_auto
  names(pop_prob_list) <- markers_auto
  return(list("pop_al_list" = pop_al_list,
              "pop_prob_list" = pop_prob_list))
}
