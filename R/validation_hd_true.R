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

    pop_al_list[[i]] <- pop_al
    pop_prob_list[[i]] <- pop_freq / sum(pop_freq)
  }

  names(pop_al_list) <- markers_use
  names(pop_prob_list) <- markers_use
  return(list("pop_al_list" = pop_al_list,
              "pop_prob_list" = pop_prob_list))
}

#' generate_randgt
#'
#' @description The function to generate random genotypes
#' @param np The number of persons generated
#' @param pop_al_list A list of the alleles in population database
#' @param pop_prob_list A list of the allele probabilities in population database
generate_randgt <- function(np, pop_al_list, pop_prob_list){
  mapply(sample, sapply(pop_al_list, as.numeric), 2 * np, prob = pop_prob_list, replace = TRUE)
}
