#' generate_randgt
#'
#' @description The function to generate random genotypes
#' @param np The number of persons generated
#' @param pop_al_list A list of the alleles in population database
#' @param pop_prob_list A list of the allele probabilities in population database
generate_randgt <- function(np, pop_al_list, pop_prob_list){
  mapply(sample, sapply(pop_al_list, as.numeric), 2 * np, prob = pop_prob_list, replace = TRUE)
}
