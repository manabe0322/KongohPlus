#' get_pae_data
#'
#' @description The function to get the proportion of the locus-specific amplification efficiency
#' @param height_list A list of peak heights
get_pae_data <- function(height_list){
  n_sample <- length(height_list)
  markers <- names(height_list[[1]])
  n_mk <- length(markers)

  pae_data <- matrix(0, n_mk, n_sample)
  for(i in 1:n_sample){
    tph_per_mk <- sapply(height_list[[i]], sum)
    sum_tph <- sum(tph_per_mk)
    if(sum_tph > 0){
      pae_data[, i] <- tph_per_mk / sum_tph
    }
  }

  return(pae_data)
}
