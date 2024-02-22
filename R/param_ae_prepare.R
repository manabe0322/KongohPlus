#' get_ae_data
#'
#' @description The function to get data of the locus-specific amplification efficiency
#' @param height_list A list of peak heights
get_ae_data <- function(height_list){
  n_sample <- length(height_list)
  markers <- names(height_list[[1]])
  n_mk <- length(markers)

  height_ae <- matrix(0, n_mk, n_sample)
  for(i in 1:n_sample){
    height_ae[, i] <- sapply(height_list[[i]], sum)
  }

  # Extract used data (i.e, removing AE = 0)
  height_ae <- height_ae[, apply(height_ae, 2, prod) > 0, drop = FALSE]
  n_sample_use <- ncol(height_ae)
  if(n_sample_use > 0){
    ae_data <- matrix(0, n_mk, n_sample_use)
    rownames(ae_data) <- markers
    colnames(ae_data) <- names(height_list)

    # Calculate the locus-specific amplification efficiency
    for(i in 1:n_sample_use){
      ae_data[, i] <- height_ae[, i] / sum(height_ae[, i])
    }
  }else{
    ae_data <- NULL
  }

  return(ae_data)
}
