#' get_ae_data
#'
#' @description The function to get data of the locus-specific amplification efficiency
#' @param height_list A list of peak heights
get_ae_data <- function(height_list){
  n_sample <- length(height_list)
  markers <- names(height_list[[1]])
  n_mk <- length(markers)

  height_ae <- matrix(0, n_sample, n_mk)
  rownames(height_ae) <- names(height_list)
  colnames(height_ae) <- markers

  for(i in 1:n_sample){
    height_ae[i, ] <- sapply(height_list[[i]], sum)
  }

  # Extract used data (i.e, removing AE = 0)
  height_ae <- height_ae[apply(height_ae, 1, prod) > 0, , drop = FALSE]
  n_sample_use <- nrow(height_ae)
  if(n_sample_use > 0){
    ae_data <- matrix(0, n_sample_use, n_mk)
    rownames(ae_data) <- rownames(height_ae)
    colnames(ae_data) <- markers

    # Calculate the locus-specific amplification efficiency
    for(i in 1:n_sample_use){
      ae_data[i, ] <- height_ae[i, ] / mean(height_ae[i, ])
    }
  }else{
    ae_data <- height_ae
  }

  ae_list <- list()
  ae_list[[1]] <- ae_data
  ae_list[[2]] <- height_ae
  names(ae_list) <- c("AE", "Height_AE")
  return(ae_list)
}
