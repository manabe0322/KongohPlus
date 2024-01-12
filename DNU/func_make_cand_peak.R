#' make_cand_peak
#'
#' @description The function to make detectable peaks in the crime stain profile under the designated alleles
#' @param alleles A numeric vector of alleles
#' @param sr_consider Information on consideration of back, forward, double-back, and minus 2-nt stutters
make_cand_peak <- function(alleles, sr_consider){

  # Define the output object
  cand_peak <- alleles

  # Whether back stutters are considered or not
  if(sr_consider[1]){
    cand_peak <- c(cand_peak, calc_bs(alleles))
  }

  # Whether forward stutters are considered or not
  if(sr_consider[2]){
    cand_peak <- c(cand_peak, calc_fs(alleles))
  }

  # whether double-back stutters are considered or not
  if(sr_consider[3]){
    cand_peak <- c(cand_peak, calc_ds(alleles))
  }

  # Whether minus 2-nt stutters are considered or not
  if(sr_consider[4]){
    cand_peak <- c(cand_peak, calc_m2s(alleles))
  }

  # Return
  return(sort(unique(cand_peak)))
}
