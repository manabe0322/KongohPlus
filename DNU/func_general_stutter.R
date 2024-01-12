#' calc_bs
#'
#' @description The function to calculate the allele call of the back stutter
#' @param al The parent allele of the back stutter
calc_bs <- function(al){
  return(round(al - 1, 1))
}

#' calc_fs
#'
#' @description The function to calculate the allele call of the forward stutter
#' @param al The parent allele of the forward stutter
calc_fs <- function(al){
  return(round(al + 1, 1))
}

#' calc_ds
#'
#' @description The function to calculate the allele call of the double-back stutter
#' @param al The parent allele of the double-back stutter
calc_ds <- function(al){
  return(round(al - 2, 1))
}

#' calc_m2s
#'
#' @description The function to calculate the allele call of the minus 2-nt stutter
#' @param al The parent allele of the minus 2-nt stutter
calc_m2s <- function(al){
  n_al <- length(al)
  m2s <- rep(0, n_al)
  for(i in 1:n_al){
    if(al[i] %% 1 < 0.15){
      m2s[i] <- round(al[i] - 0.8, 1)
    }else{
      m2s[i] <- round(al[i] - 0.2, 1)
    }
  }
  return(m2s)
}
