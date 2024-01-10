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

#' calc_par_al_bs
#'
#' @description The function to calculate the parent allele of the back stutter
#' @param bs The allele call of the back stutter
calc_par_al_bs <- function(bs){
  return(round(bs + 1, 1))
}

#' calc_par_al_fs
#'
#' @description The function to calculate the parent allele of the forward stutter
#' @param fs The allele call of the forward stutter
calc_par_al_fs <- function(fs){
  return(round(fs - 1, 1))
}

#' calc_par_al_ds
#'
#' @description The function to calculate the parent allele of the double-back stutter
#' @param ds The allele call of the double-back stutter
calc_par_al_ds <- function(ds){
  return(round(ds + 2, 1))
}

#' calc_par_al_m2s
#'
#' @description The function to calculate the parent allele of the minus 2-nt stutter
#' @param m2s The allele call of the minus 2-nt stutter
calc_par_al_m2s <- function(m2s){
  n_m2s <- length(m2s)
  al <- rep(0, n_m2s)
  for(i in 1:n_m2s){
    if(m2s[i] %% 1 < 0.15){
      al[i] <- round(m2s[i] + 0.2, 1)
    }else{
      al[i] <- round(m2s[i] + 0.8, 1)
    }
  }
  return(al)
}
