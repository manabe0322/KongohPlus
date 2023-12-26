#' make_dt_alcor
#'
#' @description The function to make a data.table for corrected alleles
#' @param markers_par The vector of markers for parameters
#' @param sr_best_model_all_mk The vector of the best BSR model for all markers
#' @param sr_alcor_all_mk The list of the information on corrected alleles of all markers
make_dt_alcor <- function(markers_par, sr_best_model_all_mk, sr_alcor_all_mk){
  vec_mk <- vec_model <- character(0)
  vec_info_al <- vec_info_ca <- numeric(0)

  n_mk <- length(markers_par)
  for(i in 1:n_mk){
    mk <- markers_par[i]
    model <- sr_best_model_all_mk[i]
    info_al <- sr_alcor_all_mk$info_al
    info_ca <- sr_alcor_all_mk$info_ca

    n_al <- length(info_al)
    vec_mk <- c(vec_mk, rep(mk, n_al))
    vec_model <- c(vec_model, rep(model, n_al))
    vec_info_al <- c(vec_info_al, info_al)
    vec_info_ca <- c(vec_info_ca, info_ca)
  }

  dt_alcor <- data.table(Marker = vec_mk,
                         Model = vec_model,
                         Allele = vec_info_al,
                         Corrected = vec_info_ca)
  return(dt_alcor)
}
