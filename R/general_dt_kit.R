#' extract_kit_auto
#'
#' @description The function to extract information on the autosomal markers from the data.table of the information on the kit
#' @param dt_kit The data.table of the information on the kit
extract_kit_auto <- function(dt_kit){
  index_auto <- which(dt_kit[, Chr] == "Auto")
  dt_kit_auto <- dt_kit[index_auto]
  return(dt_kit_auto)
}

#' extract_dye
#'
#' @description The function to extract dyes from dt_kit_auto
#' @param dt_kit_auto The data.table of the kit information on autosomal STRs
#' @param markers_par A vector of the markers for estimating parameters
extract_dye <- function(dt_kit_auto, markers_par){
  setkey(dt_kit_auto, Marker)
  n_mk <- length(markers_par)
  dyes_par <- rep("", n_mk)
  for(i in 1:n_mk){
    dt_kit_auto_mk1 <- dt_kit_auto[markers_par[i]]
    dyes_par[i] <- unique(dt_kit_auto_mk1[, Dye])
  }
  return(dyes_par)
}
