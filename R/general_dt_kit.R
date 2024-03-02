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
#' @description The function to extract dyes from dt_kit
#' @param dt_kit The data.table of the kit information
#' @param markers A vector of the markers
extract_dye <- function(dt_kit, markers){
  setkey(dt_kit, Marker)
  n_mk <- length(markers)
  dyes <- rep("", n_mk)
  for(i in 1:n_mk){
    dt_kit_mk1 <- dt_kit[markers[i]]
    dyes[i] <- unique(dt_kit_mk1[, Dye])
  }
  names(dyes) <- markers
  return(dyes)
}
