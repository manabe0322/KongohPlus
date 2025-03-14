#' extract_dye
#'
#' @description The function to extract dyes from dt_kit
#' @param dt_kit The data.table of the kit information
#' @param markers A vector of the markers
extract_dye <- function(dt_kit, markers){
  setindex(dt_kit, Marker)
  n_mk <- length(markers)
  dyes <- rep("", n_mk)
  for(i in 1:n_mk){
    dt_kit_mk1 <- dt_kit[Marker == markers[i]]
    dyes[i] <- unique(dt_kit_mk1[, Dye])
  }
  names(dyes) <- markers
  return(dyes)
}
