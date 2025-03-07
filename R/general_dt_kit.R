#' change_data_type_kit
#'
#' @description The function to change data type of dt_kit
#' @param dt_kit The data.table of the information on the kit
change_data_type_kit <- function(dt_kit){
  col_char <- c("Marker", "Dye", "Chr", "Allele")
  dt_kit[, (col_char) := lapply(.SD, as.character), .SDcols = col_char]
  col_numeric <- c("Repeat_unit_length", "Size")
  dt_kit[, (col_numeric) := lapply(.SD, as.numeric), .SDcols = col_numeric]
  return(dt_kit)
}

#' extract_kit_auto
#'
#' @description The function to extract information on the autosomal markers from the data.table of the information on the kit
#' @param dt_kit The data.table of the information on the kit
extract_kit_auto <- function(dt_kit){
  dt_kit_auto <- copy(dt_kit)
  index_auto <- which(dt_kit_auto[, Chr] == "Auto")
  dt_kit_auto <- dt_kit_auto[index_auto]
  return(dt_kit_auto)
}

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
