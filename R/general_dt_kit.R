#' change_data_type_kit
#'
#' @description The function to change data type of dt_kit
#' @param dt_kit The data.table of the information on the kit
change_data_type_kit <- function(dt_kit){
  if(all(is.element(c("Marker", "Chromosome", "Allele", "Repeat Unit Length", "Corrected Allele", "Size"), names(dt_kit)))){
    col_char <- c("Marker", "Chromosome", "Allele")
    dt_kit[, (col_char) := lapply(.SD, as.character), .SDcols = col_char]
    col_integer <- c("Repeat Unit Length")
    dt_kit[, (col_integer) := lapply(.SD, as.integer), .SDcols = col_integer]
    col_numeric <- c("Corrected Allele", "Size")
    dt_kit[, (col_numeric) := lapply(.SD, as.numeric), .SDcols = col_numeric]
    return(dt_kit)
  }else{
    stop("Error: The selected file format is invalid.")
  }
}

#' extract_kit_auto
#'
#' @description The function to extract information on the autosomal markers from the data.table of the information on the kit
#' @param dt_kit The data.table of the information on the kit
extract_kit_auto <- function(dt_kit){
  dt_kit_auto <- copy(dt_kit)
  index_auto <- which(dt_kit_auto[["Chromosome"]] == "Autosomal")
  dt_kit_auto <- dt_kit_auto[index_auto]
  return(dt_kit_auto)
}
