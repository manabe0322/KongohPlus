#' change_data_type_exp
#'
#' @description The function to change data type of dt_exp
#' @param dt_exp The data.table of experimental data
change_data_type_exp <- function(dt_exp){
  names(dt_exp) <- gsub(" ", "", names(dt_exp))

  name_dt_exp <- names(dt_exp)
  index_use <- c(which(name_dt_exp == "SampleName"), which(name_dt_exp == "Marker"), grep("Allele", name_dt_exp), grep("Height", name_dt_exp))
  dt_exp <- dt_exp[, index_use, with = FALSE]

  col_char <- c("SampleName", "Marker", name_dt_exp[grep("Allele", name_dt_exp)])
  dt_exp[, (col_char) := lapply(.SD, as.character), .SDcols = col_char]
  col_numeric <- c(name_dt_exp[grep("Height", name_dt_exp)])
  dt_exp[, (col_numeric) := lapply(.SD, as.numeric), .SDcols = col_numeric]

  dt_exp[which(dt_exp == "", arr.ind = TRUE)] <- NA

  return(dt_exp)
}

#' change_data_type_gt
#'
#' @description The function to change data type of dt_gt
#' @param dt_gt The data.table of genotype data
change_data_type_gt <- function(dt_gt){
  names(dt_gt) <- gsub(" ", "", names(dt_gt))

  name_dt_gt <- names(dt_gt)
  index_use <- c(which(name_dt_gt == "SampleName"), which(name_dt_gt == "Marker"), grep("Allele", name_dt_gt))
  dt_gt <- dt_gt[, index_use, with = FALSE]

  col_char <- c("SampleName", "Marker", name_dt_gt[grep("Allele", name_dt_gt)])
  dt_gt[, (col_char) := lapply(.SD, as.character), .SDcols = col_char]

  dt_gt[which(dt_gt == "", arr.ind = TRUE)] <- NA

  return(dt_gt)
}
