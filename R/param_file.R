#' change_data_type_exp
#'
#' @description The function to change data type of dt_exp
#' @param dt_exp The data.table of experimental data
change_data_type_exp <- function(dt_exp){
  names(dt_exp) <- gsub(" ", "", names(dt_exp))

  name_dt_exp <- names(dt_exp)
  col_char <- c("SampleName", "Marker", name_dt_exp[grep("Allele", name_dt_exp)])
  dt_exp[, (col_char) := lapply(.SD, as.character), .SDcols = col_char]
  col_numeric <- c(name_dt_exp[grep("Size", name_dt_exp)], name_dt_exp[grep("Height", name_dt_exp)])
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
  col_char <- c("SampleName", "Marker", name_dt_gt[grep("Allele", name_dt_gt)])
  dt_gt[, (col_char) := lapply(.SD, as.character), .SDcols = col_char]

  dt_gt[which(dt_gt == "", arr.ind = TRUE)] <- NA

  return(dt_gt)
}

#' change_data_type_seq
#'
#' @description The function to change data type of dt_seq
#' @param dt_seq The data.table of sequence data
change_data_type_seq <- function(dt_seq){

  names(dt_seq) <- gsub(" ", "", names(dt_seq))

  name_dt_seq <- names(dt_seq)
  col_char <- c("Marker", "Allele", "RepeatRegion", "FlankingRegion")
  dt_seq[, (col_char) := lapply(.SD, as.character), .SDcols = col_char]
  col_numeric <- c("Count")
  dt_seq[, (col_numeric) := lapply(.SD, as.numeric), .SDcols = col_numeric]

  dt_seq[which(dt_seq == "", arr.ind = TRUE)] <- NA

  return(dt_seq)
}
