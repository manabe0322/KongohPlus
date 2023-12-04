#' change_data_type_exp
#'
#' @description The function to change data type of dt_exp_init
#' @param dt_exp_init The initial data.table of experimental data
change_data_type_exp <- function(dt_exp_init){
  names(dt_exp_init) <- gsub(" ", "", names(dt_exp_init))

  name_dt_exp <- names(dt_exp_init)
  col_char <- c("SampleName", "Marker", name_dt_exp[grep("Allele", name_dt_exp)])
  dt_exp_init[, (col_char) := lapply(.SD, as.character), .SDcols = col_char]
  col_numeric <- c(name_dt_exp[grep("Size", name_dt_exp)], name_dt_exp[grep("Height", name_dt_exp)])
  dt_exp_init[, (col_numeric) := lapply(.SD, as.numeric), .SDcols = col_numeric]

  dt_exp_init[which(dt_exp_init == "", arr.ind = TRUE)] <- NA

  return(dt_exp_init)
}

#' change_data_type_gt
#'
#' @description The function to change data type of dt_gt_init
#' @param dt_gt_init The initial data.table of genotype data
change_data_type_gt <- function(dt_gt_init){
  names(dt_gt_init) <- gsub(" ", "", names(dt_gt_init))

  name_dt_gt <- names(dt_gt_init)
  col_char <- c("SampleName", "Marker", name_dt_gt[grep("Allele", name_dt_gt)])
  dt_gt_init[, (col_char) := lapply(.SD, as.character), .SDcols = col_char]

  dt_gt_init[which(dt_gt_init == "", arr.ind = TRUE)] <- NA

  return(dt_gt_init)
}

#' change_data_type_seq
#'
#' @description The function to change data type of dt_seq_init
#' @param dt_seq_init The initial data.table of sequence data
change_data_type_seq <- function(dt_seq_init){

  names(dt_seq_init) <- gsub(" ", "", names(dt_seq_init))

  name_dt_seq <- names(dt_seq_init)
  col_char <- c("Marker", "Allele", "RepeatRegion", "FlankingRegion")
  dt_seq_init[, (col_char) := lapply(.SD, as.character), .SDcols = col_char]
  col_numeric <- c("Count")
  dt_seq_init[, (col_numeric) := lapply(.SD, as.numeric), .SDcols = col_numeric]

  dt_seq_init[which(dt_seq_init == "", arr.ind = TRUE)] <- NA

  return(dt_seq_init)
}
