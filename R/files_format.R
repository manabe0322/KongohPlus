#' change_data_type_csp
#'
#' @description The function to change data type of dt_csp
#' @param dt_csp The data.table of the crime stain profile
change_data_type_csp <- function(dt_csp){
  names(dt_csp) <- gsub(" ", "", names(dt_csp))

  name_dt_csp <- names(dt_csp)
  index_use <- c(which(name_dt_csp == "SampleName"), which(name_dt_csp == "Marker"), grep("Allele", name_dt_csp), grep("Height", name_dt_csp))
  dt_csp <- dt_csp[, index_use, with = FALSE]

  name_dt_csp <- names(dt_csp)
  col_char <- c("SampleName", "Marker", name_dt_csp[grep("Allele", name_dt_csp)])
  dt_csp[, (col_char) := lapply(.SD, as.character), .SDcols = col_char]
  col_numeric <- c(name_dt_csp[grep("Size", name_dt_csp)], name_dt_csp[grep("Height", name_dt_csp)])
  dt_csp[, (col_numeric) := lapply(.SD, as.numeric), .SDcols = col_numeric]

  dt_csp[which(dt_csp == "", arr.ind = TRUE)] <- NA

  # Remove columns with all NA data
  tmp <- lapply(dt_csp, is.na)
  index_col_with_not_na <- which(!sapply(tmp, all))
  index_keep <- sort(unique(c(which(name_dt_csp == "SampleName"), which(name_dt_csp == "Marker"), index_col_with_not_na)))
  dt_csp <- dt_csp[, index_keep, with = FALSE]

  return(dt_csp)
}

#' change_data_type_ref
#'
#' @description The function to change data type of dt_ref
#' @param dt_ref The data.table of the reference profile
change_data_type_ref <- function(dt_ref){
  names(dt_ref) <- gsub(" ", "", names(dt_ref))

  name_dt_ref <- names(dt_ref)
  col_char <- c("SampleName", "Marker", name_dt_ref[grep("Allele", name_dt_ref)])
  dt_ref[, (col_char) := lapply(.SD, as.character), .SDcols = col_char]

  dt_ref[which(dt_ref == "", arr.ind = TRUE)] <- NA

  return(dt_ref)
}

#' change_data_type_af
#'
#' @description The function to change data type of dt_af
#' @param dt_csp The data.table of the allele frequency
change_data_type_af <- function(dt_af){
  names(dt_af) <- gsub(" ", "", names(dt_af))

  #name_dt_af <- names(dt_af)
  #col_char <- c("Allele")
  #dt_af[, (col_char) := lapply(.SD, as.character), .SDcols = col_char]
  col_numeric <- names(dt_af)
  dt_af[, (col_numeric) := lapply(.SD, as.numeric), .SDcols = col_numeric]

  dt_af[which(dt_af == "", arr.ind = TRUE)] <- NA

  return(dt_af)
}

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
