#' change_data_type_csp
#'
#' @description The function to change data type of dt_csp
#' @param dt_csp The data.table of the crime stain profile
change_data_type_csp <- function(dt_csp){
  names(dt_csp) <- gsub(" ", "", names(dt_csp))

  name_dt_csp <- names(dt_csp)
  col_char <- c("SampleName", "Marker", name_dt_csp[grep("Allele", name_dt_csp)])
  dt_csp[, (col_char) := lapply(.SD, as.character), .SDcols = col_char]
  col_numeric <- c(name_dt_csp[grep("Size", name_dt_csp)], name_dt_csp[grep("Height", name_dt_csp)])
  dt_csp[, (col_numeric) := lapply(.SD, as.numeric), .SDcols = col_numeric]

  dt_csp[which(dt_csp == "", arr.ind = TRUE)] <- NA

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

  name_dt_af <- names(dt_af)
  col_char <- c("Allele")
  dt_af[, (col_char) := lapply(.SD, as.character), .SDcols = col_char]
  col_numeric <- setdiff(name_dt_af, "Allele")
  dt_af[, (col_numeric) := lapply(.SD, as.numeric), .SDcols = col_numeric]

  dt_af[which(dt_af == "", arr.ind = TRUE)] <- NA

  return(dt_af)
}
