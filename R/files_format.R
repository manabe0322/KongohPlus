#' change_data_type_csp
#'
#' @description The function to change data type of dt_csp
#' @param dt_csp The data.table of the crime stain profile
change_data_type_csp <- function(dt_csp){
  names(dt_csp) <- gsub(" ", "", names(dt_csp))
  name_dt_csp <- names(dt_csp)

  # --- Extract the used columns --- #

  index_use <- c(which(name_dt_csp == "SampleName"), which(name_dt_csp == "Marker"), grep("Allele", name_dt_csp), grep("Height", name_dt_csp))
  dt_csp <- dt_csp[, index_use, with = FALSE]
  name_dt_csp <- names(dt_csp)
  n_allele <- length(grep("Allele", name_dt_csp))
  n_height <- length(grep("Height", name_dt_csp))

  if(all(is.element("SampleName", name_dt_csp), n_allele != 0, n_height != 0, n_allele == n_height)){

    # --- Change data type of each column --- #

    col_char <- c("SampleName", "Marker", name_dt_csp[grep("Allele", name_dt_csp)])
    dt_csp[, (col_char) := lapply(.SD, as.character), .SDcols = col_char]
    col_numeric <- c(name_dt_csp[grep("Size", name_dt_csp)], name_dt_csp[grep("Height", name_dt_csp)])
    dt_csp[, (col_numeric) := lapply(.SD, as.numeric), .SDcols = col_numeric]
    dt_csp[which(dt_csp == "", arr.ind = TRUE)] <- NA

    # --- Remove columns with all NA data --- #

    tmp <- lapply(dt_csp, is.na)
    index_col_with_not_na <- which(!sapply(tmp, all))
    index_keep <- sort(unique(c(which(name_dt_csp == "SampleName"), which(name_dt_csp == "Marker"), index_col_with_not_na)))
    dt_csp <- dt_csp[, index_keep, with = FALSE]

    return(dt_csp)
  }

  else{
    stop("Error: The selected file format is invalid.")
  }
}

#' change_data_type_ref
#'
#' @description The function to change data type of dt_ref
#' @param dt_ref The data.table of the reference profile
change_data_type_ref <- function(dt_ref){
  names(dt_ref) <- gsub(" ", "", names(dt_ref))
  name_dt_ref <- names(dt_ref)

  if(all(is.element(c("SampleName", "Marker", "Allele1", "Allele2"), name_dt_ref)) && length(name_dt_ref) == 4){

    # --- Change data type of each column --- #

    col_char <- c("SampleName", "Marker", name_dt_ref[grep("Allele", name_dt_ref)])
    dt_ref[, (col_char) := lapply(.SD, as.character), .SDcols = col_char]
    dt_ref[which(dt_ref == "", arr.ind = TRUE)] <- NA

    return(dt_ref)
  }

  else{
    stop("Error: The selected file format is invalid.")
  }
}

#' change_data_type_af
#'
#' @description The function to change data type of dt_af
#' @param dt_af The data.table of the allele frequency
change_data_type_af <- function(dt_af){
  names(dt_af) <- gsub(" ", "", names(dt_af))
  name_dt_af <- names(dt_af)

  if(is.element("Allele", name_dt_af)){

    # --- Change data type of each column --- #

    dt_af[, (name_dt_af) := lapply(.SD, as.numeric), .SDcols = name_dt_af]
    dt_af[which(dt_af == "", arr.ind = TRUE)] <- NA

    return(dt_af)
  }

  else{
    stop("Error: The selected file format is invalid.")
  }
}
