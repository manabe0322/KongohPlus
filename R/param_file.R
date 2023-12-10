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

#' check_mk_est_par
#'
#' @description The function to extract sample names with inappropriate marker sets
#' @param dt The input data.table for estimating parameters
check_mk_est_par <- function(dt){
  setkey(dt, SampleName)
  sample_names <- unique(dt[, SampleName])
  n_sample <- length(sample_names)

  dt_one_sample <- dt[sample_names[1]]
  markers_one_sample <- dt_one_sample[, Marker]

  bool_set_mk <- rep(TRUE, n_sample)
  for(i in 2:n_sample){
    dt_another_sample <- dt[sample_names[i]]
    markers_another_sample <- dt_another_sample[, Marker]
    if(!setequal(markers_one_sample, markers_another_sample)){
      bool_set_mk[i] <- FALSE
    }
  }

  if(all(bool_set_mk)){
    sn_mk_problem <- character(0)
  }else if(length(which(bool_set_mk)) > length(which(!bool_set_mk))){
    sn_mk_problem <- sample_names[!bool_set_mk]
  }else{
    sn_mk_problem <- sample_names[bool_set_mk]
  }

  return(sn_mk_problem)
}

#' search_kit
#'
#' @descrpition The function to search kit names with the target marker set
#' @param markers The vector of the target marker set
#' @param path_pack The package path
search_kit <- function(markers, path_pack){
  # List up all registered kits
  kit_names <- list.files(paste0(path_pack, "/extdata/kit"))
  kit_names <- gsub(".csv", "", kit_names)

  # Search kit names with the target marker set
  kit <- character(0)
  for(i in 1:length(kit_names)){
    kit_info <- fread(paste0(path_pack, "/extdata/kit/", kit_names[i], ".csv"))
    kit_markers <- kit_info[, Marker]
    if(setequal(markers, kit_markers)){
      kit <- c(kit, kit_names[i])
    }
  }
  return(kit)
}

#' extract_sex_markers
#'
#' @description The function to extract sex chromosomal markers
#' @param kit A kit name
#' @param path_pack The package path
extract_sex_markers <- function(kit, path_pack){
  kit_info <- fread(paste0(path_pack, "/extdata/kit", kit, ".csv"))
  setkey(kit_info, Sex_Chr)
  kit_info_sex_chr <- kit_info["Yes"]
  sex_markers <- kit_info_sex_chr[, Marker]
  return(sex_markers)
}

#' check_file_est_par
#'
#' @description The function to check files for estimating parameters
#' @param dt_exp The data.table of experimental data
#' @param dt_gt The data.table of genotype data
#' @param dt_seq The data.table of sequence data
#' @param path_pack The package path
check_file_est_par <- function(dt_exp, dt_gt, dt_seq, path_pack){
  error_message <- ""

  # Whether all required data is loaded or not
  bool_loaded <- all(!is.null(dt_exp), !is.null(dt_gt), !is.null(dt_seq))
  if(!bool_loaded){
    error_message <- "Load required file(s)!"
  }else{

    # Whether the marker set of the experimental data is appropriate or not
    sn_mk_problem_exp <- check_mk_est_par(dt_exp)
    if(length(sn_mk_problem_exp) > 0){
      error_message <- ""
    }else{

      # Whether the marker set is registered or not
      markers_exp <- unique(dt_exp[, Marker])
      kit_exp <- search_kit(markers_exp, path_pack)
      if(length(kit_exp) == 0){
        error_message <- "Kit information for the experimental data has not been registered!"
      }else{

        # Whether the marker set of the genotype data is appropriate or not
        sn_mk_problem_gt <- check_mk_est_par(dt_gt)
        if(length(sn_mk_problem_gt) > 0){
          error_message <- ""
        }else{

          # Whether the marker set of the experimental data is the same as that of the genotype data or not
          markers_gt <- unique(dt_gt[, Marker])
          if(!setequal(markers_exp, markers_gt)){
            error_message <- "The marker set of the experimental data is not the same as that of the genotype data!"
          }else{

            # Whether the marker set for estimating parameters is included in the marker set of the sequence data or not
            sex_markers <- extract_sex_markers(kit_exp[1], path_pack)
            markers_par <- setdiff(markers_exp, sex_markers)
            markers_seq <- unique(dt_seq[, Marker])
            if(!all(is.element(markers_par, markers_seq))){
              error_message <- "Some markers of the experimental data are not included in the sequence data!"
            }else{

              # Whether sample names of the experimental data are the same as those of the genotype data or not
              sn_exp <- unique(dt_exp[, SampleName])
              sn_gt <- unique(dt_gt[, SampleName])
              if(!setequal(sn_exp, sn_gt)){
                error_message <- "Sample names of the experimental data are not the same as those of the genotype data!"
              }
            }
          }
        }
      }
    }
  }
  return(error_message)
}
