#' check_dt_mk_est_par
#'
#' @description The function to extract sample names with inappropriate marker sets from a data.table
#' @param dt The input data.table for estimating parameters
check_dt_mk_est_par <- function(dt, markers_par){
  # Extract sample names
  setkey(dt, SampleName)
  sample_names <- unique(dt[, SampleName])
  n_sample <- length(sample_names)

  # Check sample names with inappropriate marker sets
  bool_sn_problem <- rep(FALSE, n_sample)
  for(i in 1:n_sample){
    dt_one_sample <- dt[sample_names[i]]
    markers_one_sample <- dt_one_sample[, Marker]
    bool_sn_problem[i] <- !all(is.element(markers_par, markers_one_sample))
  }
  return(sample_names[bool_sn_problem])
}

#' check_est_par
#'
#' @description The function to check files for estimating parameters
#' @param dt_exp The data.table of the experimental data
#' @param dt_gt The data.table of the genotype data
#' @param dt_seq The data.table of the sequence data
#' @param dt_kit_set The data.table of the kit-setting
check_est_par <- function(dt_exp, dt_gt, dt_seq, dt_kit_set){
  error_message <- ""

  # Whether all required data is loaded or not
  if(is.null(dt_exp)){
    error_message <- paste0(error_message, "Load experimental data!<br>")
  }
  if(is.null(dt_gt)){
    error_message <- paste0(error_message, "Load genotype data!<br>")
  }
  if(is.null(dt_seq)){
    error_message <- paste0(error_message, "Set sequence data!<br>")
  }

  # Check the analytical thresholds
  markers_par <- dt_kit_set[, Marker]
  at_all_mk <- dt_kit_set[, AT]
  bool_without_at <- is.na(at_all_mk)
  if(any(bool_without_at)){
    error_message <- paste0(error_message, "Set the analytical threshold(s) of the following marker(s): ", paste0(markers_par[bool_without_at], collapse = ", "), ".<br>")
  }

  if(all(!is.null(dt_exp), !is.null(dt_gt), !is.null(dt_seq))){
    # Check the marker set of the experimental data
    sn_mk_problem_exp <- check_dt_mk_est_par(dt_exp, markers_par)
    if(length(sn_mk_problem_exp) > 0){
      error_message <- paste0(error_message, "Inappropriate marker set of the following sample(s) in the experimental data: ", paste0(sn_mk_problem_exp, collapse = ", "), ".<br>")
    }

    # Check the marker set of the genotype data
    sn_mk_problem_gt <- check_dt_mk_est_par(dt_gt, markers_par)
    if(length(sn_mk_problem_gt) > 0){
      error_message <- paste0(error_message, "Inappropriate marker set of the following sample(s) in the genotype data: ", paste0(sn_mk_problem_gt, collapse = ", "), ".<br>")
    }

    # Check the marker set of the sequence data
    markers_seq <- unique(dt_seq[, Marker])
    markers_without_seq <- setdiff(markers_par, markers_seq)
    if(length(markers_without_seq) > 0){
      error_message <- paste0(error_message, "The following marker(s) are not included in the sequence data: ", paste0(markers_without_seq, collapse = ", "), ".<br>")
    }

    # Check the sample names
    sn_exp <- unique(dt_exp[, SampleName])
    sn_gt <- unique(dt_gt[, SampleName])
    sn_without_gt <- setdiff(sn_exp, sn_gt)
    if(length(sn_without_gt) > 0){
      error_message <- paste0(error_message, "No genotype data of the following samples: ", paste0(sn_without_gt, collapse = ", "), ".<br>")
    }
  }

  return(error_message)
}
