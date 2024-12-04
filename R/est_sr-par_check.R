#' check_dt_mk
#'
#' @description The function to extract sample names with inappropriate marker sets from a data.table
#' @param dt The input data.table for estimating average stutter ratios
#' @param markers_est_sr The target markers for estimating average stutter ratios
check_dt_mk <- function(dt, markers_est_sr){
  # Extract sample names
  setindex(dt, SampleName)
  sample_names <- unique(dt$SampleName)
  n_sample <- length(sample_names)

  # Check sample names with inappropriate marker sets
  bool_sn_problem <- rep(FALSE, n_sample)
  for(i in 1:n_sample){
    dt_one_sample <- dt[SampleName == sample_names[i]]
    markers_one_sample <- dt_one_sample$Marker
    bool_sn_problem[i] <- !all(is.element(markers_est_sr, markers_one_sample))
  }
  return(sample_names[bool_sn_problem])
}

#' check_est_sr
#'
#' @description The function to check files for estimating average stutter ratios
#' @param dt_exp The data.table of the experimental data
#' @param dt_gt The data.table of the genotype data
#' @param dt_seq The data.table of the sequence data
#' @param dt_setting_est_sr The data.table in the setting tab for estimating average stutter ratios
check_est_sr <- function(dt_exp, dt_gt, dt_seq, dt_setting_est_sr){
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
  markers_est_sr <- dt_setting_est_sr$Marker
  at_all_mk <- dt_setting_est_sr$AT
  bool_without_at <- is.na(at_all_mk)
  if(any(bool_without_at)){
    error_message <- paste0(error_message, "Set the analytical threshold(s) of the following marker(s): ", paste0(markers_est_sr[bool_without_at], collapse = ", "), ".<br>")
  }

  # Check the BSR model
  model_bsr <- dt_setting_est_sr[["BSR model"]]
  index_without_model_bsr <- which(model_bsr == "")
  if(length(index_without_model_bsr) != 0){
    error_message <- paste0(error_message, "Set the BSR model of the following marker(s): ", paste0(markers_est_sr[index_without_model_bsr], collapse = ", "), ".<br>")
  }

  # Check the FSR model
  model_fsr <- dt_setting_est_sr[["FSR model"]]
  index_without_model_fsr <- which(model_fsr == "")
  if(length(index_without_model_fsr) != 0){
    error_message <- paste0(error_message, "Set the FSR model of the following marker(s): ", paste0(markers_est_sr[index_without_model_fsr], collapse = ", "), ".<br>")
  }

  # Check the DSR model
  model_dsr <- dt_setting_est_sr[["DSR model"]]
  index_without_model_dsr <- which(model_dsr == "")
  if(length(index_without_model_dsr) != 0){
    error_message <- paste0(error_message, "Set the DSR model of the following marker(s): ", paste0(markers_est_sr[index_without_model_dsr], collapse = ", "), ".<br>")
  }

  # Check the M2SR model
  model_m2sr <- dt_setting_est_sr[["M2SR model"]]
  index_without_model_m2sr <- which(model_m2sr == "")
  if(length(index_without_model_m2sr) != 0){
    error_message <- paste0(error_message, "Set the M2SR model of the following marker(s): ", paste0(markers_est_sr[index_without_model_m2sr], collapse = ", "), ".<br>")
  }

  # Check the minimum number of repeats that affects stuttering (for multi-seq model)
  param_x <- dt_setting_est_sr[["Min. number of repeats that affects stuttering (for Multi-seq)"]]
  index_multiseq <- sort(unique(c(which(model_bsr == "Multi-seq"), which(model_fsr == "Multi-seq"), which(model_dsr == "Multi-seq"), which(model_m2sr == "Multi-seq"))))
  if(length(index_multiseq) != 0){
    index_without_x <- which(is.na(param_x))
    index_error_x <- intersect(index_multiseq, index_without_x)
    if(length(index_error_x) != 0){
      error_message <- paste0(error_message, "Set the minimum number of repeats that affects stuttering of the following marker(s): ", paste0(markers_est_sr[index_error_x], collapse = ", "), ".<br>")
    }
  }

  if(all(!is.null(dt_exp), !is.null(dt_gt), !is.null(dt_seq))){
    # Check the marker set of the experimental data
    sn_mk_problem_exp <- check_dt_mk(dt_exp, markers_est_sr)
    if(length(sn_mk_problem_exp) > 0){
      error_message <- paste0(error_message, "Inappropriate marker set of the following sample(s) in the experimental data: ", paste0(sn_mk_problem_exp, collapse = ", "), ".<br>")
    }

    # Check the marker set of the genotype data
    sn_mk_problem_gt <- check_dt_mk(dt_gt, markers_est_sr)
    if(length(sn_mk_problem_gt) > 0){
      error_message <- paste0(error_message, "Inappropriate marker set of the following sample(s) in the genotype data: ", paste0(sn_mk_problem_gt, collapse = ", "), ".<br>")
    }

    # Check the marker set of the sequence data
    markers_seq <- unique(dt_seq$Marker)
    markers_without_seq <- setdiff(markers_est_sr, markers_seq)
    if(length(markers_without_seq) > 0){
      error_message <- paste0(error_message, "The following marker(s) are not included in the sequence data: ", paste0(markers_without_seq, collapse = ", "), ".<br>")
    }

    # Check the sample names
    sn_exp <- unique(dt_exp$SampleName)
    sn_gt <- unique(dt_gt$SampleName)
    sn_without_gt <- setdiff(sn_exp, sn_gt)
    if(length(sn_without_gt) > 0){
      error_message <- paste0(error_message, "No genotype data of the following samples: ", paste0(sn_without_gt, collapse = ", "), ".<br>")
    }
  }

  return(error_message)
}

#' check_est_par
#'
#' @description The function to check files for estimating parameters
#' @param dt_exp The data.table of the experimental data
#' @param dt_gt The data.table of the genotype data
#' @param dt_setting_est_sr The data.table in the setting tab for estimating parameters
check_est_par <- function(dt_exp, dt_gt, dt_setting_est_par){
  error_message <- ""

  # Whether all required data is loaded or not
  if(is.null(dt_exp)){
    error_message <- paste0(error_message, "Load experimental data!<br>")
  }
  if(is.null(dt_gt)){
    error_message <- paste0(error_message, "Load genotype data!<br>")
  }

  # Check the analytical thresholds
  markers_est_par <- dt_setting_est_par$Marker
  at_all_mk <- dt_setting_est_par$AT
  bool_without_at <- is.na(at_all_mk)
  if(any(bool_without_at)){
    error_message <- paste0(error_message, "Set the analytical threshold(s) of the following marker(s): ", paste0(markers_est_par[bool_without_at], collapse = ", "), ".<br>")
  }

  if(all(!is.null(dt_exp), !is.null(dt_gt))){
    # Check the marker set of the experimental data
    sn_mk_problem_exp <- check_dt_mk(dt_exp, markers_est_par)
    if(length(sn_mk_problem_exp) > 0){
      error_message <- paste0(error_message, "Inappropriate marker set of the following sample(s) in the experimental data: ", paste0(sn_mk_problem_exp, collapse = ", "), ".<br>")
    }

    # Check the marker set of the genotype data
    sn_mk_problem_gt <- check_dt_mk(dt_gt, markers_est_par)
    if(length(sn_mk_problem_gt) > 0){
      error_message <- paste0(error_message, "Inappropriate marker set of the following sample(s) in the genotype data: ", paste0(sn_mk_problem_gt, collapse = ", "), ".<br>")
    }

    # Check the sample names
    sn_exp <- unique(dt_exp$SampleName)
    sn_gt <- unique(dt_gt$SampleName)
    sn_without_gt <- setdiff(sn_exp, sn_gt)
    if(length(sn_without_gt) > 0){
      error_message <- paste0(error_message, "No genotype data of the following samples: ", paste0(sn_without_gt, collapse = ", "), ".<br>")
    }
  }

  return(error_message)
}
