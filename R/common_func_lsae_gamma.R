#' create_init_dt_setting_est_par
#'
#' @description The function to create the initial data.table in the setting tab for estimating parameters
#' @param dt_kit The data.table of information on a kit
create_init_dt_setting_est_par <- function(dt_kit){
  dt_kit_auto <- extract_kit_auto(dt_kit)
  markers <- unique(dt_kit_auto$Marker)

  dt_setting_est_par <- data.table(Marker = markers,
                                   AT = rep(35, length(markers)))
  col_int <- c("AT")
  dt_setting_est_par[, (col_int) := lapply(.SD, as.integer), .SDcols = col_int]

  return(dt_setting_est_par)
}

#' check_est_par
#'
#' @description The function to check files for estimating parameters
#' @param dt_exp The data.table of the experimental data
#' @param dt_gt The data.table of the genotype data
#' @param dt_setting_est_par The data.table in the setting tab for estimating parameters
check_est_par <- function(dt_exp, dt_gt, dt_setting_est_par){
  error_message <- ""

  # --- Check whether all required data is loaded --- #

  if(is.null(dt_exp)){
    error_message <- paste0(error_message, "Load experimental data!<br>")
  }
  if(is.null(dt_gt)){
    error_message <- paste0(error_message, "Load genotype data!<br>")
  }

  # --- Check the analytical thresholds --- #

  markers_est_par <- dt_setting_est_par$Marker
  at_all_mk <- dt_setting_est_par$AT
  bool_without_at <- is.na(at_all_mk)
  if(any(bool_without_at)){
    error_message <- paste0(error_message, "Set the analytical threshold(s) of the following marker(s): ", paste0(markers_est_par[bool_without_at], collapse = ", "), ".<br>")
  }

  if(all(!is.null(dt_exp), !is.null(dt_gt))){

    # --- Check the marker set of the experimental data --- #

    sn_mk_problem_exp <- check_dt_mk(dt_exp, markers_est_par) # The function is defined in est_sr_main.R
    if(length(sn_mk_problem_exp) > 0){
      error_message <- paste0(error_message, "Inappropriate marker set of the following sample(s) in the experimental data: ", paste0(sn_mk_problem_exp, collapse = ", "), ".<br>")
    }

    # --- Check the marker set of the genotype data --- #

    sn_mk_problem_gt <- check_dt_mk(dt_gt, markers_est_par) # The function is defined in est_sr_main.R
    if(length(sn_mk_problem_gt) > 0){
      error_message <- paste0(error_message, "Inappropriate marker set of the following sample(s) in the genotype data: ", paste0(sn_mk_problem_gt, collapse = ", "), ".<br>")
    }

    # --- Check the sample names --- #

    sn_exp <- unique(dt_exp$SampleName)
    sn_gt <- unique(dt_gt$SampleName)
    sn_without_gt <- setdiff(sn_exp, sn_gt)
    if(length(sn_without_gt) > 0){
      error_message <- paste0(error_message, "No genotype data of the following samples: ", paste0(sn_without_gt, collapse = ", "), ".<br>")
    }
  }

  return(error_message)
}

#' arrange_exp_data
#'
#' @description The function to arrange the experimental data
#' @param dt_exp The data.table of the experimental data
#' @param dt_gt The data.table of the genotype data
#' @param markers_auto Autosomal markers
#' @param sr_consider_auto Information on consideration of back, forward, double-back, and minus 2-nt stutters in autosomal markers
#' @param at_auto Analytical thresholds for autosomal markers
arrange_exp_data <- function(dt_exp, dt_gt, markers_auto, sr_consider_auto, at_auto){
  peak_list <- height_list <- gt_list <- list()

  name_dt_exp <- names(dt_exp)
  index_exp_peak <- grep("Allele", name_dt_exp)
  index_exp_height <- grep("Height", name_dt_exp)
  sample_names <- unique(dt_exp[, SampleName])

  name_dt_gt <- names(dt_gt)
  index_gt_allele <- grep("Allele", name_dt_gt)

  n_mk_auto <- length(markers_auto)

  setindex(dt_exp, SampleName)
  setindex(dt_gt, SampleName)

  for(i in 1:length(sample_names)){
    dt_exp_one <- dt_exp[SampleName == sample_names[i]]
    dt_gt_one <- dt_gt[SampleName == sample_names[i]]
    setindex(dt_exp_one, Marker)
    setindex(dt_gt_one, Marker)

    # --- Define objects for one sample --- #

    peak_list_one <- height_list_one <- list()
    gt_one <- matrix(NA, nrow = n_mk_auto, ncol = 2)
    rownames(gt_one) <- markers_auto
    colnames(gt_one) <- c("Allele1", "Allele2")

    for(j in 1:n_mk_auto){

      # --- Extract data of one marker --- #

      mk <- markers_auto[j]
      dt_exp_one_mk <- dt_exp_one[Marker == mk]
      peaks <- as.character(dt_exp_one_mk[, index_exp_peak, with = FALSE])
      heights <- as.numeric(dt_exp_one_mk[, index_exp_height, with = FALSE])

      # --- Remove NA --- #

      index_not_na <- which(!is.na(peaks))
      if(length(index_not_na) == 0){
        peak_list_one[[j]] <- character(0)
        height_list_one[[j]] <- numeric(0)
      }else{
        peaks <- peaks[index_not_na]
        heights <- heights[index_not_na]

        # --- Remove peaks below AT --- #

        index_above_at <- which(heights >= at_auto[mk])
        if(length(index_above_at) == 0){
          peak_list_one[[j]] <- character(0)
          height_list_one[[j]] <- numeric(0)
        }else{
          peaks <- peaks[index_above_at]
          heights <- heights[index_above_at]

          # --- Extract the genotype of one marker --- #

          dt_gt_one_mk <- dt_gt_one[Marker == mk]
          gt <- as.character(dt_gt_one_mk[, index_gt_allele, with = FALSE])
          gt_one[j, ] <- gt

          # --- Make detectable peaks in the profile under the genotype --- #

          cand_peak <- make_cand_peak(as.numeric(gt), c(sr_consider_auto[mk, ]))

          index_use <- which(is.element(peaks, as.character(cand_peak)))
          if(length(index_use) > 0){
            peak_list_one[[j]] <- peaks_adopt <- peaks[index_use]
            height_list_one[[j]] <- heights[index_use]
          }else{
            peak_list_one[[j]] <- character(0)
            height_list_one[[j]] <- numeric(0)
          }
        }
      }
    }

    names(peak_list_one) <- names(height_list_one) <- markers_auto
    peak_list[[i]] <- peak_list_one
    height_list[[i]] <- height_list_one
    gt_list[[i]] <- gt_one
  }

  names(peak_list) <- names(height_list) <- names(gt_list) <- sample_names
  return(list("peak_list" = peak_list,
              "height_list" = height_list,
              "gt_list" = gt_list))
}
