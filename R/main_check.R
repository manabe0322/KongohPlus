#' check_main_file
#'
#' @description The function to check the main files
#' @param dt_csp The data.table of the crime stain profiles
#' @param dt_ref The data.table of the reference profiles
#' @param dt_af The data.table of the allele frequencies
#' @param markers_auto A character vector of the autosomal markers
#' @param method 'deconvo' (deconvolution) or 'lr' (likelihood ratio)
check_main_file <- function(dt_csp, dt_ref, dt_af, markers_auto, method){
  error_message <- ""
  bool_load_all <- TRUE

  # Whether all required data is loaded or not
  if(is.null(dt_csp)){
    error_message <- paste0(error_message, "Load the crime stain profile!<br>")
    bool_load_all <- FALSE
  }
  if(method == "lr"){
    if(is.null(dt_ref)){
      error_message <- paste0(error_message, "Load the reference data!<br>")
      bool_load_all <- FALSE
    }
  }
  if(is.null(dt_af)){
    error_message <- paste0(error_message, "Load the allele frequencies!<br>")
    bool_load_all <- FALSE
  }

  if(bool_load_all){
    setkey(dt_csp, SampleName)
    sn_all_csp <- unique(dt_csp[, SampleName])
    for(i in 1:length(sn_all_csp)){
      sn1_csp <- sn_all_csp[i]
      dt_csp_sn1 <- dt_csp[sn1_csp]
      if(!all(is.element(markers_auto, dt_csp_sn1[, Marker]))){
        error_message <- paste0(error_message, "Marker set of the crime stain profile '", sn1_csp, "' is not correct!<br>")
      }
    }

    if(method == "lr"){
      setkey(dt_ref, SampleName)
      sn_all_ref <- unique(dt_ref[, SampleName])
      for(i in 1:length(sn_all_ref)){
        sn1_ref <- sn_all_ref[i]
        dt_ref_sn1 <- dt_ref[sn1_ref]
        if(!all(is.element(markers_auto, dt_ref_sn1[, Marker]))){
          error_message <- paste0(error_message, "Marker set of the reference profile '", sn1_ref, "' is not correct!<br>")
        }
      }
    }
  }

  return(error_message)
}

#' check_csp_alcor
#'
#' @description The function to check whether the corrected alleles for all peaks are registered or not
#' @param dt_csp The data.table of the crime stain profiles
#' @param dt_info_mk The data.table of the information on the ignored marker and the analytical thresholds for each marker
#' @param dt_alcor_combined The data.table of the corrected alleles
check_csp_alcor <- function(dt_csp, dt_info_mk, dt_alcor_combined){
  message_alcor <- ""

  markers_use <- dt_info_mk$Marker[!dt_info_mk$Ignore]
  n_mk <- length(markers_use)

  if(n_mk == 0){
    message_alcor <- "error"
  }else{
    names_csp <- names(dt_csp)
    index_allele <- grep("Allele", names_csp)
    index_height <- grep("Height", names_csp)

    setkey(dt_csp, SampleName)
    sn_all_csp <- unique(dt_csp[, SampleName])

    setkey(dt_info_mk, Marker)
    setkey(dt_alcor_combined, Marker)

    n_warning_mk <- 0
    for(i in 1:length(sn_all_csp)){
      sn1_csp <- sn_all_csp[i]
      dt_csp_sn1 <- dt_csp[sn1_csp]
      setkey(dt_csp_sn1, Marker)
      for(j in 1:n_mk){
        mk <- markers_use[j]
        dt_csp_sn1_mk1 <- dt_csp_sn1[mk]
        peaks <- as.character(dt_csp_sn1_mk1[, index_allele, with = FALSE])
        heights <- as.numeric(dt_csp_sn1_mk1[, index_height, with = FALSE])
        peaks <- peaks[!is.na(heights)]
        heights <- heights[!is.na(heights)]

        dt_info_mk1 <- dt_info_mk[mk]
        at <- as.numeric(dt_info_mk1[, AT])
        peaks <- peaks[heights >= at]

        dt_alcor_combined_mk1 <- dt_alcor_combined[mk]
        alleles_registerd <- dt_alcor_combined_mk1$Allele

        if(!all(is.element(peaks, alleles_registerd))){
          n_warning_mk <- n_warning_mk + 1
          if(n_warning_mk == 1){
            message_alcor <- "Stutter ratios of the following alleles are regarded as zero.<br>"
          }
          message_alcor <- paste0(message_alcor, "  - ", mk, ": ", paste(peaks[!is.element(peaks, alleles_registerd)], collapse = ", "), "<br>")
        }
      }
    }
  }

  return(message_alcor)
}
