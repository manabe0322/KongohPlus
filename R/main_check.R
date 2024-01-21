#' check_csp_alcor
#'
#' @description The function to check whether the corrected alleles for all peaks are registered or not
#' @param dt_csp The data.table of the crime stain profiles
#' @param dt_at The data.table of the analytical threshold
#' @param dt_alcor The data.table of the corrected alleles
#' @param markers_auto A vector of the autosomal markers
check_csp_alcor <- function(dt_csp, dt_at, dt_alcor, markers_auto){
  error_message_alcor <- ""

  names_csp <- names(dt_csp)
  index_allele <- grep("Allele", names_csp)
  index_height <- grep("Height", names_csp)

  n_mk <- length(markers_auto)

  setkey(dt_csp, SampleName)
  sn_all_csp <- unique(dt_csp[, SampleName])

  setkey(dt_at, Marker)
  setkey(dt_alcor, Marker)

  for(i in 1:length(sn_all_csp)){
    sn1_csp <- sn_all_csp[i]
    dt_csp_sn1 <- dt_csp[sn1_csp]
    setkey(dt_csp_sn1, Marker)
    for(j in 1:n_mk){
      mk <- markers_auto[j]
      dt_csp_sn1_mk1 <- dt_csp_sn1[mk]
      peaks <- as.numeric(dt_csp_sn1_mk1[, index_allele, with = FALSE])
      peaks <- peaks[!is.na(peaks)]
      heights <- as.numeric(dt_csp_sn1_mk1[, index_height, with = FALSE])
      heights <- heights[!is.na(heights)]

      dt_at_mk1 <- dt_at[mk]
      at <- as.numeric(dt_at_mk1[, AT])
      peaks <- peaks[heights >= at]

      dt_alcor_mk1 <- dt_alcor[mk]
      sr_models <- unique(c(dt_alcor_mk1[, Model_BSR], dt_alcor_mk1[, Model_FSR], dt_alcor_mk1[, Model_DSR]))
      if(any(sr_models != "Not consider")){
        alleles <- as.numeric(dt_alcor_mk1[, Allele])
        bool_not_register <- !is.element(peaks, alleles)
        if(any(bool_not_register)){
          error_message_alcor <- paste0(error_message_alcor, "  - ", mk, ": ", paste(peaks[bool_not_register], collapse = ", "), "<br>")
        }
      }
    }
  }
  return(error_message_alcor)
}

#' check_main_file
#'
#' @description The function to check the main files
#' @param dt_csp The data.table of the crime stain profiles
#' @param dt_ref The data.table of the reference profiles
#' @param dt_af The data.table of the allele frequencies
#' @param dt_kit The data.table of the kit information
#' @param dt_at The data.table of the analytical threshold
#' @param dt_alcor The data.table of the corrected alleles
#' @param method 'deconvo' (deconvolution) or 'lr' (likelihood ratio)
check_main_file <- function(dt_csp, dt_ref, dt_af, dt_kit, dt_at, dt_alcor, method = "deconvo"){
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
    # Extract autosomal markers
    dt_kit_auto <- extract_kit_auto(dt_kit)
    markers_auto <- unique(dt_kit_auto[, Marker])

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

    # Check whether the corrected alleles for all peaks are registered or not
    if(error_message == ""){
      error_message_alcor <- check_csp_alcor(dt_csp, dt_at, dt_alcor, markers_auto)
      if(error_message_alcor != ""){
        error_message <- paste0(error_message, "Average stutter ratios of the following alleles are not registered!<br>", error_message_alcor)
      }
    }
  }

  return(error_message)
}
