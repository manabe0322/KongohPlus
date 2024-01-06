#' check_main_file
#'
#' @description The function to check the main files
#' @param dt_csp The data.table of the crime stain profiles
#' @param dt_ref The data.table of the reference profiles
#' @param dt_af The data.table of the allele frequencies
#' @param dt_kit The data.table of the kit information
#' @param method 'deconvo' (deconvolution) or 'lr' (likelihood ratio)
check_main_file <- function(dt_csp, dt_ref, dt_af, dt_kit, method = "deconvo"){
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
    # Extract markers
    setkey(dt_kit, Sex_Chr)
    dt_kit_auto <- dt_kit[""]
    markers_kit <- dt_kit_auto[, Marker]

    setkey(dt_csp, SampleName)
    sn_all_csp <- unique(dt_csp[, SampleName])
    for(i in 1:length(sn_all_csp)){
      sn1 <- sn_all_csp[i]
      dt_csp_sn1 <- dt_csp[sn1]
      if(!all(is.element(markers_kit, dt_csp_sn1[, Marker]))){
        error_message <- paste0(error_message, "Marker set of the crime stain profile '", sn1, "' is not correct!<br>")
      }
    }

    if(method == "lr"){
      setkey(dt_ref, SampleName)
      sn_all_ref <- unique(dt_ref[, SampleName])
      for(i in 1:length(sn_all_ref)){
        sn1 <- sn_all_ref[i]
        dt_ref_sn1 <- dt_ref[sn1]
        if(!all(is.element(markers_kit, dt_ref_sn1[, Marker]))){
          error_message <- paste0(error_message, "Marker set of the reference profile '", sn1, "' is not correct!<br>")
        }
      }
    }
  }

  return(error_message)
}
