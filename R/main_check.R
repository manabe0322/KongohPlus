check_main_file <- function(dt_csp, dt_ref, dt_af, dt_kit, method = "deconvo"){
  error_message <- ""

  # Whether all required data is loaded or not
  if(is.null(dt_csp)){
    error_message <- paste0(error_message, "Load the crime stain profile!<br>")
  }
  if(method == "lr"){
    if(is.null(dt_ref)){
      error_message <- paste0(error_message, "Load the reference data!<br>")
    }
  }
  if(is.null(dt_af)){
    error_message <- paste0(error_message, "Load the allele frequencies!<br>")
  }

  if(all(!is.null(dt_csp), !is.null(dt_ref), !is.null(dt_af))){
    # Extract markers
    setkey(dt_kit, Sex_Chr)
    dt_kit_auto <- dt_kit[""]
    markers_kit <- dt_kit_auto[, Marker]

    sn_all_csp <- unique(dt_csp[, SampleName])
    sn_all_ref <- unique(dt_ref[, SampleName])
  }
}
