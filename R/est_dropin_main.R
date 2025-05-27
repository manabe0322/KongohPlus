#' check_est_dropin
#'
#' @description The function to check files for estimating dropin parameters
#' @param dt_ntc The data.table of the non template control data
#' @param dt_kit The data.table for information on the kit
check_est_dropin <- function(dt_ntc, dt_kit){
  error_message <- ""

  # Whether all required data is loaded or not
  if(is.null(dt_ntc)){
    error_message <- paste0(error_message, "Load non template control data!<br>")
  }

  if(!is.null(dt_ntc)){
    # Check the marker set of the experimental data
    sn_mk_problem_dropin <- check_dt_mk(dt_ntc, unique(dt_kit$Marker)) # The function is defined in est_sr_main.R
    if(length(sn_mk_problem_dropin) > 0){
      error_message <- paste0(error_message, "Inappropriate marker set of the following sample(s) in the non template control data: ", paste0(sn_mk_problem_dropin, collapse = ", "), ".<br>")
    }
  }

  return(error_message)
}

#' generate_init_sobol_dropin
#'
#' @description The function to generate initial values based on Sobol sequence for the drop-in model
#' @param x Sobol point
#' @param min_shape_di The minimum shape
#' @param max_shape_di The maximum shape
#' @param min_rate_di The minimum rate
#' @param max_rate_di The maximum rate
generate_init_sobol_dropin <- function(x, min_shape_di, max_shape_di, min_rate_di, max_rate_di){
  # Get only one Sobol point (x-th)
  sobol_vec <- sobol(n = x, d = 2, randomize = "digital.shift")[x, ]

  # Index partitioning and scaling
  shape_di <- sobol_vec[1] * (max_shape_di - min_shape_di) + min_shape_di
  rate_di <- sobol_vec[2] * (max_rate_di - min_rate_di) + min_rate_di

  return(list(list(shape_di = shape_di, rate_di = rate_di)))
}

#' run_mle_dropin
#'
#' @description The function to perform MLE for the drop-in model
#' @param n_dropin The number of observed drop-in peaks
#' @param heights A vector of drop-in peak heights
#' @param model_dropin The stan model for the drop-in peak heights
#' @param min_shape_di The minimum shape
#' @param max_shape_di The maximum shape
#' @param min_rate_di The minimum rate
#' @param max_rate_di The maximum rate
#' @param initial_point_factor This factor determines the number of Sobol points: N = initial_point_factor × number of parameters
run_mle_dropin <- function(n_dropin, heights, model_dropin, min_shape_di, max_shape_di, min_rate_di, max_rate_di, initial_point_factor){
  # Prepare data list for the drop-in model
  data_list_dropin <- list(N = n_dropin, Height = heights)

  n_sobol <- initial_point_factor * 2
  init_list <- lapply(1:n_sobol, function(x) generate_init_sobol_dropin(x, min_shape_di, max_shape_di, min_rate_di, max_rate_di))

  tmp_param <- matrix(NA, n_sobol, 3)
  for(i in 1:n_sobol){
    fit <- model_dropin$optimize(data_list_dropin, init = init_list[[i]], refresh = 0, show_messages = FALSE, show_exceptions = FALSE)
    if(fit$return_codes() == 0){
      tmp <- fit$draws()
      tmp_param[i, ] <- tmp[1, ]
    }
  }

  index_max_lp <- which.max(tmp_param[, 1])
  if(length(index_max_lp) == 1){
    est_param <- tmp_param[index_max_lp, -1]
  }else{
    est_param <- rep(NA, 2)
  }
  names(est_param) <- c("shape_di", "rate_di")

  return(est_param)
}

#' estimate_dropin_params
#'
#' @description The function to estimate drop-in parameters
estimate_dropin_params <- function(dt_ntc, dt_kit, at_dropin, model_dropin, min_shape_di, max_shape_di, min_rate_di, max_rate_di, initial_point_factor){
  # Extract autosomal markers
  dt_kit_auto <- extract_kit_auto(dt_kit)
  markers_auto <- unique(dt_kit_auto[, Marker])
  n_mk <- length(markers_auto)

  # Extract drop-in data
  dt_ntc_auto <- dt_ntc[which(is.element(dt_ntc$Marker, markers_auto)), ]
  col_ntc <- names(dt_ntc_auto)
  heights <- as.numeric(unlist(dt_ntc_auto[, grep("Height", col_ntc), with = FALSE]))
  heights <- heights[!is.na(heights)]
  heights <- heights[heights >= at_dropin]
  n_dropin <- length(heights)
  n_sample <- length(unique(dt_ntc$SampleName))

  # Perform MLE
  tmp <- run_mle_dropin(n_dropin, heights, model_dropin, min_shape_di, max_shape_di, min_rate_di, max_rate_di, initial_point_factor)
  shape_di <- tmp["shape_di"]
  rate_di <- tmp["rate_di"]

  return(list("heights" = heights,
              "n_sample" = n_sample,
              "n_mk" = n_mk,
              "n_dropin" = n_dropin,
              "freq_dropin" = n_dropin / (n_sample * n_mk),
              "shape_di" = shape_di,
              "rate_di" = rate_di))
}
