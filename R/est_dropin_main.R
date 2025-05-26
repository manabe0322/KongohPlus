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

#' generate_init_dropin
#'
#' @description The function to generate initial values for estimating drop-in parameters
generate_init_dropin <- function(){
  shape_di <- runif(1, min = 0, max = 100)
  rate_di <- runif(1, min = 0, max = 100)
  return(list(list(shape_di = shape_di, rate_di = rate_di)))
}

#' run_mle_dropin
#'
#' @description The function to perform MLE for the drop-in model
#' @param n_dropin The number of observed drop-in peaks
#' @param heights A vector of drop-in peak heights
#' @param model_dropin The stan model for the drop-in peak heights
#' @param n_conv The number of successful convergences
#' @param n_max_init The maximum number of initializations
run_mle_dropin <- function(n_dropin, heights, model_dropin, n_conv, n_max_init){
  # Prepare data list for the drop-in model
  data_list_dropin <- list(N = n_dropin, Height = heights)

  count_all <- 0
  count_conv <- 0
  max_lp <- -Inf
  est_param <- rep(NA, 2)

  # Repeat random initialization and optimization until the condition is satisfied
  repeat{
    count_all <- count_all + 1
    init_vals <- generate_init_dropin()
    fit <- model_dropin$optimize(data_list_dropin, init = init_vals)

    # If convergence is successful, perform the following steps
    if(fit$return_codes() == 0){
      lp <- fit$lp()

      # If the difference between the current likelihood and the maximum likelihood is less than 0.1, the current run is considered convergent.
      if(abs(max_lp - lp) < 0.1){
        count_conv <- count_conv + 1
        if(lp - max_lp > 0){
          max_lp <- lp
          est_param <- fit$mle()
        }
      }

      # Reset convergence count to 1 if a new maximum likelihood is found.
      else if(lp - max_lp > 0){
        count_conv <- 1
        max_lp <- lp
        est_param <- fit$mle()
      }

      # Stop the simulation if convergence is successful for the required number of times.
      if(count_conv == n_conv){
        break
      }
    }

    # Stop the simulation if the maximum number of simulations is reached, treating it as a convergence failure.
    if(count_all == n_max_init){
      est_param <- rep(NA, 2)
      break
    }
  }

  names(est_param) <- c("shape_di", "rate_di")
  return(est_param)
}

#' estimate_dropin_params
#'
#' @description The function to estimate drop-in parameters
estimate_dropin_params <- function(dt_ntc, dt_kit, at_dropin, model_dropin, n_conv, n_max_init){
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
  tmp <- run_mle_dropin(n_dropin, heights, model_dropin, n_conv, n_max_init)
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
