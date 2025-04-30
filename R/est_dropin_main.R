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

#' estimate_dropin_params
#'
#' @description The function to estimate drop-in parameters
estimate_dropin_params <- function(dt_ntc, dt_kit, at_dropin, model_dropin, warmup, sampling, chains, n_parallel_chains = detectCores()){
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

  # Perform MCMC
  data_list_dropin <- list(N = n_dropin, Height = heights)
  fit_dropin <- model_dropin$sample(
    data_list_dropin,
    iter_warmup = warmup,
    iter_sampling = sampling,
    chains = chains,
    parallel_chains = n_parallel_chains
  )

  # Extract parameters
  tmp <- as_draws_df(fit_dropin$draws())
  lp_values <- tmp$lp__
  density_lp <- density(lp_values)
  lp_mode_estimated <- density_lp$x[which.max(density_lp$y)]
  index_hpd <- which.min(abs(lp_values - lp_mode_estimated))
  lp_dropin <- lp_values[index_hpd]

  shape_di_values <- tmp$shape_di
  shape_di <- shape_di_values[index_hpd]

  rate_di_values <- tmp$rate_di
  rate_di <- rate_di_values[index_hpd]

  fit_dropin_summary <- fit_dropin$summary()
  variables <- fit_dropin_summary$variable
  rhats <- fit_dropin_summary$rhat
  rhat_lp <- rhats[which(variables == "lp__")]
  rhat_shape_di <- rhats[which(variables == "shape_di")]
  rhat_rate_di <- rhats[which(variables == "rate_di")]

  fit_dropin$draws() %>% mcmc_trace(pars = c("shape_di")) -> traceplot_shape_di
  fit_dropin$draws() %>% mcmc_trace(pars = c("rate_di")) -> traceplot_rate_di

  return(list("heights" = heights,
              "n_sample" = n_sample,
              "n_mk" = n_mk,
              "n_dropin" = n_dropin,
              "freq_dropin" = n_dropin / (n_sample * n_mk),
              "lp_dropin" = lp_dropin,
              "shape_di" = shape_di,
              "rate_di" = rate_di,
              "rhat_lp" = rhat_lp,
              "rhat_shape_di" = rhat_shape_di,
              "rhat_rate_di" = rhat_rate_di,
              "traceplot_shape_di" = traceplot_shape_di,
              "traceplot_rate_di" = traceplot_rate_di))
}
