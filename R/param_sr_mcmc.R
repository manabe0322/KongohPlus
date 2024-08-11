#' run_mcmc_sr_gen
#'
#' @description The function to perform MCMC for the stutter ratio (locus-specific with slope)
#' @param model_sr_gen The stan model of the stutter ratio
#' @param data_list_sr The data list of the stutter ratio
#' @param warmup The number of warm-up iterations for each chain
#' @param sampling The number of sampling iterations for each chain
#' @param chains The number of chains
#' @param n_parallel_chains The number of parallel chains
run_mcmc_sr_gen <- function(model_sr_gen, data_list_sr, warmup, sampling, chains, n_parallel_chains){
  # Perform MCMC
  fit_sr <- model_sr_gen$sample(
    data_list_sr,
    iter_warmup = warmup,
    iter_sampling = sampling,
    chains = chains,
    parallel_chains = n_parallel_chains
  )

  # Summary data
  fit_sr_summary <- fit_sr$summary()
  medians <- fit_sr_summary$median
  variables <- fit_sr_summary$variable
  rhats <- fit_sr_summary$rhat
  median_lp <- medians[which(variables == "lp__")]
  median_slope <- medians[which(variables == "slope")]
  median_intercept <- medians[which(variables == "intercept")]
  median_sigma <- medians[which(variables == "sigma")]
  rhat_slope <- rhats[which(variables == "slope")]
  rhat_intercept <- rhats[which(variables == "intercept")]
  rhat_sigma <- rhats[which(variables == "sigma")]
  rhat_lp <- rhats[which(variables == "lp__")]

  # Trace plot
  fit_sr$draws() %>% mcmc_trace(pars = c("slope", "intercept", "sigma")) -> traceplot

  return(list("median_lp" = median_lp,
              "median_slope" = median_slope,
              "median_intercept" = median_intercept,
              "median_sigma" = median_sigma,
              "rhat_slope" = rhat_slope,
              "rhat_intercept" = rhat_intercept,
              "rhat_sigma" = rhat_sigma,
              "rhat_lp" = rhat_lp,
              "traceplot" = traceplot))
}

#' run_mcmc_sr_uni
#'
#' @description The function to perform MCMC for the stutter ratio (uniform)
#' @param model_sr_uni The stan model of the stutter ratio
#' @param data_list_sr The data list of the stutter ratio
#' @param warmup The number of warm-up iterations for each chain
#' @param sampling The number of sampling iterations for each chain
#' @param chains The number of chains
#' @param n_parallel_chains The number of parallel chains
run_mcmc_sr_uni <- function(model_sr_uni, data_list_sr, warmup, sampling, chains, n_parallel_chains){
  # Perform MCMC
  fit_sr <- model_sr_uni$sample(
    data_list_sr,
    iter_warmup = warmup,
    iter_sampling = sampling,
    chains = chains,
    parallel_chains = n_parallel_chains
  )

  # Summary data
  fit_sr_summary <- fit_sr$summary()
  medians <- fit_sr_summary$median
  variables <- fit_sr_summary$variable
  rhats <- fit_sr_summary$rhat
  median_lp <- medians[which(variables == "lp__")]
  median_myu <- medians[which(variables == "myu")]
  median_sigma <- medians[which(variables == "sigma")]
  rhat_myu <- rhats[which(variables == "myu")]
  rhat_sigma <- rhats[which(variables == "sigma")]
  rhat_lp <- rhats[which(variables == "lp__")]

  # Trace plot
  fit_sr$draws() %>% mcmc_trace(pars = c("myu", "sigma")) -> traceplot

  return(list("median_lp" = median_lp,
              "median_myu" = median_myu,
              "median_sigma" = median_sigma,
              "rhat_myu" = rhat_myu,
              "rhat_sigma" = rhat_sigma,
              "rhat_lp" = rhat_lp,
              "traceplot" = traceplot))
}
