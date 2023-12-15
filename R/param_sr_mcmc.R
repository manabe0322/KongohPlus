#' run_mcmc_sr_gen
#'
#' @description The function to perform MCMC for the stutter ratio (locus-specific with slope)
#' @param model_sr_gen The stan model of the stutter ratio
#' @param data_list_sr The data list of the stutter ratio
#' @param n_warmup The number of simulations for warm-up
#' @param n_sampling The number of simulations for sampling
#' @param n_chains The number of chains
#' @param n_parallel_chains The number of parallel chains
run_mcmc_sr_gen <- function(model_sr_gen, data_list_sr, n_warmup, n_sampling, n_chains, n_parallel_chains){
  # Perform MCMC
  fit_sr <- model_sr_gen$sample(
    data_list_sr,
    iter_warmup = n_warmup,
    iter_sampling = n_sampling,
    chains = n_chains,
    parallel_chains = n_parallel_chains
  )

  # Summary data
  fit_sr_summary <- fit_sr$summary()

  # Median of the likelihood
  median_loglike <- fit_sr_summary$median[which(fit_sr_summary$variable == "lp__")]

  # Median of each parameter
  median_slope <- fit_sr_summary$median[which(fit_sr_summary$variable == "slope")]
  median_intercept <- fit_sr_summary$median[which(fit_sr_summary$variable == "intercept")]
  median_sigma <- fit_sr_summary$median[which(fit_sr_summary$variable == "sigma")]

  # rhat of each parameter
  rhat_slope <- fit_sr_summary$rhat[which(fit_sr_summary$variable == "slope")]
  rhat_intercept <- fit_sr_summary$rhat[which(fit_sr_summary$variable == "intercept")]
  rhat_sigma <- fit_sr_summary$rhat[which(fit_sr_summary$variable == "sigma")]

  # Trace plot
  fit_sr$draws() %>% mcmc_trace(pars = c("myu", "sigma")) -> traceplot

  return(list("median_loglike" = median_loglike,
              "median_slope" = median_slope, "median_intercept" = median_intercept, "median_sigma" = median_sigma,
              "rhat_slope" = rhat_slope, "rhat_intercept" = rhat_intercept, "rhat_sigma" = rhat_sigma,
              "traceplot" = traceplot))
}

#' run_mcmc_sr_together
#'
#' @description The function to perform MCMC for the stutter ratio (multiple loci together)
#' @param model_sr_together The stan model of the stutter ratio
#' @param data_list_sr The data list of the stutter ratio
#' @param n_warmup The number of simulations for warm-up
#' @param n_sampling The number of simulations for sampling
#' @param n_chains The number of chains
#' @param n_parallel_chains The number of parallel chains
run_mcmc_sr_together <- function(model_sr_together, data_list_sr, n_warmup, n_sampling, n_chains, n_parallel_chains){
  # Perform MCMC
  fit_sr <- model_sr_together$sample(
    data_list_sr,
    iter_warmup = n_warmup,
    iter_sampling = n_sampling,
    chains = n_chains,
    parallel_chains = n_parallel_chains
  )

  # Summary data
  fit_sr_summary <- fit_sr$summary()

  # Median of the likelihood
  median_loglike <- fit_sr_summary$median[which(fit_sr_summary$variable == "lp__")]

  # Median of each parameter
  n_mk_together <- model_sr_together$N
  median_myus <- rep(NA, n_mk_together)
  for(i in 1:n_mk_together){
    median_myus[i] <- fit_sr_summary$median[which(fit_sr_summary$variable == paste0("myu[", i, "]"))]
  }
  median_sigma <- fit_sr_summary$median[which(fit_sr_summary$variable == "sigma")]

  # rhat of each parameter
  rhat_myus <- rep(NA, n_mk_together)
  for(i in 1:n_mk_together){
    rhat_myus[i] <- fit_sr_summary$rhat[which(fit_sr_summary$variable == paste0("myu[", i, "]"))]
  }
  rhat_sigma <- fit_sr_summary$rhat[which(fit_sr_summary$variable == "sigma")]

  # Trace plot
  fit_sr$draws() %>% mcmc_trace(pars = c(paste0("myu[", 1:n_mk_together, "]"), "sigma")) -> traceplot

  return(list("median_loglike" = median_loglike,
              "median_myus" = median_myus, "median_sigma" = median_sigma,
              "rhat_myus" = rhat_myus, "rhat_sigma" = rhat_sigma,
              "traceplot" = traceplot))
}

#' run_mcmc_sr_uni
#'
#' @description The function to perform MCMC for the stutter ratio (uniform)
#' @param model_sr_uni The stan model of the stutter ratio
#' @param data_list_sr The data list of the stutter ratio
#' @param n_warmup The number of simulations for warm-up
#' @param n_sampling The number of simulations for sampling
#' @param n_chains The number of chains
#' @param n_parallel_chains The number of parallel chains
run_mcmc_sr_uni <- function(model_sr_uni, data_list_sr, n_warmup, n_sampling, n_chains, n_parallel_chains){
  # Perform MCMC
  fit_sr <- model_sr_uni$sample(
    data_list_sr,
    iter_warmup = n_warmup,
    iter_sampling = n_sampling,
    chains = n_chains,
    parallel_chains = n_parallel_chains
  )

  # Summary data
  fit_sr_summary <- fit_sr$summary()

  # Median of the likelihood
  median_loglike <- fit_sr_summary$median[which(fit_sr_summary$variable == "lp__")]

  # Median of each parameter
  median_myu <- fit_sr_summary$median[which(fit_sr_summary$variable == "myu")]
  median_sigma <- fit_sr_summary$median[which(fit_sr_summary$variable == "sigma")]

  # rhat of each parameter
  rhat_myu <- fit_sr_summary$rhat[which(fit_sr_summary$variable == "myu")]
  rhat_sigma <- fit_sr_summary$rhat[which(fit_sr_summary$variable == "sigma")]

  # Trace plot
  fit_sr$draws() %>% mcmc_trace(pars = c("myu", "sigma")) -> traceplot

  return(list("median_loglike" = median_loglike,
              "median_myu" = median_myu, "median_sigma" = median_sigma,
              "rhat_myu" = rhat_myu, "rhat_sigma" = rhat_sigma,
              "traceplot" = traceplot))
}
