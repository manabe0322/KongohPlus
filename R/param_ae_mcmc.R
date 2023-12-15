#' run_mcmc_ae
#'
#' @description The function to perform the MCMC for the locus-specific amplification efficiency
#' @param model_ae The stan model of the locus-specific amplification efficiency
#' @param data_list_ae The data list of the locus-specific amplification efficiency
#' @param n_warmup The number of simulations for warm-up
#' @param n_sampling The number of simulations for sampling
#' @param n_chains The number of chains
#' @param n_parallel_chains The number of parallel chains
run_mcmc_ae <- function(model_ae, data_list_ae, n_warmup, n_sampling, n_chains, n_parallel_chains){
  # Perform MCMC
  fit_ae <- model_ae$sample(
    data_list_ae,
    iter_warmup = n_warmup,
    iter_sampling = n_sampling,
    chains = n_chains,
    parallel_chains = n_parallel_chains
  )

  # Summary data
  fit_ae_summary <- fit_ae$summary()

  # Median of each parameter
  median_myu <- fit_ae_summary$median[which(fit_ae_summary$variable == "myu")]
  median_sigma <- fit_ae_summary$median[which(fit_ae_summary$variable == "sigma")]

  # rhat of each parameter
  rhat_myu <- fit_ae_summary$rhat[which(fit_ae_summary$variable == "myu")]
  rhat_sigma <- fit_ae_summary$rhat[which(fit_ae_summary$variable == "sigma")]

  # Trace plot
  fit_ae$draws() %>% mcmc_trace(pars = c("myu", "sigma")) -> traceplot

  return(list("median_myu" = median_myu, "median_sigma" = median_sigma, "rhat_myu" = rhat_myu, "rhat_sigma" = rhat_sigma, "traceplot" = traceplot))
}
