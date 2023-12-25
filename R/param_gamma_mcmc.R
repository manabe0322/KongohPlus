#' run_mcmc_gamma
#'
#' @description The function to perform MCMC for the gamma distribution
#' @param model_gamma The stan model of the gamma distribution
#' @param data_list_gamma The data list of the gamma distribution
#' @param n_warmup The number of simulations for warm-up
#' @param n_sampling The number of simulations for sampling
#' @param n_chains The number of chains
#' @param n_parallel_chains The number of parallel chains
run_mcmc_gamma <- function(model_gamma, data_list_gamma, n_warmup, n_sampling, n_chains, n_parallel_chains){
  # Perform MCMC
  fit_gamma <- model_gamma$sample(
    data_list_gamma,
    iter_warmup = n_warmup,
    iter_sampling = n_sampling,
    chains = n_chains,
    parallel_chains = n_parallel_chains
  )

  # Summary data
  fit_gamma_summary <- fit_gamma$summary()

  # Median of the parameter 'beta'
  median_beta <- fit_gamma_summary$median[which(fit_gamma_summary$variable == "beta")]

  # rhat of the parameter 'beta' and the likelihood
  rhat_beta <- fit_gamma_summary$rhat[which(fit_gamma_summary$variable == "beta")]
  rhat_lp <- fit_gamma_summary$rhat[which(fit_gamma_summary$variable == "lp__")]

  # Trace plot
  fit_gamma$draws() %>% mcmc_trace(pars = c("beta")) -> traceplot

  return(list("median_beta" = median_beta,
              "rhat_beta" = rhat_beta,
              "rhat_lp" = rhat_lp,
              "traceplot" = traceplot))
}
