#' analyze_ae_main
#'
#' @description The function to analyze AE data
#' @param markers_auto The vector of autosomal markers
#' @param model_ae The stan model for AE
#' @param pae_data The data of the proportion of the AE
#' @param n_warmup The number of simulations for warm-up
#' @param n_sampling The number of simulations for sampling
#' @param n_chains The number of chains
#' @param n_parallel_chains The number of parallel chains
analyze_ae_main <- function(markers_auto, model_ae,
                            pae_data,
                            n_warmup, n_sampling, n_chains, n_parallel_chains){

  cat("MCMC for the locus-specific amplification efficiency\n")
  n_mk <- length(markers_auto)
  data_list_ae <- list(N = ncol(pae_data), N_Mk = n_mk, AE = pae_data)
  tmp <- run_mcmc_ae(model_ae, data_list_ae, n_mk, n_warmup, n_sampling, n_chains, n_parallel_chains)
  median_alpha <- tmp$median_alpha
  rhat_alpha <- tmp$rhat_alpha
  rhat_lp <- tmp$rhat_lp
  traceplot <- tmp$traceplot

  dt_result_ae <- data.table(Marker = markers_auto,
                             Alpha = median_alpha,
                             Rhat_Alpha = rhat_alpha,
                             Rhat_Likelihood = c(rhat_lp, rep(NA, n_mk)))

  return(list(dt_result = dt_result_ae,
              traceplot = traceplot,
              plotdata = pae_data))
}
