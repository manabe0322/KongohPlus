#' analyze_ae_main
#'
#' @description The function to analyze AE data
#' @param markers_par The vector of markers for parameters
#' @param model_ae The stan model for AE
#' @param ae_data The AE data of all markers
#' @param height_ae_data The height data for AE of all markers
#' @param n_warmup The number of simulations for warm-up
#' @param n_sampling The number of simulations for sampling
#' @param n_chains The number of chains
#' @param n_parallel_chains The number of parallel chains
analyze_ae_main <- function(markers_par, model_ae,
                            ae_data, height_ae_data,
                            n_warmup, n_sampling, n_chains, n_parallel_chains){
  # Define objects
  n_mk <- length(markers_par)
  n_ae_data <- nrow(ae_data)
  ae_median_myu_all_mk <- ae_median_sigma_all_mk <- ae_rhat_myu_all_mk <- ae_rhat_sigma_all_mk <- ae_rhat_lp_all_mk <- rep(NA, n_mk)
  ae_trace_all_mk <- vector(mode = "list", length = n_mk)
  for(i in 1:n_mk){
    cat(paste0("MCMC for the locus-specific amplification efficiency of the ", markers_par[i], "\n"))
    data_list_ae <- list(N = n_ae_data, AE = ae_data[, i], Height = height_ae_data[, i])
    tmp <- run_mcmc_ae(model_ae, data_list_ae, n_warmup, n_sampling, n_chains, n_parallel_chains)
    ae_median_myu_all_mk[i] <- tmp$median_myu
    ae_median_sigma_all_mk[i] <- tmp$median_sigma
    ae_rhat_myu_all_mk[i] <- tmp$rhat_myu
    ae_rhat_sigma_all_mk[i] <- tmp$rhat_sigma
    ae_rhat_lp_all_mk[i] <- tmp$rhat_lp
    ae_trace_all_mk[[i]] <- tmp$traceplot
  }

  dt_result_ae <- data.table(Marker = markers_par,
                             Myu = ae_median_myu_all_mk,
                             Sigma = ae_median_sigma_all_mk,
                             Rhat_Myu = ae_rhat_myu_all_mk,
                             Rhat_Sigma = ae_rhat_sigma_all_mk,
                             Rhat_Likelihood = ae_rhat_lp_all_mk)

  return(list(dt_result = dt_result_ae,
              traceplot = ae_trace_all_mk,
              plotdata = ae_data))
}
