#' run_mcmc_ae
#'
#' @description The function to perform the MCMC for the locus-specific amplification efficiency
run_mcmc_ae <- function(model_ae, ae_data, height_ae_data,
                        n_warmup, n_sampling, n_chains, n_parallel_chains){
  n_ae_data <- nrow(ae_data)
  n_mk <- ncol(ae_data)
  for(i in 1:n_mk){
    data_list_ae <- list(N = n_ae_data, AE = ae_data[, i], Height = height_ae_data[, i])

    fit_ae <- model_ae$sample(
      data_list_ae,
      iter_warmup = n_warmup,
      iter_sampling = n_sampling,
      chains = n_chains,
      parallel_chains = n_parallel_chains
    )
  }
}
