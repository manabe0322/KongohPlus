#' get_ae_data
#'
#' @description The function to get data of the locus-specific amplification efficiency
#' @param height_list A list of peak heights
get_ae_data <- function(height_list){
  n_sample <- length(height_list)
  markers <- names(height_list[[1]])
  n_mk <- length(markers)

  height_ae <- matrix(0, n_sample, n_mk)
  rownames(height_ae) <- names(height_list)
  colnames(height_ae) <- markers

  for(i in 1:n_sample){
    height_ae[i, ] <- sapply(height_list[[i]], sum)
  }

  # Extract used data (i.e, removing AE = 0)
  height_ae <- height_ae[apply(height_ae, 1, prod) > 0, , drop = FALSE]
  n_sample_use <- nrow(height_ae)
  if(n_sample_use > 0){
    ae_data <- matrix(0, n_sample_use, n_mk)
    rownames(ae_data) <- rownames(height_ae)
    colnames(ae_data) <- markers

    # Calculate the locus-specific amplification efficiency
    for(i in 1:n_sample_use){
      ae_data[i, ] <- height_ae[i, ] / mean(height_ae[i, ])
    }
  }else{
    ae_data <- height_ae
  }

  ae_list <- list()
  ae_list[[1]] <- ae_data
  ae_list[[2]] <- height_ae
  names(ae_list) <- c("AE", "Height_AE")
  return(ae_list)
}

#' run_mcmc_ae
#'
#' @description The function to perform the MCMC for the locus-specific amplification efficiency
#' @param model_ae The stan model of the locus-specific amplification efficiency
#' @param ae_data The data of the locus-specific amplification efficiency
#' @param height_ae_data The height data for the locus-specific amplification efficiency
#' @param n_warmup The number of simulations for warm-up
#' @param n_sampling The number of simulations for sampling
#' @param n_chains The number of chains
#' @param n_parallel_chains The number of parallel chains
run_mcmc_ae <- function(model_ae, ae_data, height_ae_data,
                        n_warmup, n_sampling, n_chains, n_parallel_chains){
  n_ae_data <- nrow(ae_data)
  n_mk <- ncol(ae_data)
  median_myu_all_mk <- median_sigma_all_mk <- rhat_myu_all_mk <- rhat_sigma_all_mk <- rep(NA, n_mk)
  traceplot_all_mk <- list()
  for(i in 1:n_mk){
    # Prepare data for MCMC
    data_list_ae <- list(N = n_ae_data, AE = ae_data[, i], Height = height_ae_data[, i])

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
    median_myu_all_mk[i] <- fit_ae_summary$median[which(fit_ae_summary$variable == "myu")]
    median_sigma_all_mk[i] <- fit_ae_summary$median[which(fit_ae_summary$variable == "sigma")]

    # rhat of each parameter
    rhat_myu_all_mk[i] <- fit_ae_summary$rhat[which(fit_ae_summary$variable == "myu")]
    rhat_sigma_all_mk[i] <- fit_ae_summary$rhat[which(fit_ae_summary$variable == "sigma")]

    # Trace plot
    fit_ae$draws() %>% mcmc_trace(pars = c("myu", "sigma")) -> traceplot_all_mk[[i]]
  }

  return(list(median_myu_all_mk, median_sigma_all_mk, rhat_myu_all_mk, rhat_sigma_all_mk, traceplot_all_mk))
}
