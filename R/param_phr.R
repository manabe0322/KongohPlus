#' get_phr_data
#'
#' @description The function to prepare data of the peak height ratio
#' @param peak_list The list of peaks
#' @param height_list The list of peak heights
#' @param gt_list The list of genotypes
get_phr_data <- function(peak_list, height_list, gt_list){
  n_sample <- length(peak_list)
  n_mk <- length(peak_list[[1]])
  phr_all <- mean_ph_all <- matrix(NA, n_mk, n_sample)

  for(i in 1:n_sample){
    peak_list_one <- peak_list[[i]]
    height_list_one <- height_list[[i]]
    gt_list_one <- gt_list[[i]]

    for(j in 1:n_mk){
      gt <- gt_list_one[j, ]
      if(gt[1] != gt[2] && abs(round(as.numeric(gt[1]) - as.numeric(gt[2]), 1)) != 1){
        peaks <- peak_list_one[[j]]
        heights <- height_list_one[[j]]

        index_al <- which(is.element(peaks, gt))
        if(length(index_al) == 2){
          oph <- heights[index_al]
          current_index <- n_mk * (i - 1) + j
          phr_all[j, i] <- oph[2] / oph[1]
          mean_ph_all[j, i] <- mean(oph)
        }
      }
    }
  }

  return(list("phr_all" = phr_all,
              "mean_ph_all" = mean_ph_all))
}

analyze_phr_main <- function(markers_auto, model_phr, phr_all, mean_ph_all,
                             n_warmup, n_sampling, n_chains, n_parallel_chains){

  cat("MCMC for the peak height ratio\n")
  n_mk <- length(markers_auto)
  phr_myu_median_all_mk <- rep(NA, n_mk)
  phr_sigma_median_all_mk <- rep(NA, n_mk)
  phr_lp_rhat_all_mk <- rep(NA, n_mk)
  phr_myu_rhat_all_mk <- rep(NA, n_mk)
  phr_sigma_rhat_all_mk <- rep(NA, n_mk)
  traceplot_phr_all_mk <- vector(mode = "list", length = n_mk)

  for(i in 1:n_mk){
    # Extract PHR data of one marker
    phr_mk <- phr_all[i, ]
    mean_ph_mk <- mean_ph_all[i, ]

    # Remove NA data
    index_phr <- which(!is.na(phr_mk))
    phr_mk <- phr_mk[index_phr]
    mean_ph_mk <- mean_ph_mk[index_phr]

    # Perform MCMC
    data_list_phr <- list(N = length(phr_mk), PHR = phr_mk, Mean_Ph = mean_ph_mk)
    tmp <- run_mcmc_phr(model_phr, data_list_phr, n_warmup, n_sampling, n_chains, n_parallel_chains)
    phr_myu_median_all_mk[i] <- tmp$median_myu
    phr_sigma_median_all_mk[i] <- tmp$median_sigma
    phr_lp_rhat_all_mk[i] <- tmp$rhat_lp
    phr_myu_rhat_all_mk[i] <- tmp$rhat_myu
    phr_sigma_rhat_all_mk[i] <- tmp$rhat_sigma
    traceplot_phr_all_mk[[i]] <- tmp$traceplot_phr
  }

  dt_result_phr <- data.table(Marker = markers_auto,
                              Myu = phr_myu_median_all_mk,
                              Sigma = phr_sigma_median_all_mk,
                              Rhat_Myu = phr_myu_rhat_all_mk,
                              Rhat_Sigma = phr_sigma_rhat_all_mk,
                              Rhat_Lp = phr_lp_rhat_all_mk)

  return(list("dt_result_phr" = dt_result_phr,
              "traceplot_phr_all_mk" = traceplot_phr_all_mk))
}

run_mcmc_phr <- function(model_phr, data_list_phr, n_warmup, n_sampling, n_chains, n_parallel_chains){
  # Perform MCMC
  fit_phr <- model_phr$sample(
    data_list_phr,
    iter_warmup = n_warmup,
    iter_sampling = n_sampling,
    chains = n_chains,
    parallel_chains = n_parallel_chains
  )

  # Summary data
  fit_phr_summary <- fit_phr$summary()
  variables <- fit_phr_summary$variable
  medians <- fit_phr_summary$median
  rhats <- fit_phr_summary$rhat
  median_myu <- medians[which(variables == "myu")]
  median_sigma <- medians[which(variables == "sigma")]
  rhat_lp <- rhats[which(variables == "lp__")]
  rhat_myu <- rhats[which(variables == "myu")]
  rhat_sigma <- rhats[which(variables == "sigma")]

  # Trace plot
  fit_phr$draws() %>% mcmc_trace(pars = c("myu", "sigma")) -> traceplot_phr

  return(list("median_myu" = median_myu,
              "median_sigma" = median_sigma,
              "rhat_lp" = rhat_lp,
              "rhat_myu" = rhat_myu,
              "rhat_sigma" = rhat_sigma,
              "traceplot_phr" = traceplot_phr))
}
