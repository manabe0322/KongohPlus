#' get_pae_data
#'
#' @description The function to get the proportion of the locus-specific amplification efficiency
#' @param peak_list A list of peaks
#' @param height_list A list of peak heights
#' @param gt_list A list of genotypes
get_pae_data <- function(peak_list, height_list, gt_list){
  n_sample <- length(height_list)
  markers <- names(height_list[[1]])
  n_mk <- length(markers)

  pae_data <- matrix(0, n_mk, n_sample)
  for(i in 1:n_sample){
    peak_list_s1 <- peak_list[[i]]
    height_list_s1 <- height_list[[i]]
    gt_list_s1 <- gt_list[[i]]
    judge_detect <- TRUE
    for(j in 1:n_mk){
      peaks <- peak_list_s1[[j]]
      heights <- height_list_s1[[j]]
      gt <- unique(gt_list_s1[j, ])
      if(!all(is.element(gt, peaks))){
        judge_detect <- FALSE
        break
      }
    }

    if(judge_detect){
      tph_per_mk <- sapply(height_list_s1, sum)
      sum_tph <- sum(tph_per_mk)
      pae_data[, i] <- tph_per_mk / sum_tph
    }
  }

  return(pae_data)
}

#' calc_tph_all_data
#'
#' @description The function to calculate the average heterozygote peak heights of all marker
#' @param height_list The list of peak heights
calc_tph_all_data <- function(height_list){
  n_data <- length(height_list)
  tph_all_data <- rep(0, n_data)
  for(i in 1:n_data){
    height_list_one_data <- height_list[[i]]
    tph_per_mk <- sapply(height_list_one_data, sum)
    tph_all_data[i] <- mean(tph_per_mk) / 2
  }
  return(tph_all_data)
}

#' analyze_ae_main
#'
#' @description The function to analyze AE data
#' @param markers_auto The vector of autosomal markers
#' @param model_ae The stan model for AE
#' @param pae_data_use The data of the proportion of the AE which is used for MCMC
#' @param n_warmup The number of simulations for warm-up
#' @param n_sampling The number of simulations for sampling
#' @param n_chains The number of chains
#' @param n_parallel_chains The number of parallel chains
analyze_ae_main <- function(markers_auto, model_ae, pae_data_use,
                            n_warmup, n_sampling, n_chains, n_parallel_chains){

  cat("MCMC for the locus-specific amplification efficiency\n")
  n_mk <- length(markers_auto)
  data_list_ae <- list(N = ncol(pae_data_use), N_Mk = nrow(pae_data_use), PAE = pae_data_use)
  tmp <- run_mcmc_ae(model_ae, data_list_ae, n_mk, n_warmup, n_sampling, n_chains, n_parallel_chains)
  median_alpha <- tmp$median_alpha
  rhat_alpha <- tmp$rhat_alpha
  rhat_lp <- tmp$rhat_lp
  traceplot_ae <- tmp$traceplot_ae

  dt_result_ae <- data.table(Marker = markers_auto,
                             Alpha = median_alpha,
                             Rhat_Alpha = rhat_alpha,
                             Rhat_Lp = c(rhat_lp, rep(NA, n_mk - 1)))

  return(list("dt_result_ae" = dt_result_ae,
              "traceplot_ae" = traceplot_ae))
}

#' run_mcmc_ae
#'
#' @description The function to perform the MCMC for the locus-specific amplification efficiency
#' @param model_ae The stan model of the locus-specific amplification efficiency
#' @param data_list_ae The data list of the locus-specific amplification efficiency
#' @param n_mk The number of markers
#' @param n_warmup The number of simulations for warm-up
#' @param n_sampling The number of simulations for sampling
#' @param n_chains The number of chains
#' @param n_parallel_chains The number of parallel chains
run_mcmc_ae <- function(model_ae, data_list_ae, n_mk, n_warmup, n_sampling, n_chains, n_parallel_chains){
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
  variables <- fit_ae_summary$variable
  medians <- fit_ae_summary$median
  rhats <- fit_ae_summary$rhat
  median_lp <- medians[which(variables == "lp__")]
  median_alpha <- medians[which(is.element(variables, paste0("alpha[", 1:n_mk, "]")))]
  rhat_lp <- rhats[which(variables == "lp__")]
  rhat_alpha <- rhats[which(is.element(variables, paste0("alpha[", 1:n_mk, "]")))]

  # Trace plot
  fit_ae$draws() %>% mcmc_trace(pars = paste0("alpha[", 1:n_mk, "]")) -> traceplot_ae

  return(list("median_alpha" = median_alpha,
              "rhat_lp" = rhat_lp,
              "rhat_alpha" = rhat_alpha,
              "traceplot_ae" = traceplot_ae))
}
