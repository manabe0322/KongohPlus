calc_rphs <- function(gt, peaks, ae_actual, sr_consider, bsr_allele, bsr_average, fsr_allele, fsr_average, dsr_allele, dsr_average, m2sr_myu){
  # Change data type
  gt <- as.numeric(gt)
  peaks <- as.numeric(peaks)

  al <- unique(gt)
  n_al <- length(al)

  if(sr_consider["BS"]){
    bs <- calc_bs(al)
    bsr <- bsr_average[match(al, bsr_allele)]
  }else{
    bs <- bsr <- rep(0, n_al)
  }

  if(sr_consider["FS"]){
    fs <- calc_fs(al)
    fsr <- fsr_average[match(al, fsr_allele)]
  }else{
    fs <- fsr <- rep(0, n_al)
  }

  if(sr_consider["DS"]){
    ds <- calc_ds(al)
    dsr <- dsr_average[match(al, dsr_allele)]
  }else{
    ds <- dsr <- rep(0, n_al)
  }

  if(sr_consider["M2S"]){
    m2s <- calc_m2s(al)
    m2sr <- rep(m2sr_myu, n_al)
  }else{
    m2s <- m2sr <- rep(0, n_al)
  }

  # Base relative peak height
  if(n_al == 1){
    rph_base <- ae_actual
  }else{
    rph_base <- 0.5 * ae_actual
  }

  # Indices of each allele and stutter
  index_al <- match(al, peaks)
  index_bs <- match(bs, peaks)
  index_fs <- match(fs, peaks)
  index_ds <- match(ds, peaks)
  index_m2s <- match(m2s, peaks)

  # Calculate relative peak heights
  n_peak <- length(peaks)
  rphs <- rep(0, n_peak)
  for(i in 1:n_al){
    if(!is.na(index_al[i])){
      rphs[index_al[i]] <- rphs[index_al[i]] + rph_base / (1 + bsr[i] + fsr[i] + dsr[i] + m2sr[i])
    }
    if(!is.na(index_bs[i])){
      rphs[index_bs[i]] <- rphs[index_bs[i]] + rph_base * bsr[i] / (1 + bsr[i] + fsr[i] + dsr[i] + m2sr[i])
    }
    if(!is.na(index_fs[i])){
      rphs[index_fs[i]] <- rphs[index_fs[i]] + rph_base * fsr[i] / (1 + bsr[i] + fsr[i] + dsr[i] + m2sr[i])
    }
    if(!is.na(index_ds[i])){
      rphs[index_ds[i]] <- rphs[index_ds[i]] + rph_base * dsr[i] / (1 + bsr[i] + fsr[i] + dsr[i] + m2sr[i])
    }
    if(!is.na(index_m2s[i])){
      rphs[index_m2s[i]] <- rphs[index_m2s[i]] + rph_base * m2sr[i] / (1 + bsr[i] + fsr[i] + dsr[i] + m2sr[i])
    }
  }
  return(rphs)
}

make_data_gamma <- function(markers_auto, gt_use, peak_use, pae_data_use, tph_data_use, sr_consider_auto, dt_alcor_bsr, dt_alcor_fsr, dt_alcor_dsr, dt_result_m2sr){
  n_gt <- length(gt_use)
  n_mk <- length(markers_auto)

  setindex(dt_alcor_bsr, Marker)
  setindex(dt_alcor_fsr, Marker)
  setindex(dt_alcor_dsr, Marker)
  setindex(dt_result_m2sr, Marker)

  rphs_gamma <- tphs_gamma <- vector(mode = "list", length = n_gt)

  for(i in 1:n_gt){
    gt_s1 <- gt_use[[i]]
    peak_s1 <- peak_use[[i]]
    pae_s1 <- pae_data_use[, i]
    tph_s1 <- tph_data_use[i]

    rphs_gamma_s1 <- tphs_gamma_s1 <- vector(mode = "list", length = n_mk)
    names(rphs_gamma_s1) <- names(tphs_gamma_s1) <- markers_auto

    for(j in 1:n_mk){
      peaks <- peak_s1[[j]]
      ae_actual <- n_mk * pae_s1[j]
      mk <- markers_auto[j]
      dt_alcor_bsr_mk <- dt_alcor_bsr[Marker == mk]
      bsr_allele <- dt_alcor_bsr_mk[, Allele]
      bsr_average <- dt_alcor_bsr_mk[, Average_SR]
      dt_alcor_fsr_mk <- dt_alcor_fsr[Marker == mk]
      fsr_allele <- dt_alcor_fsr_mk[, Allele]
      fsr_average <- dt_alcor_fsr_mk[, Average_SR]
      dt_alcor_dsr_mk <- dt_alcor_dsr[Marker == mk]
      dsr_allele <- dt_alcor_dsr_mk[, Allele]
      dsr_average <- dt_alcor_dsr_mk[, Average_SR]
      dt_result_m2sr_mk <- dt_result_m2sr[Marker == mk]
      m2sr_myu <- as.numeric(dt_result_m2sr_mk[, Myu])

      rphs_gamma_s1[[j]] <- calc_rphs(gt_s1[j, ], peaks, ae_actual, sr_consider_auto[j, ], bsr_allele, bsr_average, fsr_allele, fsr_average, dsr_allele, dsr_average, m2sr_myu)
      tphs_gamma_s1[[j]] <- rep(tph_s1, length(peaks))
    }

    rphs_gamma[[i]] <- rphs_gamma_s1
    tphs_gamma[[i]] <- tphs_gamma_s1
  }

  return(list("rphs_gamma" = rphs_gamma,
              "tphs_gamma" = tphs_gamma))
}

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
  median_cv <- fit_gamma_summary$median[which(fit_gamma_summary$variable == "cv")]

  # rhat of the parameter 'beta' and the likelihood
  rhat_cv <- fit_gamma_summary$rhat[which(fit_gamma_summary$variable == "cv")]
  rhat_lp <- fit_gamma_summary$rhat[which(fit_gamma_summary$variable == "lp__")]

  # Trace plot
  fit_gamma$draws() %>% mcmc_trace(pars = c("cv")) -> traceplot

  return(list("median_cv" = median_cv,
              "rhat_cv" = rhat_cv,
              "rhat_lp" = rhat_lp,
              "traceplot" = traceplot))
}
