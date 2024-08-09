#' set_consider_sr_model
#'
#' @description The function to set consideration of each stutter ratio model
#' @param meth_model_sr The modeling method of the stutter in one marker
#' @param reps_one_mk Repeat motifs of one marker
#' @param min_n_motif_lus The minimum number of repeat motifs for considering the LUS model
#' @param min_n_motif_multi The minimum number of repeat motifs for considering the multi-sequence model
set_consider_sr_model <- function(meth_model_sr, reps_one_mk, min_n_motif_lus, min_n_motif_multi){
  consider_sr_model <- rep(FALSE, 3)
  names(consider_sr_model) <- c("Allele", "LUS", "Multi-sequence")
  if(meth_model_sr == "Allele"){
    consider_sr_model["Allele"] <- TRUE
  }else if(meth_model_sr == "LUS"){
    consider_sr_model["LUS"] <- TRUE
  }else if(meth_model_sr == "Multi-sequence"){
    consider_sr_model["Multi-sequence"] <- TRUE
  }else{
    consider_sr_model["Allele"] <- TRUE
    max_motif <- max(sapply(reps_one_mk, length))
    if(max_motif >= min_n_motif_lus){
      consider_sr_model["LUS"] <- TRUE
    }
    if(max_motif >= min_n_motif_multi){
      consider_sr_model["Multi-sequence"] <- TRUE
    }
  }
  return(consider_sr_model)
}

#' calc_average_sr
#'
#' @description The function to calculate the average stutter ratio
#' @param allele A vector of alleles or corrected alleles
#' @param slope The slope parameter
#' @param intercept The intercept parameter
calc_average_sr <- function(allele, slope, intercept){
  return(allele * slope + intercept)
}

#' make_dt_alcor
#'
#' @description The function to make a data.table for corrected alleles
#' @param markers_auto The vector of autosomal markers
#' @param sr_best_model_all_mk The vector of the best SR model for all markers
#' @param sr_best_median_myu_all_mk The vector of the myu values for all markers under the best SR model
#' @param sr_best_median_slope_all_mk The vector of the slope values for all markers under the best SR model
#' @param sr_best_median_intercept_all_mk The vector of the intercept values for all markers under the best SR model
#' @param sr_alcor_all_mk The list of the information on corrected alleles of all markers
make_dt_alcor <- function(markers_auto, sr_best_model_all_mk,
                          sr_best_median_myu_all_mk, sr_best_median_slope_all_mk, sr_best_median_intercept_all_mk,
                          sr_alcor_all_mk){
  vec_mk <- vec_model <- character(0)
  vec_info_al <- vec_info_ca <- vec_average_sr <- numeric(0)

  n_mk <- length(markers_auto)
  for(i in 1:n_mk){
    mk <- markers_auto[i]
    model <- sr_best_model_all_mk[i]
    myu <- sr_best_median_myu_all_mk[i]
    slope <- sr_best_median_slope_all_mk[i]
    intercept <- sr_best_median_intercept_all_mk[i]
    sr_alcor_one_mk <- sr_alcor_all_mk[[i]]
    info_al <- sr_alcor_one_mk$info_al
    info_ca <- sr_alcor_one_mk$info_ca

    n_al <- length(info_al)
    vec_mk <- c(vec_mk, rep(mk, n_al))
    vec_model <- c(vec_model, rep(model, n_al))
    vec_info_al <- c(vec_info_al, info_al)
    vec_info_ca <- c(vec_info_ca, info_ca)

    if(is.element(model, c("Allele", "LUS", "Multi-sequence"))){
      average_sr <- calc_average_sr(info_ca, slope, intercept)
      vec_average_sr <- c(vec_average_sr, average_sr)
    }else{
      vec_average_sr <- c(vec_average_sr, rep(myu, n_al))
    }
  }

  dt_alcor <- data.table(Marker = vec_mk,
                         Model = vec_model,
                         Allele = vec_info_al,
                         Corrected = vec_info_ca,
                         Average_SR = vec_average_sr)
  return(dt_alcor)
}

#' make_dt_alcor_combined
#'
#' @description The function to combine each dt_alcor
#' @param dt_alcor_bsr The data.table of corrected alleles for BSR
#' @param dt_alcor_fsr The data.table of corrected alleles for FSR
#' @param dt_alcor_dsr The data.table of corrected alleles for DSR
make_dt_alcor_combined <- function(dt_alcor_bsr, dt_alcor_fsr, dt_alcor_dsr){
  names(dt_alcor_bsr) <- c("Marker", "Model_BSR", "Allele", "Corrected_BSR", "Average_BSR")
  names(dt_alcor_fsr) <- c("Marker", "Model_FSR", "Allele", "Corrected_FSR", "Average_FSR")
  names(dt_alcor_dsr) <- c("Marker", "Model_DSR", "Allele", "Corrected_DSR", "Average_DSR")
  dt_alcor_combined <- copy(dt_alcor_bsr)
  dt_alcor_combined <- full_join(dt_alcor_combined, dt_alcor_fsr, by = c("Marker", "Allele"))
  dt_alcor_combined <- full_join(dt_alcor_combined, dt_alcor_dsr, by = c("Marker", "Allele"))

  # Change the data type in each column
  col_char <- c("Marker", "Model_BSR", "Allele", "Model_FSR", "Model_DSR")
  dt_alcor_combined[, (col_char) := lapply(.SD, as.character), .SDcols = col_char]
  col_numeric <- c("Corrected_BSR", "Average_BSR", "Corrected_FSR", "Average_FSR", "Corrected_DSR", "Average_DSR")
  dt_alcor_combined[, (col_numeric) := lapply(.SD, as.numeric), .SDcols = col_numeric]

  return(dt_alcor_combined)
}

#' make_dt_sr
#'
#' @description The function to make a data.table for the stutter ratio
#' @param dt_alcor_combined The data.table of the corrected alleles of BSR, FSR, and DSR
#' @param dt_result_m2sr The data.table of the results for M2SR
make_dt_sr <- function(dt_alcor_combined, dt_result_m2sr){
  index_consider_m2sr <- which(dt_result_m2sr$Model != "Not consider")
  markers_m2sr_consider <- dt_result_m2sr$Marker[index_consider_m2sr]
  models_m2sr_consider <- dt_result_m2sr$Model[index_consider_m2sr]
  myus_m2sr_consider <- dt_result_m2sr$Myu[index_consider_m2sr]

  n_data_alcor <- nrow(dt_alcor_combined)
  model_m2sr_alcor <- rep("Not consider", n_data_alcor)
  value_m2sr_alcor <- rep(0, n_data_alcor)
  for(i in 1:length(markers_m2sr_consider)){
    index_mk_m2sr <- which(dt_alcor_combined$Marker == markers_m2sr_consider[i])
    model_m2sr_alcor[index_mk_m2sr] <- models_m2sr_consider[i]
    value_m2sr_alcor[index_mk_m2sr] <- myus_m2sr_consider[i]
  }

  dt_sr <- data.table(Marker = dt_alcor_combined$Marker,
                      Allele = dt_alcor_combined$Allele,
                      Model_BSR = dt_alcor_combined$Model_BSR,
                      Corrected_BSR = dt_alcor_combined$Corrected_BSR,
                      Average_BSR = dt_alcor_combined$Average_BSR,
                      Model_FSR = dt_alcor_combined$Model_FSR,
                      Corrected_FSR = dt_alcor_combined$Corrected_FSR,
                      Average_FSR = dt_alcor_combined$Average_FSR,
                      Model_DSR = dt_alcor_combined$Model_DSR,
                      Corrected_DSR = dt_alcor_combined$Corrected_DSR,
                      Average_DSR = dt_alcor_combined$Average_DSR,
                      Model_M2SR = model_m2sr_alcor,
                      Average_M2SR = value_m2sr_alcor)
  return(dt_sr)
}

#' make_dt_aic
#'
#' @description The function to make a data.table for information on AIC values and parameters of the candidate models
#' @param markers_auto The vector of autosomal markers
#' @param meth_model_sr_auto The vector of modeling methods of all markers
#' @param aics_sr_all_mk The list of AIC values of all markers
#' @param sr_allele_median_slope_all_mk The slope parameters of all markers under the Allele model
#' @param sr_allele_median_intercept_all_mk The intercept parameters of all markers under the Allele model
#' @param sr_allele_median_sigma_all_mk The sigma parameters of all markers under the Allele model
#' @param sr_allele_rhat_slope_all_mk The Rhat values for the slope parameters of all markers under the Allele model
#' @param sr_allele_rhat_intercept_all_mk The Rhat values for the intercept parameters of all markers under the Allele model
#' @param sr_allele_rhat_sigma_all_mk The Rhat values for the sigma parameters of all markers under the Allele model
#' @param sr_allele_rhat_lp_all_mk The Rhat values for the likelihoods of all markers under the Allele model
#' @param sr_lus_median_slope_all_mk The slope parameters of all markers under the LUS model
#' @param sr_lus_median_intercept_all_mk The intercept parameters of all markers under the LUS model
#' @param sr_lus_median_sigma_all_mk The sigma parameters of all markers under the LUS model
#' @param sr_lus_rhat_slope_all_mk The Rhat values for the slope parameters of all markers under the LUS model
#' @param sr_lus_rhat_intercept_all_mk The Rhat values for the intercept parameters of all markers under the LUS model
#' @param sr_lus_rhat_sigma_all_mk The Rhat values for the sigma parameters of all markers under the LUS model
#' @param sr_lus_rhat_lp_all_mk The Rhat values for the likelihoods of all markers under the LUS model
#' @param sr_multi_median_slope_all_mk The slope parameters of all markers under the Multi-sequence model
#' @param sr_multi_median_intercept_all_mk The intercept parameters of all markers under the Multi-sequence model
#' @param sr_multi_median_sigma_all_mk The sigma parameters of all markers under the Multi-sequence model
#' @param sr_multi_min_us_all_mk The MinUS parameters of all markers under the Multi-sequence model
#' @param sr_multi_rhat_slope_all_mk The Rhat values for the slope parameters of all markers under the Multi-sequence model
#' @param sr_multi_rhat_intercept_all_mk The Rhat values for the intercept parameters of all markers under the Multi-sequence model
#' @param sr_multi_rhat_sigma_all_mk The Rhat values for the sigma parameters of all markers under the Multi-sequence model
#' @param sr_multi_rhat_lp_all_mk The Rhat values for the likelihoods of all markers under the Multi-sequence model
make_dt_aic <- function(markers_auto, meth_model_sr_auto, aics_sr_all_mk,
                        sr_allele_median_slope_all_mk, sr_allele_median_intercept_all_mk, sr_allele_median_sigma_all_mk,
                        sr_allele_rhat_slope_all_mk, sr_allele_rhat_intercept_all_mk, sr_allele_rhat_sigma_all_mk, sr_allele_rhat_lp_all_mk,
                        sr_lus_median_slope_all_mk, sr_lus_median_intercept_all_mk, sr_lus_median_sigma_all_mk,
                        sr_lus_rhat_slope_all_mk, sr_lus_rhat_intercept_all_mk, sr_lus_rhat_sigma_all_mk, sr_lus_rhat_lp_all_mk,
                        sr_multi_median_slope_all_mk, sr_multi_median_intercept_all_mk, sr_multi_median_sigma_all_mk, sr_multi_min_us_all_mk,
                        sr_multi_rhat_slope_all_mk, sr_multi_rhat_intercept_all_mk, sr_multi_rhat_sigma_all_mk, sr_multi_rhat_lp_all_mk){
  n_mk <- length(markers_auto)
  sr_markers <- sr_models <- character(0)
  sr_aics <- sr_slopes <- sr_intercepts <- sr_sigmas <- sr_rhat_slopes <- sr_rhat_intercepts <- sr_rhat_sigmas <- sr_rhat_lps <- numeric(0)
  sr_min_us <- integer(0)
  for(i in 1:n_mk){
    if(meth_model_sr_auto[i] == "Best"){
      mk <- markers_auto[i]
      aics <- aics_sr_all_mk[[i]]
      if(!is.null(aics)){
        models_cand <- names(aics)
        sr_aics <- c(sr_aics, aics)

        if(is.element("Allele", models_cand)){
          sr_markers <- c(sr_markers, mk)
          sr_models <- c(sr_models, "Allele")
          sr_slopes <- c(sr_slopes, sr_allele_median_slope_all_mk[i])
          sr_intercepts <- c(sr_intercepts, sr_allele_median_intercept_all_mk[i])
          sr_sigmas <- c(sr_sigmas, sr_allele_median_sigma_all_mk[i])
          sr_min_us <- c(sr_min_us, NA)
          sr_rhat_slopes <- c(sr_rhat_slopes, sr_allele_rhat_slope_all_mk[i])
          sr_rhat_intercepts <- c(sr_rhat_intercepts, sr_allele_rhat_intercept_all_mk[i])
          sr_rhat_sigmas <- c(sr_rhat_sigmas, sr_allele_rhat_sigma_all_mk[i])
          sr_rhat_lps <- c(sr_rhat_lps, sr_allele_rhat_lp_all_mk[i])
        }

        if(is.element("LUS", models_cand)){
          sr_markers <- c(sr_markers, mk)
          sr_models <- c(sr_models, "LUS")
          sr_slopes <- c(sr_slopes, sr_lus_median_slope_all_mk[i])
          sr_intercepts <- c(sr_intercepts, sr_lus_median_intercept_all_mk[i])
          sr_sigmas <- c(sr_sigmas, sr_lus_median_sigma_all_mk[i])
          sr_min_us <- c(sr_min_us, NA)
          sr_rhat_slopes <- c(sr_rhat_slopes, sr_lus_rhat_slope_all_mk[i])
          sr_rhat_intercepts <- c(sr_rhat_intercepts, sr_lus_rhat_intercept_all_mk[i])
          sr_rhat_sigmas <- c(sr_rhat_sigmas, sr_lus_rhat_sigma_all_mk[i])
          sr_rhat_lps <- c(sr_rhat_lps, sr_lus_rhat_lp_all_mk[i])
        }

        if(is.element("Multi-sequence", models_cand)){
          sr_markers <- c(sr_markers, mk)
          sr_models <- c(sr_models, "Multi-sequence")
          sr_slopes <- c(sr_slopes, sr_multi_median_slope_all_mk[i])
          sr_intercepts <- c(sr_intercepts, sr_multi_median_intercept_all_mk[i])
          sr_sigmas <- c(sr_sigmas, sr_multi_median_sigma_all_mk[i])
          sr_min_us <- c(sr_min_us, sr_multi_min_us_all_mk[i])
          sr_rhat_slopes <- c(sr_rhat_slopes, sr_multi_rhat_slope_all_mk[i])
          sr_rhat_intercepts <- c(sr_rhat_intercepts, sr_multi_rhat_intercept_all_mk[i])
          sr_rhat_sigmas <- c(sr_rhat_sigmas, sr_multi_rhat_sigma_all_mk[i])
          sr_rhat_lps <- c(sr_rhat_lps, sr_multi_rhat_lp_all_mk[i])
        }
      }
    }
  }

  dt_aic <- data.table(Marker = sr_markers,
                       Model = sr_models,
                       AIC = sr_aics,
                       Slope = sr_slopes,
                       Intercept = sr_intercepts,
                       Sigma = sr_sigmas,
                       MinUS = sr_min_us,
                       Rhat_Slope = sr_rhat_slopes,
                       Rhat_Intercept = sr_rhat_intercepts,
                       Rhat_Sigma = sr_rhat_sigmas,
                       Rhat_Lp = sr_rhat_lps
                       )
  return(dt_aic)
}

#' analyze_bfdstr_main
#'
#' @description The function to analyze BSR, FSR, or DSR data
#' @param sr_type "BSR", "FSR", or "DSR"
#' @param markers_auto The vector of autosomal markers
#' @param meth_model_sr_auto The vector of modeling methods of autosomal markers
#' @param model_sr_gen The stan model for locus-specific models with the slope parameter
#' @param model_sr_uni The stan model for the uniform model
#' @param sr_data The stutter ratios of autosomal markers
#' @param parent_al_data The parent alleles for each stutter ratios of autosomal markers
#' @param total_al_prod_data The total allelic products of autosomal markers
#' @param lus_auto The information on the LUS values of autosomal markers
#' @param reps_auto The information on the repeat motifs of all alleles
#' @param seq_al_auto The information on the alleles for each sequence
#' @param seq_count_auto The information on the number of counts for each sequence
#' @param n_warmup The number of simulations for warm-up
#' @param n_sampling The number of simulations for sampling
#' @param n_chains The number of chains
#' @param n_parallel_chains The number of parallel chains
#' @param min_n_motif_lus The minimum number of repeat motifs for considering the LUS model
#' @param min_n_motif_multi The minimum number of repeat motifs for considering the multi-sequence model
#' @param max_n_motif_multi The maximum number of repeat motifs for considering the multi-sequence model
analyze_bfdstr_main <- function(sr_type, markers_auto, meth_model_sr_auto, model_sr_gen, model_sr_uni,
                                sr_data, parent_al_data, total_al_prod_data,
                                lus_auto, reps_auto, seq_al_auto, seq_count_auto,
                                n_warmup, n_sampling, n_chains, n_parallel_chains, min_n_motif_lus, min_n_motif_multi, max_n_motif_multi){
  # Define objects
  n_mk <- length(markers_auto)
  sr_data_all_mk <- vector(mode = "list", length = n_mk)
  names(sr_data_all_mk) <- markers_auto
  sr_alcor_all_mk <- vector(mode = "list", length = n_mk)
  names(sr_alcor_all_mk) <- markers_auto
  aics_sr_all_mk <- vector(mode = "list", length = n_mk)
  names(aics_sr_all_mk) <- markers_auto
  sr_allele_median_slope_all_mk <- sr_allele_median_intercept_all_mk <- sr_allele_median_sigma_all_mk <- rep(NA, n_mk)
  sr_allele_rhat_slope_all_mk <- sr_allele_rhat_intercept_all_mk <- sr_allele_rhat_sigma_all_mk <- sr_allele_rhat_lp_all_mk <- rep(NA, n_mk)
  sr_allele_trace_all_mk <- vector(mode = "list", length = n_mk)
  sr_lus_median_slope_all_mk <- sr_lus_median_intercept_all_mk <- sr_lus_median_sigma_all_mk <- rep(NA, n_mk)
  sr_lus_rhat_slope_all_mk <- sr_lus_rhat_intercept_all_mk <- sr_lus_rhat_sigma_all_mk <- sr_lus_rhat_lp_all_mk <- rep(NA, n_mk)
  sr_lus_trace_all_mk <- vector(mode = "list", length = n_mk)
  sr_multi_median_slope_all_mk <- sr_multi_median_intercept_all_mk <- sr_multi_median_sigma_all_mk <- sr_multi_min_us_all_mk <- rep(NA, n_mk)
  sr_multi_rhat_slope_all_mk <- sr_multi_rhat_intercept_all_mk <- sr_multi_rhat_sigma_all_mk <- sr_multi_rhat_lp_all_mk <- rep(NA, n_mk)
  sr_multi_trace_all_mk <- vector(mode = "list", length = n_mk)
  sr_best_model_all_mk <- rep(NA, n_mk)
  sr_best_median_myu_all_mk <- sr_best_median_slope_all_mk <- sr_best_median_intercept_all_mk <- sr_best_median_sigma_all_mk <- sr_best_min_us_all_mk <- rep(NA, n_mk)
  sr_best_rhat_myu_all_mk <- sr_best_rhat_slope_all_mk <- sr_best_rhat_intercept_all_mk <- sr_best_rhat_sigma_all_mk <- sr_best_rhat_lp_all_mk <- rep(NA, n_mk)
  sr_best_trace_all_mk <- vector(mode = "list", length = n_mk)

  # locus-specific model
  index_mk_specific_sr <- which(is.element(meth_model_sr_auto, c("Best", "Allele", "LUS", "Multi-sequence")))
  n_mk_specific_sr <- length(index_mk_specific_sr)
  if(n_mk_specific_sr > 0){
    for(i in 1:n_mk_specific_sr){
      id_mk <- index_mk_specific_sr[i]

      # Extract data of one marker
      tmp <- extract_sr_use(sr_data[, id_mk], parent_al_data[, id_mk], total_al_prod_data[, id_mk], lus_auto[[id_mk]])
      sr_one_mk <- tmp[[1]]
      sr_al_one_mk <- tmp[[2]]
      sr_h_one_mk <- tmp[[3]]
      sr_lus_one_mk <- tmp[[4]]
      reps_one_mk <- reps_auto[[id_mk]]
      seq_al_one_mk <- seq_al_auto[[id_mk]]
      seq_count_one_mk <- seq_count_auto[[id_mk]]

      # Set consideration of each model
      model_sr <- meth_model_sr_auto[id_mk]
      consider_sr <- set_consider_sr_model(model_sr, reps_one_mk, min_n_motif_lus, min_n_motif_multi)

      # Define an object to assign log-likelihoods of each model
      log_likes <- rep(NA, 3)
      names(log_likes) <- c("Allele", "LUS", "Multi-sequence")

      # Allele model
      if(consider_sr["Allele"]){
        cat(paste0("MCMC for the ", sr_type, " of the ", markers_auto[id_mk], ": Allele model\n"))
        data_list_sr_allele <- list(N = length(sr_one_mk), SR = sr_one_mk, XVal = sr_al_one_mk, Height = sr_h_one_mk)
        tmp <- run_mcmc_sr_gen(model_sr_gen, data_list_sr_allele, n_warmup, n_sampling, n_chains, n_parallel_chains)
        log_likes[1] <- tmp$median_lp
        sr_allele_median_slope_all_mk[id_mk] <- tmp$median_slope
        sr_allele_median_intercept_all_mk[id_mk] <- tmp$median_intercept
        sr_allele_median_sigma_all_mk[id_mk] <- tmp$median_sigma
        sr_allele_rhat_slope_all_mk[id_mk] <- tmp$rhat_slope
        sr_allele_rhat_intercept_all_mk[id_mk] <- tmp$rhat_intercept
        sr_allele_rhat_sigma_all_mk[id_mk] <- tmp$rhat_sigma
        sr_allele_rhat_lp_all_mk[id_mk] <- tmp$rhat_lp
        sr_allele_trace_all_mk[[id_mk]] <- tmp$traceplot
      }

      # LUS model
      if(consider_sr["LUS"]){
        cat(paste0("MCMC for the ", sr_type, " of the ", markers_auto[id_mk], ": LUS model\n"))
        data_list_sr_lus <- list(N = length(sr_one_mk), SR = sr_one_mk, XVal = sr_lus_one_mk, Height = sr_h_one_mk)
        tmp <- run_mcmc_sr_gen(model_sr_gen, data_list_sr_lus, n_warmup, n_sampling, n_chains, n_parallel_chains)
        log_likes[2] <- tmp$median_lp
        sr_lus_median_slope_all_mk[id_mk] <- tmp$median_slope
        sr_lus_median_intercept_all_mk[id_mk] <- tmp$median_intercept
        sr_lus_median_sigma_all_mk[id_mk] <- tmp$median_sigma
        sr_lus_rhat_slope_all_mk[id_mk] <- tmp$rhat_slope
        sr_lus_rhat_intercept_all_mk[id_mk] <- tmp$rhat_intercept
        sr_lus_rhat_sigma_all_mk[id_mk] <- tmp$rhat_sigma
        sr_lus_rhat_lp_all_mk[id_mk] <- tmp$rhat_lp
        sr_lus_trace_all_mk[[id_mk]] <- tmp$traceplot
      }

      # Multi-sequence model
      if(consider_sr["Multi-sequence"]){
        cat(paste0("MCMC for the ", sr_type, " of the ", markers_auto[id_mk], ": Multi-sequence model\n"))
        cand_min_us <- sort(unique(unlist(reps_one_mk)))
        cand_min_us <- cand_min_us[cand_min_us <= max_n_motif_multi]
        n_cand_min_us <- length(cand_min_us)
        loglike_multi <- rep(NA, n_cand_min_us)
        sr_ca_one_mk_cand <- result_sr_multi_cand <- list()
        for(j in 1:n_cand_min_us){
          sr_ca_one_mk_cand[[j]] <- sr_ca_one_mk <- make_sr_ca_one_mk(cand_min_us[j], reps_one_mk, seq_al_one_mk, seq_count_one_mk, sr_al_one_mk)
          data_list_sr_multi <- list(N = length(sr_one_mk), SR = sr_one_mk, XVal = sr_ca_one_mk, Height = sr_h_one_mk)
          result_sr_multi_cand[[j]] <- run_mcmc_sr_gen(model_sr_gen, data_list_sr_multi, n_warmup, n_sampling, n_chains, n_parallel_chains)
          loglike_multi[j] <- result_sr_multi_cand[[j]]$median_lp
        }
        index_max_like_multi <- which.max(loglike_multi)
        sr_ca_one_mk <- sr_ca_one_mk_cand[[index_max_like_multi]]
        tmp <- result_sr_multi_cand[[index_max_like_multi]]
        log_likes[3] <- tmp$median_lp
        sr_multi_median_slope_all_mk[id_mk] <- tmp$median_slope
        sr_multi_median_intercept_all_mk[id_mk] <- tmp$median_intercept
        sr_multi_median_sigma_all_mk[id_mk] <- tmp$median_sigma
        sr_multi_min_us_all_mk[id_mk] <- cand_min_us[index_max_like_multi]
        sr_multi_rhat_slope_all_mk[id_mk] <- tmp$rhat_slope
        sr_multi_rhat_intercept_all_mk[id_mk] <- tmp$rhat_intercept
        sr_multi_rhat_sigma_all_mk[id_mk] <- tmp$rhat_sigma
        sr_multi_rhat_lp_all_mk[id_mk] <- tmp$rhat_lp
        sr_multi_trace_all_mk[[id_mk]] <- tmp$traceplot
      }

      # Best model
      index_cand <- which(!is.na(log_likes))
      log_likes <- log_likes[index_cand]
      n_pars <- c(3, 3, 4)
      n_pars <- n_pars[index_cand]
      aics <- mapply(calc_aic, log_likes, n_pars)
      names(aics) <- names(log_likes)
      index_best_model <- which.min(aics)
      if(names(aics)[index_best_model] == "Allele"){
        sr_best_model_all_mk[id_mk] <- "Allele"
        sr_best_median_slope_all_mk[id_mk] <- sr_allele_median_slope_all_mk[id_mk]
        sr_best_median_intercept_all_mk[id_mk] <- sr_allele_median_intercept_all_mk[id_mk]
        sr_best_median_sigma_all_mk[id_mk] <- sr_allele_median_sigma_all_mk[id_mk]
        sr_best_rhat_slope_all_mk[id_mk] <- sr_allele_rhat_slope_all_mk[id_mk]
        sr_best_rhat_intercept_all_mk[id_mk] <- sr_allele_rhat_intercept_all_mk[id_mk]
        sr_best_rhat_sigma_all_mk[id_mk] <- sr_allele_rhat_sigma_all_mk[id_mk]
        sr_best_rhat_lp_all_mk[id_mk] <- sr_allele_rhat_lp_all_mk[id_mk]
        sr_best_trace_all_mk[[id_mk]] <- sr_allele_trace_all_mk[[id_mk]]
      }else if(names(aics)[index_best_model] == "LUS"){
        sr_best_model_all_mk[id_mk] <- "LUS"
        sr_best_median_slope_all_mk[id_mk] <- sr_lus_median_slope_all_mk[id_mk]
        sr_best_median_intercept_all_mk[id_mk] <- sr_lus_median_intercept_all_mk[id_mk]
        sr_best_median_sigma_all_mk[id_mk] <- sr_lus_median_sigma_all_mk[id_mk]
        sr_best_rhat_slope_all_mk[id_mk] <- sr_lus_rhat_slope_all_mk[id_mk]
        sr_best_rhat_intercept_all_mk[id_mk] <- sr_lus_rhat_intercept_all_mk[id_mk]
        sr_best_rhat_sigma_all_mk[id_mk] <- sr_lus_rhat_sigma_all_mk[id_mk]
        sr_best_rhat_lp_all_mk[id_mk] <- sr_lus_rhat_lp_all_mk[id_mk]
        sr_best_trace_all_mk[[id_mk]] <- sr_lus_trace_all_mk[[id_mk]]
      }else{
        sr_best_model_all_mk[id_mk] <- "Multi-sequence"
        sr_best_median_slope_all_mk[id_mk] <- sr_multi_median_slope_all_mk[id_mk]
        sr_best_median_intercept_all_mk[id_mk] <- sr_multi_median_intercept_all_mk[id_mk]
        sr_best_median_sigma_all_mk[id_mk] <- sr_multi_median_sigma_all_mk[id_mk]
        sr_best_min_us_all_mk[id_mk] <- sr_multi_min_us_all_mk[id_mk]
        sr_best_rhat_slope_all_mk[id_mk] <- sr_multi_rhat_slope_all_mk[id_mk]
        sr_best_rhat_intercept_all_mk[id_mk] <- sr_multi_rhat_intercept_all_mk[id_mk]
        sr_best_rhat_sigma_all_mk[id_mk] <- sr_multi_rhat_sigma_all_mk[id_mk]
        sr_best_rhat_lp_all_mk[id_mk] <- sr_multi_rhat_lp_all_mk[id_mk]
        sr_best_trace_all_mk[[id_mk]] <- sr_multi_trace_all_mk[[id_mk]]
      }

      # Assign data
      sr_data_one_mk <- vector(mode = "list", length = 3)
      names(sr_data_one_mk) <- c("sr", "sr_al", "sr_ca")
      sr_data_one_mk[[1]] <- sr_one_mk
      if(consider_sr["Allele"]){
        sr_data_one_mk[[2]] <- sr_al_one_mk
      }
      if(consider_sr["LUS"]){
        sr_data_one_mk[[3]] <- sr_lus_one_mk
      }
      if(consider_sr["Multi-sequence"]){
        sr_data_one_mk[[3]] <- sr_ca_one_mk
      }
      sr_data_all_mk[[id_mk]] <- sr_data_one_mk
      aics_sr_all_mk[[id_mk]] <- aics

      # Allele repeat correction
      info_al_one_mk <- as.numeric(names(lus_auto[[id_mk]]))
      if(sr_best_model_all_mk[id_mk] == "Allele"){
        info_ca_one_mk <- info_al_one_mk
      }else if(sr_best_model_all_mk[id_mk] == "LUS"){
        info_ca_one_mk <- lus_auto[[id_mk]]
      }else if(sr_best_model_all_mk[id_mk] == "Multi-sequence"){
        info_ca_one_mk <- make_sr_ca_one_mk(cand_min_us[index_max_like_multi], reps_one_mk, seq_al_one_mk, seq_count_one_mk, info_al_one_mk)
      }else{
        info_ca_one_mk <- rep(NA, length(info_al_one_mk))
      }
      sr_alcor_all_mk[[id_mk]] <- list(info_al = info_al_one_mk,
                                       info_ca = info_ca_one_mk)
    }
  }

  # Multiple loci together
  index_mk_together_sr <- which(meth_model_sr_auto == "Multiple loci together")
  n_mk_together_sr <- length(index_mk_together_sr)
  sr_together <- sr_h_together <- id_mk_together <- numeric(0)
  id_mk_count <- 0
  if(n_mk_together_sr > 0){
    cat(paste0("MCMC for the ", sr_type, ": Multiple loci together\n"))

    # Prepare data of the target loci
    for(i in 1:n_mk_together_sr){
      # Extract data of one marker
      id_mk <- index_mk_together_sr[i]
      tmp <- extract_sr_use(sr_data[, id_mk], parent_al_data[, id_mk], total_al_prod_data[, id_mk], lus_auto[[id_mk]])
      sr_one_mk <- tmp[[1]]
      sr_al_one_mk <- tmp[[2]]
      sr_h_one_mk <- tmp[[3]]

      # Combine data
      sr_together <- c(sr_together, sr_one_mk)
      sr_h_together <- c(sr_h_together, sr_h_one_mk)

      # Make id_mk_together
      id_mk_count <- id_mk_count + 1
      id_mk_together <- c(id_mk_together, rep(id_mk_count, length(sr_one_mk)))

      # Assign data
      sr_data_one_mk <- vector(mode = "list", length = 3)
      names(sr_data_one_mk) <- c("sr", "sr_al", "sr_ca")
      sr_data_one_mk[[1]] <- sr_one_mk
      sr_data_one_mk[[2]] <- sr_al_one_mk
      sr_data_all_mk[[id_mk]] <- sr_data_one_mk

      # Allele repeat correction
      info_al_one_mk <- as.numeric(names(lus_auto[[id_mk]]))
      info_ca_one_mk <- rep(NA, length(info_al_one_mk))
      sr_alcor_all_mk[[id_mk]] <- list(info_al = info_al_one_mk,
                                        info_ca = info_ca_one_mk)
    }
    data_list_sr_together <- list(N = length(sr_together), SR = sr_together, Height = sr_h_together)
    tmp <- run_mcmc_sr_uni(model_sr_uni, data_list_sr_together, n_warmup, n_sampling, n_chains, n_parallel_chains)
    sr_best_model_all_mk[index_mk_together_sr] <- "Multiple loci together"
    sr_best_median_myu_all_mk[index_mk_together_sr] <- tmp$median_myu
    sr_best_median_sigma_all_mk[index_mk_together_sr] <- tmp$median_sigma
    sr_best_rhat_myu_all_mk[index_mk_together_sr] <- tmp$rhat_myu
    sr_best_rhat_sigma_all_mk[index_mk_together_sr] <- tmp$rhat_sigma
    sr_best_rhat_lp_all_mk[index_mk_together_sr] <- tmp$rhat_lp
    for(i in 1:n_mk_together_sr){
      sr_best_trace_all_mk[[index_mk_together_sr[i]]] <- tmp$traceplot
    }
  }

  # Not consider
  index_mk_not_consider_sr <- which(meth_model_sr_auto == "Not consider")
  if(length(index_mk_not_consider_sr) > 0){
    sr_best_model_all_mk[index_mk_not_consider_sr] <- "Not consider"
    for(i in 1:length(index_mk_not_consider_sr)){
      id_mk <- index_mk_not_consider_sr[i]
      info_al_one_mk <- as.numeric(names(lus_auto[[id_mk]]))
      info_ca_one_mk <- rep(NA, length(info_al_one_mk))
      sr_alcor_all_mk[[id_mk]] <- list(info_al = info_al_one_mk,
                                       info_ca = info_ca_one_mk)
    }
  }

  # Make a data.table for results
  dt_result <- data.table(Marker = markers_auto,
                          Model = sr_best_model_all_mk,
                          Myu = sr_best_median_myu_all_mk,
                          Slope = sr_best_median_slope_all_mk,
                          Intercept = sr_best_median_intercept_all_mk,
                          Sigma = sr_best_median_sigma_all_mk,
                          MinUS = sr_best_min_us_all_mk,
                          Rhat_Myu = sr_best_rhat_myu_all_mk,
                          Rhat_Slope = sr_best_rhat_slope_all_mk,
                          Rhat_Intercept = sr_best_rhat_intercept_all_mk,
                          Rhat_Sigma = sr_best_rhat_sigma_all_mk,
                          Rhat_Lp = sr_best_rhat_lp_all_mk)

  ## sr: Make a data.table for corrected alleles
  dt_alcor <- make_dt_alcor(markers_auto, sr_best_model_all_mk,
                            sr_best_median_myu_all_mk, sr_best_median_slope_all_mk, sr_best_median_intercept_all_mk,
                            sr_alcor_all_mk)

  # Make a data.table for AIC
  dt_aic <- make_dt_aic(markers_auto, meth_model_sr_auto, aics_sr_all_mk,
                        sr_allele_median_slope_all_mk, sr_allele_median_intercept_all_mk, sr_allele_median_sigma_all_mk,
                        sr_allele_rhat_slope_all_mk, sr_allele_rhat_intercept_all_mk, sr_allele_rhat_sigma_all_mk, sr_allele_rhat_lp_all_mk,
                        sr_lus_median_slope_all_mk, sr_lus_median_intercept_all_mk, sr_lus_median_sigma_all_mk,
                        sr_lus_rhat_slope_all_mk, sr_lus_rhat_intercept_all_mk, sr_lus_rhat_sigma_all_mk, sr_lus_rhat_lp_all_mk,
                        sr_multi_median_slope_all_mk, sr_multi_median_intercept_all_mk, sr_multi_median_sigma_all_mk, sr_multi_min_us_all_mk,
                        sr_multi_rhat_slope_all_mk, sr_multi_rhat_intercept_all_mk, sr_multi_rhat_sigma_all_mk, sr_multi_rhat_lp_all_mk)

  return(list(dt_result = dt_result,
              dt_alcor = dt_alcor,
              dt_aic = dt_aic,
              traceplot = sr_best_trace_all_mk,
              plotdata = sr_data_all_mk))
}

#' analyze_m2sr_main
#'
#' @description The function to analyze M2SR data
#' @param markers_auto The vector of autosomal markers
#' @param meth_model_sr_auto The vector of modeling methods of autosomal markers
#' @param model_sr_uni The stan model for the uniform model
#' @param sr_data The minus 2-nt stutter ratios of autosomal markers
#' @param parent_al_data The parent alleles for each stutter ratios of autosomal markers
#' @param total_al_prod_data The total allelic products of all autosomal markers
#' @param lus_auto The information on the LUS values of autosomal markers
#' @param n_warmup The number of simulations for warm-up
#' @param n_sampling The number of simulations for sampling
#' @param n_chains The number of chains
#' @param n_parallel_chains The number of parallel chains
analyze_m2sr_main <- function(markers_auto, meth_model_sr_auto, model_sr_uni,
                              sr_data, parent_al_data, total_al_prod_data,
                              lus_auto,
                              n_warmup, n_sampling, n_chains, n_parallel_chains){
  # Define objects
  n_mk <- length(markers_auto)
  sr_best_model_all_mk <- rep(NA, n_mk)
  sr_best_median_myu_all_mk <- sr_best_median_sigma_all_mk <- rep(NA, n_mk)
  sr_best_rhat_myu_all_mk <- sr_best_rhat_sigma_all_mk <- sr_best_rhat_lp_all_mk <- rep(NA, n_mk)
  sr_best_trace_all_mk <- vector(mode = "list", length = n_mk)

  # Locus-specific
  sr_data_all_mk <- vector(mode = "list", length = n_mk)
  names(sr_data_all_mk) <- markers_auto
  index_mk_specific_sr <- which(meth_model_sr_auto == "Uniform")
  n_mk_specific_sr <- length(index_mk_specific_sr)
  if(n_mk_specific_sr > 0){
    for(i in 1:n_mk_specific_sr){
      id_mk <- index_mk_specific_sr[i]
      cat(paste0("MCMC for the M2SR of the ", markers_auto[id_mk], ": Uniform model\n"))

      # Extract data of one marker
      tmp <- extract_sr_use(sr_data[, id_mk], parent_al_data[, id_mk], total_al_prod_data[, id_mk], lus_auto[[id_mk]])
      sr_one_mk <- tmp[[1]]
      sr_al_one_mk <- tmp[[2]]
      sr_h_one_mk <- tmp[[3]]

      data_list_sr_uniform <- list(N = length(sr_one_mk), SR = sr_one_mk, Height = sr_h_one_mk)
      tmp <- run_mcmc_sr_uni(model_sr_uni, data_list_sr_uniform, n_warmup, n_sampling, n_chains, n_parallel_chains)
      sr_best_model_all_mk[id_mk] <- "Uniform"
      sr_best_median_myu_all_mk[id_mk] <- tmp$median_myu
      sr_best_median_sigma_all_mk[id_mk] <- tmp$median_sigma
      sr_best_rhat_myu_all_mk[id_mk] <- tmp$rhat_myu
      sr_best_rhat_sigma_all_mk[id_mk] <- tmp$rhat_sigma
      sr_best_rhat_lp_all_mk[id_mk] <- tmp$rhat_lp
      sr_best_trace_all_mk[[id_mk]] <- tmp$traceplot

      # Assign data
      sr_data_one_mk <- vector(mode = "list", length = 3)
      names(sr_data_one_mk) <- c("sr", "sr_al", "sr_ca")
      sr_data_one_mk[[1]] <- sr_one_mk
      sr_data_one_mk[[2]] <- sr_al_one_mk
      sr_data_all_mk[[id_mk]] <- sr_data_one_mk
    }
  }

  # Multiple loci together
  index_mk_together_sr <- which(meth_model_sr_auto == "Multiple loci together")
  n_mk_together_sr <- length(index_mk_together_sr)
  sr_together <- sr_h_together <- id_mk_together <- numeric(0)
  id_mk_count <- 0
  if(n_mk_together_sr > 0){
    cat("MCMC for the M2SR: Multiple loci together\n")

    # Prepare data of the target loci
    for(i in 1:n_mk_together_sr){
      # Extract data of one marker
      id_mk <- index_mk_together_sr[i]
      tmp <- extract_sr_use(sr_data[, id_mk], parent_al_data[, id_mk], total_al_prod_data[, id_mk], lus_auto[[id_mk]])
      sr_one_mk <- tmp[[1]]
      sr_al_one_mk <- tmp[[2]]
      sr_h_one_mk <- tmp[[3]]

      # Combine data
      sr_together <- c(sr_together, sr_one_mk)
      sr_h_together <- c(sr_h_together, sr_h_one_mk)

      # Make id_mk_together
      id_mk_count <- id_mk_count + 1
      id_mk_together <- c(id_mk_together, rep(id_mk_count, length(sr_one_mk)))

      # Assign data
      sr_data_one_mk <- vector(mode = "list", length = 3)
      names(sr_data_one_mk) <- c("sr", "sr_al", "sr_ca")
      sr_data_one_mk[[1]] <- sr_one_mk
      sr_data_one_mk[[2]] <- sr_al_one_mk
      sr_data_all_mk[[id_mk]] <- sr_data_one_mk
    }
    data_list_sr_together <- list(N = length(sr_together), SR = sr_together, Height = sr_h_together)
    tmp <- run_mcmc_sr_uni(model_sr_uni, data_list_sr_together, n_warmup, n_sampling, n_chains, n_parallel_chains)
    sr_best_model_all_mk[index_mk_together_sr] <- "Multiple loci together"
    sr_best_median_myu_all_mk[index_mk_together_sr] <- tmp$median_myu
    sr_best_median_sigma_all_mk[index_mk_together_sr] <- tmp$median_sigma
    sr_best_rhat_myu_all_mk[index_mk_together_sr] <- tmp$rhat_myu
    sr_best_rhat_sigma_all_mk[index_mk_together_sr] <- tmp$rhat_sigma
    sr_best_rhat_lp_all_mk[index_mk_together_sr] <- tmp$rhat_lp
    for(i in 1:n_mk_together_sr){
      sr_best_trace_all_mk[[index_mk_together_sr[i]]] <- tmp$traceplot
    }
  }

  # Not consider
  index_mk_not_consider_sr <- which(meth_model_sr_auto == "Not consider")
  if(length(index_mk_not_consider_sr) > 0){
    sr_best_model_all_mk[index_mk_not_consider_sr] <- "Not consider"
  }

  # Make a data.table for results
  dt_result <- data.table(Marker = markers_auto,
                          Model = sr_best_model_all_mk,
                          Myu = sr_best_median_myu_all_mk,
                          Sigma = sr_best_median_sigma_all_mk,
                          Rhat_Myu = sr_best_rhat_myu_all_mk,
                          Rhat_Sigma = sr_best_rhat_sigma_all_mk,
                          Rhat_Lp = sr_best_rhat_lp_all_mk)

  return(list(dt_result = dt_result,
              traceplot = sr_best_trace_all_mk,
              plotdata = sr_data_all_mk))
}
