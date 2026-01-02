#' calc_ave_het_ph_all_data
#'
#' @description The function to calculate the average heterozygote allele peak heights of all marker
#' @param height_list The list of peak heights
calc_ave_het_ph_all_data <- function(height_list){
  n_data <- length(height_list)
  ave_het_ph_all_data <- rep(0, n_data)
  for(i in 1:n_data){
    height_list_one_data <- height_list[[i]]
    ave_het_ph_per_mk <- sapply(height_list_one_data, sum)
    ave_het_ph_all_data[i] <- mean(ave_het_ph_per_mk) / 2
  }
  return(ave_het_ph_all_data)
}

#' calc_relative_ph
#'
#' @description The function to calculate the relative peak heights in one locus
#' @param cand_peaks A vector of the candidate detectable peaks
#' @param gt A vector of the genotype in one locus
#' @param ae The amplification efficiency in one locus
#' @param allele_sr A vector of the alleles for stutter ratios
#' @param bsr_average A vector of the average back stutter ratios
#' @param fsr_average A vector of the average forward stutter ratios
#' @param dsr_average A vector of the average double-back stutter ratios
#' @param m2sr_average A vector of the average minus 2-nt stutter ratios
calc_relative_ph <- function(cand_peaks, gt, ae, allele_sr, bsr_average, fsr_average, dsr_average, m2sr_average){
  rph_al_1 <- rph_bs_1 <- rph_fs_1 <- rph_ds_1 <- rph_m2s_1 <- rep(0, length(cand_peaks))
  rph_al_2 <- rph_bs_2 <- rph_fs_2 <- rph_ds_2 <- rph_m2s_2 <- rep(0, length(cand_peaks))

  for(i in 1:2){
    al <- gt[i]
    index_sr_al <- which(is.element(allele_sr, al))
    if(length(index_sr_al) == 1){
      bsr <- bsr_average[index_sr_al]
      fsr <- fsr_average[index_sr_al]
      dsr <- dsr_average[index_sr_al]
      m2sr <- m2sr_average[index_sr_al]
    }else{
      bsr <- 0
      fsr <- 0
      dsr <- 0
      m2sr <- 0
    }

    index_al <- which(is.element(cand_peaks, al))
    bs <- calc_bs_one(al)
    index_bs <- which(is.element(cand_peaks, bs))
    fs <- calc_fs_one(al)
    index_fs <- which(is.element(cand_peaks, fs))
    ds <- calc_ds_one(al)
    index_ds <- which(is.element(cand_peaks, ds))
    m2s <- calc_m2s_one(al)
    index_m2s <- which(is.element(cand_peaks, m2s))

    if(i == 1){
      rph_al_1[index_al] <- ae / (1 + bsr + fsr + dsr + m2sr)
      if(length(index_bs) == 1){
        rph_bs_1[index_bs] <- ae * bsr / (1 + bsr + fsr + dsr + m2sr)
      }
      if(length(index_fs) == 1){
        rph_fs_1[index_fs] <- ae * fsr / (1 + bsr + fsr + dsr + m2sr)
      }
      if(length(index_ds) == 1){
        rph_ds_1[index_ds] <- ae * dsr / (1 + bsr + fsr + dsr + m2sr)
      }
      if(length(index_m2s) == 1){
        rph_m2s_1[index_m2s] <- ae * m2sr / (1 + bsr + fsr + dsr + m2sr)
      }
    }else{
      rph_al_2[index_al] <- ae / (1 + bsr + fsr + dsr + m2sr)
      if(length(index_bs) == 1){
        rph_bs_2[index_bs] <- ae * bsr / (1 + bsr + fsr + dsr + m2sr)
      }
      if(length(index_fs) == 1){
        rph_fs_2[index_fs] <- ae * fsr / (1 + bsr + fsr + dsr + m2sr)
      }
      if(length(index_ds) == 1){
        rph_ds_2[index_ds] <- ae * dsr / (1 + bsr + fsr + dsr + m2sr)
      }
      if(length(index_m2s) == 1){
        rph_m2s_2[index_m2s] <- ae * m2sr / (1 + bsr + fsr + dsr + m2sr)
      }
    }
  }

  return(rbind(rph_al_1, rph_al_2, rph_bs_1, rph_bs_2, rph_fs_1, rph_fs_2, rph_ds_1, rph_ds_2, rph_m2s_1, rph_m2s_2))
}

#' prepare_gamma_data
#'
#' @description The function to prepare data for estimating the gamma distribution
#' @param markers_auto A vector of the autosomal markers
#' @param peak_list A list of the peaks
#' @param height_list A list of the heights
#' @param gt_list A list of the genotypes
#' @param expect_ae A vector of the expected locus-specific amplification efficiencies
#' @param ave_het_ph_all_data A vector of the average heterozygote allele peak heights
#' @param list_dt_average_sr A list of the average stutter ratios
#' @param sr_consider_auto Information on consideration of back, forward, double-back, and minus 2-nt stutters in autosomal markers
#' @param at_auto A vector of ATs in each autosomal marker
prepare_gamma_data <- function(markers_auto, peak_list, height_list, gt_list, expect_ae, ave_het_ph_all_data, list_dt_average_sr, sr_consider_auto, at_auto){
  samples <- names(gt_list)
  n_sample <- length(samples)
  n_mk <- length(markers_auto)

  sum_peaks_all <- length(unlist(peak_list))
  rphs_gamma_list <- vector(mode = "list", length = n_sample)
  tphs_gamma_list <- vector(mode = "list", length = n_sample)
  ophs_gamma_list <- vector(mode = "list", length = n_sample)
  ats_gamma_list <- vector(mode = "list", length = n_sample)
  samples_gamma_list <- vector(mode = "list", length = n_sample)
  peaks_gamma_list <- vector(mode = "list", length = n_sample)
  markers_gamma_list <- vector(mode = "list", length = n_sample)
  n_adopt_all <- rep(0, n_sample * n_mk)

  count <- 0
  for(i in 1:n_sample){
    peak_s1 <- peak_list[[i]]
    height_s1 <- height_list[[i]]
    gt_s1 <- gt_list[[i]]
    ave_het_ph_s1 <- ave_het_ph_all_data[i]

    rphs_gamma_list_s1 <- vector(mode = "list", length = n_mk)
    names(rphs_gamma_list_s1) <- markers_auto
    tphs_gamma_list_s1 <- vector(mode = "list", length = n_mk)
    names(tphs_gamma_list_s1) <- markers_auto
    ophs_gamma_list_s1 <- vector(mode = "list", length = n_mk)
    names(ophs_gamma_list_s1) <- markers_auto
    ats_gamma_list_s1 <- vector(mode = "list", length = n_mk)
    names(ats_gamma_list_s1) <- markers_auto
    samples_gamma_list_s1 <- vector(mode = "list", length = n_mk)
    names(samples_gamma_list_s1) <- markers_auto
    markers_gamma_list_s1 <- vector(mode = "list", length = n_mk)
    names(markers_gamma_list_s1) <- markers_auto
    peaks_gamma_list_s1 <- vector(mode = "list", length = n_mk)
    names(peaks_gamma_list_s1) <- markers_auto

    for(j in 1:n_mk){
      count <- count + 1
      mk <- markers_auto[j]

      peaks <- peak_s1[[mk]]
      peaks <- round(as.numeric(peaks), 1)
      n_peak <- length(peaks)

      if(n_peak != 0){
        gt <- gt_s1[mk, ]
        gt <- round(as.numeric(gt), 1)

        heights <- height_s1[[mk]]
        ae <- expect_ae[mk]

        dt_average_sr <- list_dt_average_sr[[mk]]
        allele_sr <- round(as.numeric(dt_average_sr$Allele), 1)
        bsr_average <- dt_average_sr[["Average BSR"]]
        fsr_average <- dt_average_sr[["Average FSR"]]
        dsr_average <- dt_average_sr[["Average DSR"]]
        m2sr_average <- dt_average_sr[["Average M2SR"]]
        at <- at_auto[mk]

        cand_peaks <- make_cand_peak(gt, sr_consider_auto[mk, ])

        rphs_per_component <- calc_relative_ph(cand_peaks, gt, ae, allele_sr, bsr_average, fsr_average, dsr_average, m2sr_average)
        rphs <- apply(rphs_per_component, 2, sum)

        index_rph <- which(rphs != 0) # If the stutter ratio = 0, rph for the stutter product is zero.
        if(length(index_rph) == 0){
          rphs_gamma_list_s1[[j]] <- NULL
          tphs_gamma_list_s1[[j]] <- numeric(0)
          ophs_gamma_list_s1[[j]] <- numeric(0)
          ats_gamma_list_s1[[j]] <- numeric(0)
          samples_gamma_list_s1[[j]] <- character(0)
          markers_gamma_list_s1[[j]] <- character(0)
          peaks_gamma_list_s1[[j]] <- numeric(0)
        }else{
          rphs_per_component <- rphs_per_component[, index_rph, drop = FALSE]
          rphs <- rphs[index_rph]
          cand_peaks_tmp <- cand_peaks[index_rph]

          index_gt <- which(is.element(cand_peaks_tmp, gt))
          index_height <- which(is.element(cand_peaks_tmp, peaks))
          index_adopt <- sort(unique(c(index_gt, index_height)))
          n_adopt <- length(index_adopt)

          if(n_adopt == 0){
            rphs_gamma_list_s1[[j]] <- NULL
            tphs_gamma_list_s1[[j]] <- numeric(0)
            ophs_gamma_list_s1[[j]] <- numeric(0)
            ats_gamma_list_s1[[j]] <- numeric(0)
            samples_gamma_list_s1[[j]] <- character(0)
            markers_gamma_list_s1[[j]] <- character(0)
            peaks_gamma_list_s1[[j]] <- numeric(0)
          }else{
            if(all(rphs_per_component[, index_adopt, drop = FALSE] == 0)){
              print(paste0("Sample: ", i, ", Marker: ", j))
            }

            rphs_gamma_list_s1[[j]] <- rphs_per_component[, index_adopt, drop = FALSE]
            tphs_gamma_list_s1[[j]] <- rep(ave_het_ph_s1, n_adopt)
            cand_peaks_adopt <- cand_peaks_tmp[index_adopt]

            ophs <- rep(0, n_adopt)
            index_oph <- which(is.element(peaks, cand_peaks_adopt))
            if(length(index_oph) != 0){
              peaks_adopt <- peaks[index_oph]
              heights_adopt <- heights[index_oph]
              ophs[match(peaks_adopt, cand_peaks_adopt)] <- heights_adopt
            }
            ophs_gamma_list_s1[[j]] <- ophs

            ats_gamma_list_s1[[j]] <- rep(at, n_adopt)
            samples_gamma_list_s1[[j]] <- rep(samples[i], n_adopt)
            markers_gamma_list_s1[[j]] <- rep(mk, n_adopt)
            peaks_gamma_list_s1[[j]] <- cand_peaks_adopt

            n_adopt_all[count] <- n_adopt
          }
        }
      }else{
        rphs_gamma_list_s1[[j]] <- NULL
        tphs_gamma_list_s1[[j]] <- numeric(0)
        ophs_gamma_list_s1[[j]] <- numeric(0)
        ats_gamma_list_s1[[j]] <- numeric(0)
        samples_gamma_list_s1[[j]] <- character(0)
        markers_gamma_list_s1[[j]] <- character(0)
        peaks_gamma_list_s1[[j]] <- numeric(0)
      }
    }
    rphs_gamma_list[[i]] <- rphs_gamma_list_s1
    tphs_gamma_list[[i]] <- tphs_gamma_list_s1
    ophs_gamma_list[[i]] <- ophs_gamma_list_s1
    ats_gamma_list[[i]] <- ats_gamma_list_s1
    samples_gamma_list[[i]] <- samples_gamma_list_s1
    markers_gamma_list[[i]] <- markers_gamma_list_s1
    peaks_gamma_list[[i]] <- peaks_gamma_list_s1
  }

  # Bind rphs
  rphs_gamma <- matrix(NA, 10, sum(n_adopt_all))
  index_start <- 1
  for(i in 1:n_sample){
    for(j in 1:n_mk){
      n_adopt <- n_adopt_all[n_mk * (i - 1) + j]
      if(n_adopt != 0){
        rphs_gamma[, index_start:(index_start + n_adopt - 1)] <- rphs_gamma_list[[i]][[j]]
        index_start <- index_start + n_adopt
      }
    }
  }

  tphs_gamma <- unlist(tphs_gamma_list)
  ophs_gamma <- unlist(ophs_gamma_list)
  ats_gamma <- unlist(ats_gamma_list)
  samples_gamma <- unlist(samples_gamma_list)
  markers_gamma <- unlist(markers_gamma_list)
  peaks_gamma <- unlist(peaks_gamma_list)

  return(list("rphs_gamma" = rphs_gamma,
              "tphs_gamma" = tphs_gamma,
              "ophs_gamma" = ophs_gamma,
              "ats_gamma" = ats_gamma,
              "samples_gamma" = samples_gamma,
              "markers_gamma" = markers_gamma,
              "peaks_gamma" = peaks_gamma))
}

#' run_mle_gamma
#'
#' @description The function to perform MLE for the gamma model
#' @param rphs_gamma The relative peak heights for the gamma distribution
#' @param tphs_gamma The total peak heights in one heterozygote allele for the gamma distribution
#' @param ophs_gamma The observed peak heights for the gamma distribution
#' @param ats_gamma The analytical thresholds for the gamma distribution
#' @param min_scale_allele The minimum scale parameter for alleles
#' @param max_scale_allele The maximum scale parameter for alleles
#' @param min_scale_stutter The minimum scale parameter for stutters
#' @param max_scale_stutter The maximum scale parameter for stutters
#' @param stop_iter Number of iterations without improvement before stopping one MLE run
#' @param max_iter Maximum number of iterations in one MLE run
#' @param num_converged Number of converged MLE runs before stopping
#' @param max_mle_runs Maximum number of MLE attempts
run_mle_gamma <- function(rphs_gamma, tphs_gamma, ophs_gamma, ats_gamma, min_scale_allele, max_scale_allele, min_scale_stutter, max_scale_stutter, stop_iter, max_iter, num_converged, max_mle_runs){

  like_func_gamma <- function(par){
    scale_allele <- par[1] * (max_scale_allele - min_scale_allele) + min_scale_allele
    scale_stutter <- par[2] * (max_scale_stutter - min_scale_stutter) + min_scale_stutter

    n_data <- length(tphs_gamma)
    dens <- rep(Inf, n_data)

    for(i in 1:n_data){
      tmp1 <- rphs_gamma[, i] * tphs_gamma[i]
      tmp2 <- tmp1
      for(j in 1:2){
        tmp2[j] <- tmp2[j] * scale_allele
      }
      for(j in 3:10){
        tmp2[j] <- tmp2[j] * scale_stutter
      }

      tmp1_sum <- sum(tmp1)
      tmp2_sum <- sum(tmp2)
      gamma_shape <- tmp1_sum^2 / tmp2_sum
      gamma_scale <- tmp2_sum / tmp1_sum

      at <- ats_gamma[i]
      oph <- ophs_gamma[i]

      if(oph == 0){
        dens[i] <- pgamma(at, shape = gamma_shape, scale = gamma_scale)
      }else{
        dens[i] <- dgamma(oph, shape = gamma_shape, scale = gamma_scale)
      }
    }

    loglike <- - sum(log(dens))

    if(is.na(loglike) || is.nan(loglike) || is.infinite(loglike)){
      loglike <- Inf
    }

    return(loglike)
  }

  lower_bounds <- c(0, 0)
  upper_bounds <- c(1, 1)

  # Perform simulated annealing
  current_opt <- NULL
  current_score <- Inf
  current_score_round <- Inf
  count_mle <- 0
  count_solution <- 0
  init_list <- list()
  opt_list <- list()
  repeat{
    count_mle <- count_mle + 1

    init_list[[count_mle]] <- init <- c(runif(1, min = min_scale_allele, max = max_scale_allele), runif(1, min = min_scale_stutter, max = max_scale_stutter))

    opt_list[[count_mle]] <- opt <- GenSA(
      par = init,
      fn = like_func_gamma,
      lower = lower_bounds,
      upper = upper_bounds,
      control = list(verbose = TRUE, nb.stop.improvement = stop_iter, maxit = max_iter)
    )

    score <- opt$value
    if(score < current_score){
      current_opt <- opt
      current_score <- score
    }

    score_round <- round(score, 3)
    if(score_round < current_score_round){
      current_score_round <- score_round
      count_solution <- 1
    }else if(score_round == current_score_round){
      count_solution <- count_solution + 1
    }

    if(count_solution == num_converged){
      break
    }else if(count_mle == max_mle_runs){
      current_opt <- NULL
      break
    }
  }

  # Extract parameters
  est_param <- rep(NA, 2)
  names(est_param) <- c("scale_allele", "scale_stutter")
  if(!is.null(current_opt)){
    tmp <- current_opt$par
    est_param[1] <- tmp[1] * (max_scale_allele - min_scale_allele) + min_scale_allele
    est_param[2] <- tmp[2] * (max_scale_stutter - min_scale_stutter) + min_scale_stutter
  }

  return(list(est_param = est_param, init_list = init_list, opt_list = opt_list))
}

#' create_dt_gamma_data
#'
#' @description The function to create the data.table for gamma data
#' @param scale_allele The scale parameter for allelels
#' @param scale_stutter The scale parameter for stutters
#' @param samples_gamma A vector of the sample names for gamma data
#' @param markers_gamma A vector of the markers for gamma data
#' @param peaks_gamma A vector of the peaks for gamma data
#' @param rphs_gamma A matrix of the relative peak heights for gamma data
#' @param tphs_gamma A vector of the total peak heights for gamma data
#' @param ophs_gamma A vector of the observed peak heights for gamma data
#' @param ats_gamma A vector of the ATs for gamma data
create_dt_gamma_data <- function(scale_allele, scale_stutter, samples_gamma, markers_gamma, peaks_gamma, rphs_gamma, tphs_gamma, ophs_gamma, ats_gamma){
  scale_base <- c(rep(scale_allele, 2), rep(scale_stutter, 8))
  index_cpp <- 0:9
  n_data <- length(samples_gamma)
  gamma_shape <- gamma_scale <- quantile95 <- quantile5 <- rep(NA, n_data)
  peak_type <- rep("", n_data)
  for(i in 1:n_data){
    rphs <- rphs_gamma[, i]

    # Categorize peak type
    if(rphs[1] == 0 && rphs[2] == 0){
      peak_type[i] <- "Stutter"
    }else if(any(rphs[3:10] != 0)){
      peak_type[i] <- "Allele + Stutter"
    }else{
      peak_type[i] <- "Allele"
    }

    # Calculate gamma parameters
    gamma_param <- calc_gamma_param_cpp(rphs, scale_base, index_cpp, tphs_gamma[i])
    gamma_shape[i] <- gamma_param[1]
    gamma_scale[i] <- gamma_param[2]

    # Calculate 95% and 5% quantiles
    quantile95[i] <- qgamma(0.95, shape = gamma_param[1], scale = gamma_param[2])
    quantile5[i] <- qgamma(0.05, shape = gamma_param[1], scale = gamma_param[2])
  }

  # Calculate relative peak heights
  rphs_gamma_sum <- apply(rphs_gamma, 2, sum)

  # Calculate expected peak heights
  ephs_gamma <- rphs_gamma_sum * tphs_gamma

  # Create data.table
  dt_gamma_data <- data.table("SampleName" = samples_gamma,
                              "Marker" = markers_gamma,
                              "Peak" = peaks_gamma,
                              "Observed peak height" = ophs_gamma,
                              "Expected peak height" = ephs_gamma,
                              "Relative peak height" = rphs_gamma_sum,
                              "Total peak height" = tphs_gamma,
                              "AT" = ats_gamma,
                              "Shape" = gamma_shape,
                              "Scale" = gamma_scale,
                              "95% quantile" = quantile95,
                              "5% quantile" = quantile5,
                              "Peak type" = peak_type)

  return(dt_gamma_data)
}

#' estimate_gamma
#'
#' @description The function to estimate parameters of the locus-specific amplification efficiencies and the gamma distribution
#' @param dt_exp The data.table for experimental data
#' @param dt_gt The data.table for genotype data
#' @param list_dt_average_sr A list of the data.table for average stutter ratios in all markers
#' @param dt_setting_est_par The data.table for the parameter settings
#' @param dt_expect_ae The data.table for the expected LSAE
#' @param min_scale_allele The minimum scale parameter for alleles
#' @param max_scale_allele The maximum scale parameter for alleles
#' @param min_scale_stutter The minimum scale parameter for stutters
#' @param max_scale_stutter The maximum scale parameter for stutters
#' @param stop_iter Number of iterations without improvement before stopping one MLE run
#' @param max_iter Maximum number of iterations in one MLE run
#' @param num_converged Number of converged MLE runs before stopping
#' @param max_mle_runs Maximum number of MLE attempts
#' @param waiter_ui Information on UI for waiter
estimate_gamma <- function(dt_exp, dt_gt, list_dt_average_sr, dt_setting_est_par, dt_expect_ae,
                           min_scale_allele, max_scale_allele, min_scale_stutter, max_scale_stutter,
                           stop_iter, max_iter, num_converged, max_mle_runs,
                           waiter_ui = NULL){

  # Extract information on setting
  markers_auto <- dt_setting_est_par$Marker
  n_mk_auto <- length(markers_auto)
  at_auto <- dt_setting_est_par$AT
  names(at_auto) <- markers_auto

  # Information on considering stutter ratios
  sr_consider_auto <- matrix(TRUE, n_mk_auto, 4)
  rownames(sr_consider_auto) <- markers_auto
  colnames(sr_consider_auto) <- c("BS", "FS", "DS", "M2S")
  for(i in 1:n_mk_auto){
    mk <- markers_auto[i]
    dt_average_sr <- list_dt_average_sr[[mk]]
    sr_consider_auto[i, "BS"] <- any(dt_average_sr[["Average BSR"]] != 0)
    sr_consider_auto[i, "FS"] <- any(dt_average_sr[["Average FSR"]] != 0)
    sr_consider_auto[i, "DS"] <- any(dt_average_sr[["Average DSR"]] != 0)
    sr_consider_auto[i, "M2S"] <- any(dt_average_sr[["Average M2SR"]] != 0)
  }

  # --- Arrange experimental data --- #

  if(!is.null(waiter_ui)){
    waiter_ui$update(
      html = tagList(
        spin_3k(),
        h3("Preparing data...")
      )
    )
  }

  tmp <- arrange_exp_data(dt_exp, dt_gt, markers_auto, sr_consider_auto, at_auto)
  peak_list <- tmp$peak_list
  height_list <- tmp$height_list
  gt_list <- tmp$gt_list

  # --- Prepare data for the expected LSAE --- #

  expect_ae <- dt_expect_ae$LSAE
  names(expect_ae) <- dt_expect_ae$Marker

  # --- Prepare data for estimating gamma distribution --- #

  ave_het_ph_all_data <- calc_ave_het_ph_all_data(height_list)

  mp <- prepare_gamma_data(markers_auto, peak_list, height_list, gt_list, expect_ae, ave_het_ph_all_data, list_dt_average_sr, sr_consider_auto, at_auto)
  rphs_gamma <- tmp$rphs_gamma
  tphs_gamma <- tmp$tphs_gamma
  ophs_gamma <- tmp$ophs_gamma
  ats_gamma <- tmp$ats_gamma
  samples_gamma <- tmp$samples_gamma
  markers_gamma <- tmp$markers_gamma
  peaks_gamma <- tmp$peaks_gamma

  # --- Perform MLE for the gamma scale parameter --- #

  if(!is.null(waiter_ui)){
    waiter_ui$update(
      html = tagList(
        spin_3k(),
        h3("Performing MLE for peak height parameters...")
      )
    )
  }

  tmp <- run_mle_gamma(rphs_gamma, tphs_gamma, ophs_gamma, ats_gamma, min_scale_allele, max_scale_allele, min_scale_stutter, max_scale_stutter, stop_iter, max_iter, num_converged, max_mle_runs)
  est_param <- tmp$est_param
  scale_allele <- est_param["scale_allele"]
  scale_stutter <- est_param["scale_stutter"]
  init_list_gamma <- tmp$init_list
  opt_list_gamma <- tmp$opt_list

  # --- Create the data.tables --- #

  dt_gamma_data <- create_dt_gamma_data(scale_allele, scale_stutter, samples_gamma, markers_gamma, peaks_gamma, rphs_gamma, tphs_gamma, ophs_gamma, ats_gamma)
  dt_est_gamma <- data.table(Parameter = c("Scale parameter (allele)", "Scale parameter (stutter)"), Value = c(scale_allele, scale_stutter))

  return(list("dt_est_gamma" = dt_est_gamma,
              "init_list_gamma" = init_list_gamma,
              "opt_list_gamma" = opt_list_gamma,
              "dt_gamma_data" = dt_gamma_data))
}
