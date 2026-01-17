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

#' calc_eph_per_component
#'
#' @description The function to calculate the relative peak heights in one locus
#' @param cand_peaks A vector of the candidate detectable peaks
#' @param gt A vector of the genotype in one locus
#' @param temp The template amount
#' @param ae The amplification efficiency in one locus
#' @param allele_sr A vector of the alleles for stutter ratios
#' @param bsr_average A vector of the average back stutter ratios
#' @param fsr_average A vector of the average forward stutter ratios
#' @param dsr_average A vector of the average double-back stutter ratios
#' @param m2sr_average A vector of the average minus 2-nt stutter ratios
calc_eph_per_component <- function(cand_peaks, gt, temp, ae, allele_sr, bsr_average, fsr_average, dsr_average, m2sr_average){
  eph_al_1 <- eph_bs_1 <- eph_fs_1 <- eph_ds_1 <- eph_m2s_1 <- rep(0, length(cand_peaks))
  eph_al_2 <- eph_bs_2 <- eph_fs_2 <- eph_ds_2 <- eph_m2s_2 <- rep(0, length(cand_peaks))

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
      eph_al_1[index_al] <- temp * ae / (1 + bsr + fsr + dsr + m2sr)
      if(length(index_bs) == 1){
        eph_bs_1[index_bs] <- temp * ae * bsr / (1 + bsr + fsr + dsr + m2sr)
      }
      if(length(index_fs) == 1){
        eph_fs_1[index_fs] <- temp * ae * fsr / (1 + bsr + fsr + dsr + m2sr)
      }
      if(length(index_ds) == 1){
        eph_ds_1[index_ds] <- temp * ae * dsr / (1 + bsr + fsr + dsr + m2sr)
      }
      if(length(index_m2s) == 1){
        eph_m2s_1[index_m2s] <- temp * ae * m2sr / (1 + bsr + fsr + dsr + m2sr)
      }
    }else{
      eph_al_2[index_al] <- temp * ae / (1 + bsr + fsr + dsr + m2sr)
      if(length(index_bs) == 1){
        eph_bs_2[index_bs] <- temp * ae * bsr / (1 + bsr + fsr + dsr + m2sr)
      }
      if(length(index_fs) == 1){
        eph_fs_2[index_fs] <- temp * ae * fsr / (1 + bsr + fsr + dsr + m2sr)
      }
      if(length(index_ds) == 1){
        eph_ds_2[index_ds] <- temp * ae * dsr / (1 + bsr + fsr + dsr + m2sr)
      }
      if(length(index_m2s) == 1){
        eph_m2s_2[index_m2s] <- temp * ae * m2sr / (1 + bsr + fsr + dsr + m2sr)
      }
    }
  }

  return(rbind(eph_al_1, eph_al_2, eph_bs_1, eph_bs_2, eph_fs_1, eph_fs_2, eph_ds_1, eph_ds_2, eph_m2s_1, eph_m2s_2))
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
  ephs_gamma_list <- vector(mode = "list", length = n_sample)
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
    temp <- ave_het_ph_all_data[i]

    ephs_gamma_list_s1 <- vector(mode = "list", length = n_mk)
    names(ephs_gamma_list_s1) <- markers_auto
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

        # --- List all peaks that can be observed under the assumed genotype --- #

        cand_peaks <- make_cand_peak(gt, sr_consider_auto[mk, ])

        # --- Calculate the expected peak heights --- #

        eph_per_component <- calc_eph_per_component(cand_peaks, gt, temp, ae, allele_sr, bsr_average, fsr_average, dsr_average, m2sr_average)
        ephs <- apply(eph_per_component, 2, sum)

        index_eph <- which(ephs != 0) # If the stutter ratio = 0, rph for the stutter product is zero.
        if(length(index_eph) == 0){
          ephs_gamma_list_s1[[j]] <- NULL
          ophs_gamma_list_s1[[j]] <- numeric(0)
          ats_gamma_list_s1[[j]] <- numeric(0)
          samples_gamma_list_s1[[j]] <- character(0)
          markers_gamma_list_s1[[j]] <- character(0)
          peaks_gamma_list_s1[[j]] <- numeric(0)
        }else{

          # --- Keep data for the allele in the genotype and the stutter with the observed peak heights --- #

          eph_per_component <- eph_per_component[, index_eph, drop = FALSE]
          ephs <- ephs[index_eph]
          cand_peaks_tmp <- cand_peaks[index_eph]

          index_gt <- which(is.element(cand_peaks_tmp, gt))
          index_height <- which(is.element(cand_peaks_tmp, peaks))
          index_adopt <- sort(unique(c(index_gt, index_height)))
          n_adopt <- length(index_adopt)

          if(n_adopt == 0){
            ephs_gamma_list_s1[[j]] <- NULL
            ophs_gamma_list_s1[[j]] <- numeric(0)
            ats_gamma_list_s1[[j]] <- numeric(0)
            samples_gamma_list_s1[[j]] <- character(0)
            markers_gamma_list_s1[[j]] <- character(0)
            peaks_gamma_list_s1[[j]] <- numeric(0)
          }else{
            ephs_gamma_list_s1[[j]] <- eph_per_component[, index_adopt, drop = FALSE]
            cand_peaks_adopt <- cand_peaks_tmp[index_adopt]

            # --- Assign the observed peak heights --- #

            ophs <- rep(0, n_adopt)
            index_oph <- which(is.element(peaks, cand_peaks_adopt))
            if(length(index_oph) != 0){
              peaks_adopt <- peaks[index_oph]
              heights_adopt <- heights[index_oph]
              ophs[match(peaks_adopt, cand_peaks_adopt)] <- heights_adopt
            }
            ophs_gamma_list_s1[[j]] <- ophs

            # --- Assign other data --- #

            ats_gamma_list_s1[[j]] <- rep(at, n_adopt)
            samples_gamma_list_s1[[j]] <- rep(samples[i], n_adopt)
            markers_gamma_list_s1[[j]] <- rep(mk, n_adopt)
            peaks_gamma_list_s1[[j]] <- cand_peaks_adopt

            n_adopt_all[count] <- n_adopt
          }
        }
      }else{
        ephs_gamma_list_s1[[j]] <- NULL
        ophs_gamma_list_s1[[j]] <- numeric(0)
        ats_gamma_list_s1[[j]] <- numeric(0)
        samples_gamma_list_s1[[j]] <- character(0)
        markers_gamma_list_s1[[j]] <- character(0)
        peaks_gamma_list_s1[[j]] <- numeric(0)
      }
    }
    ephs_gamma_list[[i]] <- ephs_gamma_list_s1
    ophs_gamma_list[[i]] <- ophs_gamma_list_s1
    ats_gamma_list[[i]] <- ats_gamma_list_s1
    samples_gamma_list[[i]] <- samples_gamma_list_s1
    markers_gamma_list[[i]] <- markers_gamma_list_s1
    peaks_gamma_list[[i]] <- peaks_gamma_list_s1
  }

  # --- Bind all eph data --- #

  ephs_gamma <- matrix(NA, 10, sum(n_adopt_all))
  index_start <- 1
  for(i in 1:n_sample){
    for(j in 1:n_mk){
      n_adopt <- n_adopt_all[n_mk * (i - 1) + j]
      if(n_adopt != 0){
        ephs_gamma[, index_start:(index_start + n_adopt - 1)] <- ephs_gamma_list[[i]][[j]]
        index_start <- index_start + n_adopt
      }
    }
  }

  ophs_gamma <- unlist(ophs_gamma_list)
  ats_gamma <- unlist(ats_gamma_list)
  samples_gamma <- unlist(samples_gamma_list)
  markers_gamma <- unlist(markers_gamma_list)
  peaks_gamma <- unlist(peaks_gamma_list)

  return(list("ephs_gamma" = ephs_gamma,
              "ophs_gamma" = ophs_gamma,
              "ats_gamma" = ats_gamma,
              "samples_gamma" = samples_gamma,
              "markers_gamma" = markers_gamma,
              "peaks_gamma" = peaks_gamma))
}

#' run_mle_gamma
#'
#' @description The function to perform MLE for the gamma model
#' @param ephs_gamma A matrix of the expected peak heights (per component) for the gamma distribution
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
#' @param max_diff_like Reproducibility threshold (log10-likelihood difference)
#' @param waiter_ui Information on UI for waiter
run_mle_gamma <- function(ephs_gamma, ophs_gamma, ats_gamma, min_scale_allele, max_scale_allele, min_scale_stutter, max_scale_stutter, stop_iter, max_iter, num_converged, max_mle_runs, max_diff_like, waiter_ui = NULL){

  # --- Define the likelihood function --- #

  like_func_gamma <- function(par){
    scale_allele <- par[1] * (max_scale_allele - min_scale_allele) + min_scale_allele
    scale_stutter <- par[2] * (max_scale_stutter - min_scale_stutter) + min_scale_stutter

    n_data <- length(ophs_gamma)
    dens <- rep(Inf, n_data)

    for(i in 1:n_data){
      tmp1 <- ephs_gamma[, i]
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

    log10like <- - sum(log(dens)) / log(10)

    if(is.na(log10like) || is.nan(log10like) || is.infinite(log10like)){
      log10like <- Inf
    }

    return(log10like)
  }

  lower_bounds <- c(0, 0)
  upper_bounds <- c(1, 1)

  # --- Perform MLE --- #

  current_opt <- NULL
  current_score <- Inf
  count_mle <- 0
  count_solution <- 0
  init_list <- list()
  opt_list <- list()
  index_converged <- integer(0)
  repeat{
    count_mle <- count_mle + 1

    if(!is.null(waiter_ui)){
      waiter_ui$update(
        html = tagList(
          spin_3k(),
          h3("Performing MLE for the gamma scale parameter..."),
          h3(paste0("Running MLE #", count_mle, " of ", max_mle_runs)),
          h3(paste0("Number of convergences so far: ", count_solution))
        )
      )
    }

    init_list[[count_mle]] <- init <- c(runif(1, min = min_scale_allele, max = max_scale_allele), runif(1, min = min_scale_stutter, max = max_scale_stutter))

    opt_list[[count_mle]] <- opt <- GenSA(
      par = init,
      fn = like_func_gamma,
      lower = lower_bounds,
      upper = upper_bounds,
      control = list(verbose = TRUE, nb.stop.improvement = stop_iter, maxit = max_iter)
    )

    score <- opt$value

    if(score <= current_score){

      if(current_score - score > max_diff_like){
        count_solution <- 1
        index_converged <- count_mle
      }else{
        count_solution <- count_solution + 1
        index_converged <- c(index_converged, count_mle)
      }

      current_opt <- opt
      current_score <- score
    }

    else if(score - current_score <= max_diff_like){
      count_solution <- count_solution + 1
      index_converged <- c(index_converged, count_mle)
    }

    if(count_solution == num_converged){
      break
    }else if(count_mle == max_mle_runs){
      current_opt <- NULL
      break
    }
  }

  # --- Create the warning message --- #

  if(count_solution == num_converged){
    warning_message <- NULL
  }else{
    warning_message <- paste0("- A total of ", count_mle, " MLE runs were performed, of which ", count_solution, " converged. This number is below the specified threshold of ", num_converged, " converged MLE runs , which may compromise the reliability of the results.<br>")
  }

  # --- Calculate the estimated parameters --- #

  est_param <- rep(NA, 2)
  names(est_param) <- c("scale_allele", "scale_stutter")
  if(!is.null(current_opt)){
    tmp <- current_opt$par
    est_param[1] <- tmp[1] * (max_scale_allele - min_scale_allele) + min_scale_allele
    est_param[2] <- tmp[2] * (max_scale_stutter - min_scale_stutter) + min_scale_stutter
  }

  return(list(est_param = est_param, init_list = init_list, opt_list = opt_list, index_converged = index_converged, warning_message = warning_message))
}

#' create_dt_gamma_data
#'
#' @description The function to create the data.table for gamma data
#' @param scale_allele The scale parameter for allelels
#' @param scale_stutter The scale parameter for stutters
#' @param samples_gamma A vector of the sample names for gamma data
#' @param markers_gamma A vector of the markers for gamma data
#' @param peaks_gamma A vector of the peaks for gamma data
#' @param ephs_gamma A matrix of the expected peak heights (per component) for the gamma distribution
#' @param ophs_gamma A vector of the observed peak heights for gamma data
#' @param ats_gamma A vector of the ATs for gamma data
create_dt_gamma_data <- function(scale_allele, scale_stutter, samples_gamma, markers_gamma, peaks_gamma, ephs_gamma, ophs_gamma, ats_gamma){

  scale_base <- c(rep(scale_allele, 2), rep(scale_stutter, 8))
  index_cpp <- 0:9
  n_data <- length(samples_gamma)
  gamma_shape <- gamma_scale <- quantile95 <- quantile5 <- rep(NA, n_data)
  peak_type <- rep("", n_data)

  for(i in 1:n_data){
    ephs <- ephs_gamma[, i]

    # --- Categorize the peak type --- #

    if(ephs[1] == 0 && ephs[2] == 0){
      peak_type[i] <- "Stutter"
    }else if(any(ephs[3:10] != 0)){
      peak_type[i] <- "Allele + Stutter"
    }else{
      peak_type[i] <- "Allele"
    }

    # --- Calculate the gamma parameters --- #

    gamma_param <- calc_gamma_param_cpp(ephs, scale_base, index_cpp)
    gamma_shape[i] <- gamma_param[1]
    gamma_scale[i] <- gamma_param[2]

    # --- Calculate the 95% and 5% quantiles --- #

    quantile95[i] <- qgamma(0.95, shape = gamma_param[1], scale = gamma_param[2])
    quantile5[i] <- qgamma(0.05, shape = gamma_param[1], scale = gamma_param[2])
  }

  # --- Sum up all expected allelic and stutter products in each peak --- #

  ephs_gamma <- apply(ephs_gamma, 2, sum)

  # --- Create the data.table for data of the gamma distribution --- #

  dt_gamma_data <- data.table("SampleName" = samples_gamma,
                              "Marker" = markers_gamma,
                              "Peak" = peaks_gamma,
                              "Observed peak height" = ophs_gamma,
                              "Expected peak height" = ephs_gamma,
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
#' @param max_diff_like Reproducibility threshold (log10-likelihood difference)
#' @param waiter_ui Information on UI for waiter
estimate_gamma <- function(dt_exp, dt_gt, list_dt_average_sr, dt_setting_est_par, dt_expect_ae,
                           min_scale_allele, max_scale_allele, min_scale_stutter, max_scale_stutter,
                           stop_iter, max_iter, num_converged, max_mle_runs, max_diff_like,
                           waiter_ui = NULL){

  # --- Extract the specified conditions --- #

  markers_auto <- dt_setting_est_par$Marker
  n_mk_auto <- length(markers_auto)
  at_auto <- dt_setting_est_par$AT
  names(at_auto) <- markers_auto

  # --- Information on considering the stutter ratios --- #

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

  tmp <- prepare_gamma_data(markers_auto, peak_list, height_list, gt_list, expect_ae, ave_het_ph_all_data, list_dt_average_sr, sr_consider_auto, at_auto)
  ephs_gamma <- tmp$ephs_gamma
  ophs_gamma <- tmp$ophs_gamma
  ats_gamma <- tmp$ats_gamma
  samples_gamma <- tmp$samples_gamma
  markers_gamma <- tmp$markers_gamma
  peaks_gamma <- tmp$peaks_gamma

  # --- Perform MLE for the gamma scale parameter --- #

  tmp <- run_mle_gamma(ephs_gamma, ophs_gamma, ats_gamma, min_scale_allele, max_scale_allele, min_scale_stutter, max_scale_stutter, stop_iter, max_iter, num_converged, max_mle_runs, max_diff_like, waiter_ui)
  est_param <- tmp$est_param
  scale_allele <- est_param["scale_allele"]
  scale_stutter <- est_param["scale_stutter"]
  init_list <- tmp$init_list
  opt_list <- tmp$opt_list
  index_converged <- tmp$index_converged
  warning_message <- tmp$warning_message

  # --- Create the data.tables --- #

  dt_gamma_data <- create_dt_gamma_data(scale_allele, scale_stutter, samples_gamma, markers_gamma, peaks_gamma, ephs_gamma, ophs_gamma, ats_gamma)
  dt_est_gamma <- data.table(Parameter = c("Scale parameter (allele)", "Scale parameter (stutter)"), Value = c(scale_allele, scale_stutter))

  return(list("dt_est_gamma" = dt_est_gamma,
              "init_list" = init_list,
              "opt_list" = opt_list,
              "index_converged" = index_converged,
              "warning_message" = warning_message,
              "dt_gamma_data" = dt_gamma_data))
}
