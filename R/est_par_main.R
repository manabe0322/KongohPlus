#' check_est_par
#'
#' @description The function to check files for estimating parameters
#' @param dt_exp The data.table of the experimental data
#' @param dt_gt The data.table of the genotype data
#' @param dt_setting_est_par The data.table in the setting tab for estimating parameters
check_est_par <- function(dt_exp, dt_gt, dt_setting_est_par){
  error_message <- ""

  # Whether all required data is loaded or not
  if(is.null(dt_exp)){
    error_message <- paste0(error_message, "Load experimental data!<br>")
  }
  if(is.null(dt_gt)){
    error_message <- paste0(error_message, "Load genotype data!<br>")
  }

  # Check the analytical thresholds
  markers_est_par <- dt_setting_est_par$Marker
  at_all_mk <- dt_setting_est_par$AT
  bool_without_at <- is.na(at_all_mk)
  if(any(bool_without_at)){
    error_message <- paste0(error_message, "Set the analytical threshold(s) of the following marker(s): ", paste0(markers_est_par[bool_without_at], collapse = ", "), ".<br>")
  }

  if(all(!is.null(dt_exp), !is.null(dt_gt))){
    # Check the marker set of the experimental data
    sn_mk_problem_exp <- check_dt_mk(dt_exp, markers_est_par) # The function is defined in est_sr_main.R
    if(length(sn_mk_problem_exp) > 0){
      error_message <- paste0(error_message, "Inappropriate marker set of the following sample(s) in the experimental data: ", paste0(sn_mk_problem_exp, collapse = ", "), ".<br>")
    }

    # Check the marker set of the genotype data
    sn_mk_problem_gt <- check_dt_mk(dt_gt, markers_est_par) # The function is defined in est_sr_main.R
    if(length(sn_mk_problem_gt) > 0){
      error_message <- paste0(error_message, "Inappropriate marker set of the following sample(s) in the genotype data: ", paste0(sn_mk_problem_gt, collapse = ", "), ".<br>")
    }

    # Check the sample names
    sn_exp <- unique(dt_exp$SampleName)
    sn_gt <- unique(dt_gt$SampleName)
    sn_without_gt <- setdiff(sn_exp, sn_gt)
    if(length(sn_without_gt) > 0){
      error_message <- paste0(error_message, "No genotype data of the following samples: ", paste0(sn_without_gt, collapse = ", "), ".<br>")
    }
  }

  return(error_message)
}

#' create_dt_ae_data
#'
#' @description The function to create a data.table of the locus-specific amplification efficiencies
#' @param peak_list A list of peaks
#' @param height_list A list of heights
create_dt_ae_data <- function(peak_list, height_list){
  sample_names <- names(peak_list)
  n_sample <- length(sample_names)
  markers <- names(peak_list[[1]])
  n_mk <- length(markers)

  n_data <- n_sample * n_mk
  sample_name_all <- rep("", n_data)
  marker_all <- rep("", n_data)
  temp_amount_all <- rep(NA, n_data)
  prop_ae_all <- rep(NA, n_data)
  ae_all <- rep(NA, n_data)

  current_index <- 1
  for(i in 1:n_sample){
    peak_list_one <- peak_list[[i]]
    height_list_one <- height_list[[i]]

    index_assign <- current_index:(current_index + n_mk - 1)
    sample_name_all[index_assign] <- sample_names[i]
    marker_all[index_assign] <- markers

    temp_amount <- sapply(height_list_one, sum)
    sum_temp_amount <- sum(temp_amount)
    temp_amount_all[index_assign] <- temp_amount
    prop_ae_all[index_assign] <- temp_amount / sum_temp_amount
    ae_all[index_assign] <- n_mk * temp_amount / sum_temp_amount

    current_index <- current_index + n_mk
  }

  dt_ae_data <- data.table(SampleName = sample_name_all,
                           Marker = marker_all,
                           TempAmount = temp_amount_all,
                           PropAE = prop_ae_all,
                           AE = ae_all)

  return(dt_ae_data)
}

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

#' run_mle_ae
#'
#' @description The function to perform MLE for the locus-specific amplification efficiency
#' @param prop_ae_data_use A matrix of the proportion of the amplification efficiencies (only data used for MCMC)
#' @param min_alpha The minimum alpha
#' @param max_alpha The maximum alpha
#' @param max_call The max.call parameter for GenSA
#' @param num_converged_solutions The number of converged solutions required to terminate the MLE
#' @param max_mle_runs The maximum number of MLE attempts with different initial values
run_mle_ae <- function(prop_ae_data_use, min_alpha, max_alpha, max_call, num_converged_solutions, max_mle_runs){

  # Likelihood function for the locus-specific amplification efficiency
  like_fun_ae <- function(par){

    n_data <- ncol(prop_ae_data_use)
    dens <- rep(Inf, n_data)

    for(i in 1:n_data){
      dens[i] <- ddirichlet(prop_ae_data_use[, i], par)
    }

    loglike <- - sum(log(dens))

    if(is.na(loglike) || is.nan(loglike) || is.infinite(loglike)){
      loglike <- Inf
    }

    return(loglike)
  }

  n_mk <- nrow(prop_ae_data_use)
  lower_bounds <- rep(min_alpha, n_mk)
  upper_bounds <- rep(max_alpha, n_mk)

  # Perform simulated annealing
  current_opt <- NULL
  current_score <- Inf
  count_mle <- 0
  count_solution <- 0
  init_list <- list()
  record_accept <- rep(NA, max_mle_runs)
  repeat{
    count_mle <- count_mle + 1

    init_list[[count_mle]] <- init <- rdirichlet(1, rep(1, n_mk))

    opt <- GenSA(
      par = init,
      fn = like_fun_ae,
      lower = lower_bounds,
      upper = upper_bounds,
      control = list(verbose = TRUE, max.call = max_call)
    )

    score <- round(opt$value, 3)
    if(score < current_score){
      current_opt <- opt
      current_score <- score
      count_solution <- 1
    }else if(score == current_score){
      count_solution <- count_solution + 1
    }
    record_accept[count_mle] <- count_solution

    if(count_solution == num_converged_solutions){
      break
    }else if(count_mle == max_mle_runs){
      current_opt <- NULL
      break
    }
  }

  # Extract parameters
  if(!is.null(current_opt)){
    est_param <- current_opt$par
  }else{
    est_param <- rep(NA, n_mk)
  }
  names(est_param) <- paste0("alpha[", 1:n_mk, "]")

  return(est_param)
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
  n_sample <- length(gt_list)
  n_mk <- length(markers_auto)

  sum_peaks_all <- length(unlist(peak_list))
  rphs_gamma_list <- vector(mode = "list", length = n_sample)
  tphs_gamma_list <- vector(mode = "list", length = n_sample)
  ophs_gamma_list <- vector(mode = "list", length = n_sample)
  ats_gamma_list <- vector(mode = "list", length = n_sample)
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

    for(j in 1:n_mk){
      count <- count + 1

      peaks <- peak_s1[[j]]
      peaks <- round(as.numeric(peaks), 1)
      n_peak <- length(peaks)
      if(n_peak != 0){
        gt <- gt_s1[j, ]
        gt <- round(as.numeric(gt), 1)

        heights <- height_s1[[j]]
        ae <- expect_ae[j]
        mk <- markers_auto[j]
        dt_average_sr <- list_dt_average_sr[[mk]]
        allele_sr <- round(as.numeric(dt_average_sr$Allele), 1)
        bsr_average <- dt_average_sr[["Average BSR"]]
        fsr_average <- dt_average_sr[["Average FSR"]]
        dsr_average <- dt_average_sr[["Average DSR"]]
        m2sr_average <- dt_average_sr[["Average M2SR"]]
        at <- at_auto[mk]

        cand_peaks <- make_cand_peak(gt, sr_consider_auto[j, ])

        rphs_per_component <- calc_relative_ph(cand_peaks, gt, ae, allele_sr, bsr_average, fsr_average, dsr_average, m2sr_average)
        rphs <- apply(rphs_per_component, 2, sum)

        index_rph <- which(rphs != 0) # If the stutter ratio = 0, rph for the stutter product is zero.
        if(length(index_rph) == 0){
          rphs_gamma_list_s1[[j]] <- NULL
          tphs_gamma_list_s1[[j]] <- numeric(0)
          ophs_gamma_list_s1[[j]] <- numeric(0)
          ats_gamma_list_s1[[j]] <- numeric(0)
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

            n_adopt_all[count] <- n_adopt
          }
        }
      }else{
        rphs_gamma_list_s1[[j]] <- NULL
        tphs_gamma_list_s1[[j]] <- numeric(0)
        ophs_gamma_list_s1[[j]] <- numeric(0)
        ats_gamma_list_s1[[j]] <- numeric(0)
      }
    }
    rphs_gamma_list[[i]] <- rphs_gamma_list_s1
    tphs_gamma_list[[i]] <- tphs_gamma_list_s1
    ophs_gamma_list[[i]] <- ophs_gamma_list_s1
    ats_gamma_list[[i]] <- ats_gamma_list_s1
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

  return(list("rphs_gamma" = rphs_gamma,
              "tphs_gamma" = tphs_gamma,
              "ophs_gamma" = ophs_gamma,
              "ats_gamma" = ats_gamma))
}

#' run_mle_gamma
#'
#' @description The function to perform MLE for the gamma model
#' @param rphs_gamma The relative peak heights for the gamma distribution
#' @param tphs_gamma The total peak heights in one heterozygote allele for the gamma distribution
#' @param ophs_gamma The observed peak heights for the gamma distribution
#' @param ats_gamma The analytical thresholds for the gamma distribution
#' @param min_cv_allele The minimum CV for alleles
#' @param max_cv_allele The maximum CV for alleles
#' @param min_cv_stutter The minimum CV for stutters
#' @param max_cv_stutter The maximum CV for stutters
#' @param max_call The max.call parameter for GenSA
#' @param num_converged_solutions The number of converged solutions required to terminate the MLE
#' @param max_mle_runs The maximum number of MLE attempts with different initial values
run_mle_gamma <- function(rphs_gamma, tphs_gamma, ophs_gamma, ats_gamma, min_cv_allele, max_cv_allele, min_cv_stutter, max_cv_stutter, max_call, num_converged_solutions, max_mle_runs){

  like_func_gamma <- function(par){
    cv_allele_squared <- par[1]^2
    cv_stutter_squared <- par[2]^2

    n_data <- length(tphs_gamma)
    dens <- rep(Inf, n_data)

    for(i in 1:n_data){
      tmp1 <- rphs_gamma[, i] * tphs_gamma[i]
      tmp2 <- tmp1
      for(j in 1:2){
        tmp2[j] <- tmp2[j] * cv_allele_squared
      }
      for(j in 3:10){
        tmp2[j] <- tmp2[j] * cv_stutter_squared
      }

      tmp1_sum <- sum(tmp1)
      tmp2_sum <- sum(tmp2)
      gamma_alpha <- tmp1_sum^2 / tmp2_sum
      gamma_beta <- tmp1_sum / tmp2_sum

      at <- ats_gamma[i]
      oph <- ophs_gamma[i]

      if(oph == 0){
        dens[i] <- pgamma(at, shape = gamma_alpha, rate = gamma_beta)
      }else{
        dens[i] <- dgamma(oph, shape = gamma_alpha, rate = gamma_beta)
      }
    }

    loglike <- - sum(log(dens))

    if(is.na(loglike) || is.nan(loglike) || is.infinite(loglike)){
      loglike <- Inf
    }

    return(loglike)
  }

  lower_bounds <- c(min_cv_allele, min_cv_stutter)
  upper_bounds <- c(max_cv_allele, max_cv_stutter)

  # Perform simulated annealing
  current_opt <- NULL
  current_score <- Inf
  count_mle <- 0
  count_solution <- 0
  init_list <- list()
  record_accept <- rep(NA, max_mle_runs)
  repeat{
    count_mle <- count_mle + 1

    init_list[[count_mle]] <- init <- c(runif(1, min = min_cv_allele, max = max_cv_allele), runif(1, min = min_cv_stutter, max = max_cv_stutter))

    opt <- GenSA(
      par = init,
      fn = like_func_gamma,
      lower = lower_bounds,
      upper = upper_bounds,
      control = list(verbose = TRUE, max.call = max_call)
    )

    score <- round(opt$value, 3)
    if(score < current_score){
      current_opt <- opt
      current_score <- score
      count_solution <- 1
    }else if(score == current_score){
      count_solution <- count_solution + 1
    }
    record_accept[count_mle] <- count_solution

    if(count_solution == num_converged_solutions){
      break
    }else if(count_mle == max_mle_runs){
      current_opt <- NULL
      break
    }
  }

  # Extract parameters
  if(!is.null(current_opt)){
    est_param <- current_opt$par
  }else{
    est_param <- rep(NA, 2)
  }
  names(est_param) <- c("cv_allele", "cv_stutter")

  return(est_param)
}

#' estimate_parameters
#'
#' @description The function to estimate parameters of the locus-specific amplification efficiencies and the gamma distribution
#' @param dt_exp The data.table for experimental data
#' @param dt_gt The data.table for genotype data
#' @param list_dt_average_sr A list of the data.table for average stutter ratios in all markers
#' @param dt_setting_est_par The data.table for the parameter settings
#' @param min_alpha The minimum alpha
#' @param max_alpha The maximum alpha
#' @param min_cv_allele The minimum CV for alleles
#' @param max_cv_allele The maximum CV for alleles
#' @param min_cv_stutter The minimum CV for stutters
#' @param max_cv_stutter The maximum CV for stutters
#' @param max_call The max.call parameter for GenSA
#' @param num_converged_solutions The number of converged solutions required to terminate the MLE
#' @param max_mle_runs The maximum number of MLE attempts with different initial values
estimate_parameters <- function(dt_exp, dt_gt, list_dt_average_sr, dt_setting_est_par,
                                min_alpha, max_alpha,
                                min_cv_allele, max_cv_allele, min_cv_stutter, max_cv_stutter,
                                max_call, num_converged_solutions, max_mle_runs){
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

  # Arrange experimental data
  tmp <- arrange_exp_data(dt_exp, dt_gt, markers_auto, sr_consider_auto, at_auto)
  peak_list <- tmp$peak_list
  height_list <- tmp$height_list
  gt_list <- tmp$gt_list

  # Prepare data for AE
  dt_ae_data <- create_dt_ae_data(peak_list, height_list)
  prop_ae_data <- matrix(dt_ae_data$PropAE, nrow = n_mk_auto)
  prop_ae_data_use <- prop_ae_data[, which(apply(prop_ae_data, 2, prod) > 0), drop = FALSE]

  # Perform MCMC for AE
  alpha <- run_mle_ae(prop_ae_data_use, min_alpha, max_alpha, max_call, num_converged_solutions, max_mle_runs)

  if(any(is.na(alpha))){
    cv_allele <- cv_stutter <- NA

  }else{
    # Calculate expected LSAE
    expect_ae <- n_mk_auto * (alpha + 1) / (sum(alpha) + n_mk_auto)
    names(expect_ae) <- markers_auto

    # Prepare data for average heterozygote allele peak heights
    ave_het_ph_all_data <- calc_ave_het_ph_all_data(height_list)

    # Prepare peak height data for estimating gamma distribution
    tmp <- prepare_gamma_data(markers_auto, peak_list, height_list, gt_list, expect_ae, ave_het_ph_all_data, list_dt_average_sr, sr_consider_auto, at_auto)
    rphs_gamma <- tmp$rphs_gamma
    tphs_gamma <- tmp$tphs_gamma
    ophs_gamma <- tmp$ophs_gamma
    ats_gamma <- tmp$ats_gamma

    # Perform MCMC for gamma distribution
    tmp <- run_mle_gamma(rphs_gamma, tphs_gamma, ophs_gamma, ats_gamma, min_cv_allele, max_cv_allele, min_cv_stutter, max_cv_stutter, max_call, num_converged_solutions, max_mle_runs)
    cv_allele <- tmp["cv_allele"]
    cv_stutter <- tmp["cv_stutter"]
  }

  return(list("alpha" = alpha,
              "expect_ae" = expect_ae,
              "cv_allele" = cv_allele,
              "cv_stutter" = cv_stutter))
}
