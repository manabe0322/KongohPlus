#' check_dt_mk
#'
#' @description The function to extract sample names with inappropriate marker sets from a data.table
#' @param dt The input data.table
#' @param markers The target markers
check_dt_mk <- function(dt, markers){
  # Extract sample names
  setindex(dt, SampleName)
  sample_names <- unique(dt$SampleName)
  n_sample <- length(sample_names)

  # Check sample names with inappropriate marker sets
  bool_sn_problem <- rep(FALSE, n_sample)
  for(i in 1:n_sample){
    dt_one_sample <- dt[SampleName == sample_names[i]]
    markers_one_sample <- dt_one_sample$Marker
    bool_sn_problem[i] <- !all(is.element(markers, markers_one_sample))
  }
  return(sample_names[bool_sn_problem])
}

#' check_est_sr
#'
#' @description The function to check files for estimating average stutter ratios
#' @param dt_exp The data.table of the experimental data
#' @param dt_gt The data.table of the genotype data
#' @param dt_kit The data.table of information on the kit
#' @param dt_setting_est_sr The data.table in the setting tab for estimating average stutter ratios
check_est_sr <- function(dt_exp, dt_gt, dt_kit, dt_setting_est_sr){
  error_message <- ""
  markers_est_sr <- dt_setting_est_sr$Marker

  # Check the marker set of the experimental data
  sn_mk_problem_exp <- check_dt_mk(dt_exp, markers_est_sr)
  if(length(sn_mk_problem_exp) > 0){
    error_message <- paste0(error_message, "Inappropriate marker set of the following sample(s) in the experimental data: ", paste0(sn_mk_problem_exp, collapse = ", "), ".<br>")
  }

  # Check the marker set of the genotype data
  sn_mk_problem_gt <- check_dt_mk(dt_gt, markers_est_sr)
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

  # Check the analytical thresholds
  at_all_mk <- dt_setting_est_sr$AT
  bool_without_at <- is.na(at_all_mk)
  if(any(bool_without_at)){
    error_message <- paste0(error_message, "Set the analytical threshold(s) of the following marker(s): ", paste0(markers_est_sr[bool_without_at], collapse = ", "), ".<br>")
  }

  # Check the BSR model
  model_bsr <- dt_setting_est_sr[["BSR model"]]
  index_without_model_bsr <- which(model_bsr == "")
  if(length(index_without_model_bsr) != 0){
    error_message <- paste0(error_message, "Set the BSR model of the following marker(s): ", paste0(markers_est_sr[index_without_model_bsr], collapse = ", "), ".<br>")
  }

  # Check the FSR model
  model_fsr <- dt_setting_est_sr[["FSR model"]]
  index_without_model_fsr <- which(model_fsr == "")
  if(length(index_without_model_fsr) != 0){
    error_message <- paste0(error_message, "Set the FSR model of the following marker(s): ", paste0(markers_est_sr[index_without_model_fsr], collapse = ", "), ".<br>")
  }

  # Check the DSR model
  model_dsr <- dt_setting_est_sr[["DSR model"]]
  index_without_model_dsr <- which(model_dsr == "")
  if(length(index_without_model_dsr) != 0){
    error_message <- paste0(error_message, "Set the DSR model of the following marker(s): ", paste0(markers_est_sr[index_without_model_dsr], collapse = ", "), ".<br>")
  }

  # Check the M2SR model
  model_m2sr <- dt_setting_est_sr[["M2SR model"]]
  index_without_model_m2sr <- which(model_m2sr == "")
  if(length(index_without_model_m2sr) != 0){
    error_message <- paste0(error_message, "Set the M2SR model of the following marker(s): ", paste0(markers_est_sr[index_without_model_m2sr], collapse = ", "), ".<br>")
  }

  # Check the registration of all alleles
  n_mk <- length(markers_est_sr)
  dt_gt_tmp <- copy(dt_gt)
  setindex(dt_gt_tmp, Marker)
  setindex(dt_kit, Marker)
  error_message_info_al <- character(0)
  for(i in 1:n_mk){
    mk <- markers_est_sr[i]
    dt_gt_tmp_mk <- dt_gt_tmp[Marker == mk]
    gt_alleles <- unique(c(dt_gt_tmp_mk$Allele1, dt_gt_tmp_mk$Allele2))

    dt_kit_mk <- dt_kit[Marker == mk]
    kit_alleles <- dt_kit_mk[["Allele"]]

    index_without_info_al <- which(!is.element(gt_alleles, kit_alleles))
    if(length(index_without_info_al) != 0){
      if(length(error_message_info_al) == 0){
        error_message_info_al <- "The following alleles are not included in the kit information: <br>"
      }
      error_message_info_al <- paste0(error_message_info_al, "  - ", paste(gt_alleles[index_without_info_al], collapse = ", "), " in ", mk, "<br>")
    }
  }
  error_message <- paste0(error_message, error_message_info_al)

  # Check corrected alleles
  error_message_info_ca <- character(0)
  for(i in 1:n_mk){
    if(is.element("Linear (Corrected Allele)", c(model_bsr[i], model_fsr[i], model_dsr[i], model_m2sr[i]))){
      mk <- markers_est_sr[i]
      dt_kit_mk <- dt_kit[Marker == mk]
      kit_alleles <- dt_kit_mk[["Allele"]]
      kit_cas <- dt_kit_mk[["Corrected Allele"]]
      index_without_info_ca <- which(is.na(kit_cas))
      if(length(index_without_info_ca) != 0){
        if(length(error_message_info_ca) == 0){
          error_message_info_ca <- "The corrected alleles for the following alleles are not included in the kit information: <br>"
        }
        error_message_info_ca <- paste0(error_message_info_ca, "  - ", paste(kit_alleles[index_without_info_ca], collapse = ", "), " in ", mk, "<br>")
      }
    }
  }
  error_message <- paste0(error_message, error_message_info_ca)

  return(error_message)
}

#' extract_sr_data
#'
#' @description The function to extract stutter ratios
#' @param peaks A character vector of peaks
#' @param heights A numeric vector of heights
#' @param gt A character vector of the genotype in one marker
#' @param sr_consider Information on consideration of back, forward, double-back, and minus 2-nt stutters in one marker
extract_sr_data <- function(peaks, heights, gt, sr_consider){
  peaks <- as.numeric(peaks)
  n_peak <- length(peaks)
  gt <- as.numeric(gt)
  gt <- sort(unique(gt[!is.na(gt)]))
  n_al <- length(gt)

  parent_al <- stutter <- rep("", n_peak)
  stutter_ph <- parent_ph <- total_al_prod <- stutter_ratio <- rep(NA, n_peak)
  bool_al <- bool_bs <- bool_fs <- bool_ds <- bool_m2s <- bool_shared_allele <- bool_shared_stutter <- rep(FALSE, n_peak)

  # At least one peak
  if(n_peak > 0){
    bool_al <- is.element(peaks, gt)
    if(sr_consider["BS"]){
      bs <- calc_bs(gt)
      bool_bs <- is.element(peaks, bs)
    }
    if(sr_consider["FS"]){
      fs <- calc_fs(gt)
      bool_fs <- is.element(peaks, fs)
    }
    if(sr_consider["DS"]){
      ds <- calc_ds(gt)
      bool_ds <- is.element(peaks, ds)
    }
    if(sr_consider["M2S"]){
      m2s <- calc_m2s(gt)
      bool_m2s <- is.element(peaks, m2s)
    }

    bool_shared_allele_base <- bool_al & bool_bs + bool_fs + bool_ds + bool_m2s > 0
    bool_shared_stutter <- bool_bs + bool_fs + bool_ds + bool_m2s > 1

    for(i in 1:n_peak){
      if(!bool_al[i]){
        peak <- peaks[i]
        if(!bool_shared_stutter[i]){
          if(bool_bs[i]){
            al <- calc_parent_al_bs_one(peak)
          }else if(bool_fs[i]){
            al <- calc_parent_al_fs_one(peak)
          }else if(bool_ds[i]){
            al <- calc_parent_al_ds_one(peak)
          }else if(bool_m2s[i]){
            al <- calc_parent_al_m2s_one(peak)
          }else{
            al <- NA
          }

          index_al <- which(is.element(peaks, al))

          parent_al[i] <- as.character(al)
          stutter[i] <- as.character(peak)
          stutter_ph[i] <- heights[i]

          # The parent allele is above AT
          if(length(index_al) == 1){
            parent_ph[i] <- heights[index_al]
            total_al_prod[i] <- stutter_ph[i] + parent_ph[i]
            stutter_ratio[i] <- stutter_ph[i] / parent_ph[i]
            bool_shared_allele[i] <- bool_shared_allele_base[index_al]
          }
        }

        # Shared stutters
        else{
          parent_al[i] <- "Multiple"
          stutter[i] <- as.character(peak)
          stutter_ph[i] <- heights[i]
        }
      }
    }

    parent_al <- parent_al[!bool_al]
    stutter <- stutter[!bool_al]
    stutter_ph <- stutter_ph[!bool_al]
    parent_ph <- parent_ph[!bool_al]
    total_al_prod <- total_al_prod[!bool_al]
    stutter_ratio <- stutter_ratio[!bool_al]
    bool_bs <- bool_bs[!bool_al]
    bool_fs <- bool_fs[!bool_al]
    bool_ds <- bool_ds[!bool_al]
    bool_m2s <- bool_m2s[!bool_al]
    bool_shared_allele <- bool_shared_allele[!bool_al]
  }

  return(list("parent_al" = parent_al,
              "stutter" = stutter,
              "stutter_ph" = stutter_ph,
              "parent_ph" = parent_ph,
              "total_al_prod" = total_al_prod,
              "stutter_ratio" = stutter_ratio,
              "bool_bs" = bool_bs,
              "bool_fs" = bool_fs,
              "bool_ds" = bool_ds,
              "bool_m2s" = bool_m2s,
              "bool_shared_allele" = bool_shared_allele))
}

#' create_dt_sr_data
#'
#' @description The function to create a data.table of the stutter ratios
#' @param peak_list A list of peaks
#' @param height_list A list of heights
#' @param gt_list A list of genotypes
#' @param sr_consider_auto Information on consideration of back, forward, double-back, and minus 2-nt stutters in autosomal markers
#' @param dt_kit The data.table of information on the kit
create_dt_sr_data <- function(peak_list, height_list, gt_list, sr_consider_auto, dt_kit){
  sample_names <- names(peak_list)
  n_sample <- length(sample_names)
  markers <- names(peak_list[[1]])
  n_mk <- length(markers)
  setindex(dt_kit, Marker)

  n_max_data <- 8 * n_sample * n_mk
  sample_name_all <- rep("", n_max_data)
  marker_all <- rep("", n_max_data)
  parent_al_all <- rep("", n_max_data)
  stutter_all <- rep("", n_max_data)
  stutter_ratio_all <- rep(NA, n_max_data)
  stutter_ph_all <- rep(NA, n_max_data)
  parent_ph_all <- rep(NA, n_max_data)
  total_al_prod_all <- rep(NA, n_max_data)
  bool_bs_all <- rep(FALSE, n_max_data)
  bool_fs_all <- rep(FALSE, n_max_data)
  bool_ds_all <- rep(FALSE, n_max_data)
  bool_m2s_all <- rep(FALSE, n_max_data)
  bool_shared_allele_all <- rep(FALSE, n_max_data)
  bool_info_al_all <- rep(FALSE, n_max_data)

  current_index <- 1
  for(i in 1:n_sample){
    peak_list_one <- peak_list[[i]]
    height_list_one <- height_list[[i]]
    gt_list_one <- gt_list[[i]]

    for(j in 1:n_mk){
      peaks <- peak_list_one[[j]]
      if(length(peaks) > 0){
        heights <- height_list_one[[j]]
        gt <- gt_list_one[j, ]
        sr_consider <- sr_consider_auto[j, ]
        mk <- markers[j]
        dt_kit_mk <- dt_kit[Marker == mk]
        register_info_al <- dt_kit_mk[["Allele"]]

        tmp <- extract_sr_data(peaks, heights, gt, sr_consider)
        parent_al <- tmp$parent_al
        n_data <- length(parent_al)

        if(n_data > 0){
          index_assign <- current_index:(current_index + n_data - 1)
          sample_name_all[index_assign] <- sample_names[i]
          marker_all[index_assign] <- mk
          parent_al_all[index_assign] <- parent_al
          stutter_all[index_assign] <- tmp$stutter
          stutter_ratio_all[index_assign] <- tmp$stutter_ratio
          stutter_ph_all[index_assign] <- tmp$stutter_ph
          parent_ph_all[index_assign] <- tmp$parent_ph
          total_al_prod_all[index_assign] <- tmp$total_al_prod
          bool_bs_all[index_assign] <- tmp$bool_bs
          bool_fs_all[index_assign] <- tmp$bool_fs
          bool_ds_all[index_assign] <- tmp$bool_ds
          bool_m2s_all[index_assign] <- tmp$bool_m2s
          bool_shared_allele_all[index_assign] <- tmp$bool_shared_allele
          bool_info_al_all[index_assign] <- is.element(parent_al, register_info_al)

          current_index <- current_index + n_data
        }
      }
    }
  }

  sample_name_all <- sample_name_all[1:(current_index - 1)]
  marker_all <- marker_all[1:(current_index - 1)]
  parent_al_all <- parent_al_all[1:(current_index - 1)]
  stutter_all <- stutter_all[1:(current_index - 1)]
  stutter_ratio_all <- stutter_ratio_all[1:(current_index - 1)]
  stutter_ph_all <- stutter_ph_all[1:(current_index - 1)]
  parent_ph_all <- parent_ph_all[1:(current_index - 1)]
  total_al_prod_all <- total_al_prod_all[1:(current_index - 1)]
  bool_bs_all <- bool_bs_all[1:(current_index - 1)]
  bool_fs_all <- bool_fs_all[1:(current_index - 1)]
  bool_ds_all <- bool_ds_all[1:(current_index - 1)]
  bool_m2s_all <- bool_m2s_all[1:(current_index - 1)]
  bool_shared_allele_all <- bool_shared_allele_all[1:(current_index - 1)]
  bool_info_al_all <- bool_info_al_all[1:(current_index - 1)]

  dt_sr_data <- data.table(SampleName = sample_name_all,
                           Marker = marker_all,
                           Allele = parent_al_all,
                           Stutter = stutter_all,
                           SR = stutter_ratio_all,
                           PH_Stutter = stutter_ph_all,
                           PH_Allele = parent_ph_all,
                           TAP = total_al_prod_all,
                           Bool_BS = bool_bs_all,
                           Bool_FS = bool_fs_all,
                           Bool_DS = bool_ds_all,
                           Bool_M2S = bool_m2s_all,
                           Bool_SharedAllele = bool_shared_allele_all,
                           Bool_InfoAllele = bool_info_al_all)

  return(dt_sr_data)
}

#' generate_init_sr_gen
#'
#' @description The function to generate initial values for estimating SR parameters in general model
generate_init_sr_gen <- function(){
  slope <- runif(1)
  intercept <- runif(1, min = -1, max = 1)
  sigma <- runif(1, min = 1, max = 10000)
  return(list(list(slope = slope, intercept = intercept, sigma = sigma)))
}

#' run_mle_sr_gen
#'
#' @description The function to perform MLE for the stutter ratio (locus-specific with slope)
#' @param model_sr_gen The stan model of the stutter ratio
#' @param data_list_sr The data list of the stutter ratio
#' @param n_conv The number of successful convergences
#' @param n_max_init The maximum number of initializations
run_mle_sr_gen <- function(model_sr_gen, data_list_sr, n_conv, n_max_init){
  count_all <- 0
  count_conv <- 0
  max_lp <- -Inf
  est_param <- rep(NA, 3)

  # Repeat random initialization and optimization until the condition is satisfied
  repeat{
    count_all <- count_all + 1
    init_vals <- generate_init_sr_gen()
    fit <- model_sr_gen$optimize(data_list_sr, init = init_vals)

    # If convergence is successful, perform the following steps
    if(fit$return_codes() == 0){
      lp <- fit$lp()

      # If the difference between the current likelihood and the maximum likelihood is less than 1e-6, the current run is considered convergent.
      if(abs(max_lp - lp) < 1e-6){
        count_conv <- count_conv + 1
        if(lp - max_lp > 0){
          max_lp <- lp
          est_param <- fit$mle()
        }
      }

      # Reset convergence count to 1 if a new maximum likelihood is found.
      else if(lp - max_lp > 0){
        count_conv <- 1
        max_lp <- lp
        est_param <- fit$mle()
      }

      # Stop the simulation if convergence is successful for the required number of times.
      if(count_conv == n_conv){
        break
      }

      # Stop the simulation if the maximum number of simulations is reached, treating it as a convergence failure.
      if(count_all == n_max_init){
        est_param <- rep(NA, 3)
        break
      }
    }
  }

  names(est_param) <- c("slope", "intercept", "sigma")
  return(est_param)
}

#' generate_init_sr_uni
#'
#' @description The function to generate initial values for estimating SR parameters in uniform model
generate_init_sr_uni <- function(){
  myu <- runif(1)
  sigma <- runif(1, min = 1, max = 10000)
  return(list(list(myu = myu, sigma = sigma)))
}

#' run_mle_sr_uni
#'
#' @description The function to perform MLE for the stutter ratio (uniform)
#' @param model_sr_uni The stan model of the stutter ratio
#' @param data_list_sr The data list of the stutter ratio
#' @param n_conv The number of successful convergences
#' @param n_max_init The maximum number of initializations
run_mle_sr_uni <- function(model_sr_uni, data_list_sr, n_conv, n_max_init){
  count_all <- 0
  count_conv <- 0
  max_lp <- -Inf
  est_param <- rep(NA, 2)

  # Repeat random initialization and optimization until the condition is satisfied
  repeat{
    count_all <- count_all + 1
    init_vals <- generate_init_sr_uni()
    fit <- model_sr_uni$optimize(data_list_sr, init = init_vals)

    # If convergence is successful, perform the following steps
    if(fit$return_codes() == 0){
      lp <- fit$lp()

      # If the difference between the current likelihood and the maximum likelihood is less than 1e-6, the current run is considered convergent.
      if(abs(max_lp - lp) < 1e-6){
        count_conv <- count_conv + 1
        if(lp - max_lp > 0){
          max_lp <- lp
          est_param <- fit$mle()
        }
      }

      # Reset convergence count to 1 if a new maximum likelihood is found.
      else if(lp - max_lp > 0){
        count_conv <- 1
        max_lp <- lp
        est_param <- fit$mle()
      }

      # Stop the simulation if convergence is successful for the required number of times.
      if(count_conv == n_conv){
        break
      }

      # Stop the simulation if the maximum number of simulations is reached, treating it as a convergence failure.
      if(count_all == n_max_init){
        est_param <- rep(NA, 2)
        break
      }
    }
  }

  names(est_param) <- c("myu", "sigma")
  return(est_param)
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

#' perform_mle_sr_all_mk
#'
#' @description The function to perform MLE for SR of all markers
#' @param n_conv The number of successful convergences
#' @param n_max_init The maximum number of initializations
perform_mle_sr_all_mk <- function(dt_sr_data, dt_kit, stutter_type, markers_auto, meth_model_sr, model_sr_gen, model_sr_uni, n_conv, n_max_init){
  # Extract data of the targeted stutter type
  if(stutter_type == "BS"){
    dt_st_type <- dt_sr_data[dt_sr_data$Bool_BS]
  }else if(stutter_type == "FS"){
    dt_st_type <- dt_sr_data[dt_sr_data$Bool_FS]
  }else if(stutter_type == "DS"){
    dt_st_type <- dt_sr_data[dt_sr_data$Bool_DS]
  }else if(stutter_type == "M2S"){
    dt_st_type <- dt_sr_data[dt_sr_data$Bool_M2S]
  }else{
    dt_st_type <- dt_sr_data[rep(FALSE, nrow(dt_sr_data))]
  }

  # Define objects
  n_mk <- length(markers_auto)
  setindex(dt_st_type, Marker, Allele)
  setindex(dt_kit, Marker)
  myu_all_mk <- rep(NA, n_mk)
  slope_all_mk <- rep(NA, n_mk)
  intercept_all_mk <- rep(NA, n_mk)
  sigma_all_mk <- rep(NA, n_mk)
  dt_average_sr_all_mk <- vector(mode = "list", length = n_mk)
  names(dt_average_sr_all_mk) <- markers_auto

  index_multi_mk <- integer(0)
  y_sr_multi_mk <- numeric(0)
  tap_multi_mk <- numeric(0)
  tmp_alname_multi_mk <- vector(mode = "list", length = n_mk)
  tmp_alcor_multi_mk <- vector(mode = "list", length = n_mk)

  # Repeat for the number of markers
  for(i in 1:n_mk){
    mk <- markers_auto[i]
    meth <- meth_model_sr[i]
    if(meth != "Not consider"){

      # Extract stutter data for the current marker
      dt_st_mk <- dt_st_type[Marker == mk]
      dt_st_mk <- dt_st_mk[Allele != "Multiple"]
      dt_st_mk <- dt_st_mk[!dt_st_mk$Bool_SharedAllele]
      dt_kit_mk <- dt_kit[Marker == mk]

      if(nrow(dt_st_mk) != 0){

        # Extract (corrected) allele repeat numbers, stutter ratios, and total allelic products
        if(meth == "Linear (Corrected Allele)"){
          dt_st_mk <- dt_st_mk[dt_st_mk$Bool_InfoAllele]
          st_alleles <- dt_st_mk$Allele

          alname <- dt_kit_mk[["Allele"]]
          alcor <- dt_kit_mk[["Corrected Allele"]]

          index_adopt <- which(is.element(st_alleles, alname))
          st_alleles <- st_alleles[index_adopt]

          index_alcor <- match(st_alleles, alname)
          x_al <- alcor[index_alcor]
          y_sr <- dt_st_mk$SR[index_adopt]
          tap <- dt_st_mk$TAP[index_adopt]
        }else{
          st_alleles <- dt_st_mk$Allele
          x_al <- as.numeric(st_alleles)
          y_sr <- dt_st_mk$SR
          tap <- dt_st_mk$TAP

          tmp_alname <- c(dt_kit_mk[["Allele"]], st_alleles)
          alcor <- sort(unique(as.numeric(tmp_alname)))
          alname <- as.character(alcor)
        }

        # Merge the data of the current marker with data from other markers where Uniform (Multiple loci) is selected.
        if(meth == "Uniform (Multiple loci)"){
          index_multi_mk <- c(index_multi_mk, i)
          y_sr_multi_mk <- c(y_sr_multi_mk, y_sr)
          tap_multi_mk <- c(tap_multi_mk, tap)
          tmp_alname_multi_mk[[i]] <- alname
          tmp_alcor_multi_mk[[i]] <- alcor
        }

        # Perform MLE for uniform model
        else if(meth == "Uniform"){
          data_list_sr <- list(N = length(x_al), SR = y_sr, Height = tap)
          est_param <- run_mle_sr_uni(model_sr_uni, data_list_sr, n_conv, n_max_init)
          myu_all_mk[i] <- myu <- est_param["myu"]
          sigma_all_mk[i] <- est_param["sigma"]

          # Set the average stutter ratios for each allele
          average_sr <- rep(myu, length(alcor))
          if(!is.na(myu)){
            average_sr[average_sr < 0] <- 0
          }
          dt_average_sr_all_mk[[i]] <- data.table(Allele = alname,
                                                  AlCor = alcor,
                                                  AveSR = average_sr)
        }

        # Perform MLE for linear model
        else{
          data_list_sr <- list(N = length(x_al), SR = y_sr, XVal = x_al, Height = tap)
          est_param <- run_mle_sr_gen(model_sr_gen, data_list_sr, n_conv, n_max_init)
          slope_all_mk[i] <- slope <- est_param["slope"]
          intercept_all_mk[i] <- intercept <- est_param["intercept"]
          sigma_all_mk[i] <- est_param["sigma"]

          # Set the average stutter ratios for each allele
          if(any(is.na(slope), is.na(intercept))){
            average_sr <- rep(NA, length(alcor))
          }else{
            average_sr <- calc_average_sr(alcor, slope, intercept)
            average_sr[average_sr < 0] <- 0
          }
          dt_average_sr_all_mk[[i]] <- data.table(Allele = alname,
                                                  AlCor = alcor,
                                                  AveSR = average_sr)
        }
      }else{
        if(meth == "Uniform (Multiple loci)"){
          index_multi_mk <- c(index_multi_mk, i)
          tmp_alname <- dt_kit_mk[["Allele"]]
          alcor <- sort(unique(as.numeric(tmp_alname)))
          alname <- as.character(alcor)
          tmp_alname_multi_mk[[i]] <- alname
          tmp_alcor_multi_mk[[i]] <- alcor
        }
      }
    }else{
      dt_average_sr_all_mk[[i]] <- data.table(Allele = character(0),
                                              AlCor = numeric(0),
                                              AveSR = numeric(0))
    }
  }

  # Perform MLE for uniform model (multiple loci together)
  if(length(y_sr_multi_mk) != 0){
    data_list_sr <- list(N = length(y_sr_multi_mk), SR = y_sr_multi_mk, Height = tap_multi_mk)
    est_param <- run_mle_sr_uni(model_sr_uni, data_list_sr, n_conv, n_max_init)
    myu_all_mk[index_multi_mk] <- est_param["myu"]
    sigma_all_mk[index_multi_mk] <- est_param["sigma"]

    # Set the average stutter ratios for each allele
    for(i in 1:length(index_multi_mk)){
      index_mk <- index_multi_mk[i]
      myu <- myu_all_mk[index_mk]
      alname <- tmp_alname_multi_mk[[index_mk]]
      alcor <- tmp_alcor_multi_mk[[index_mk]]

      average_sr <- rep(myu, length(alcor))
      if(!is.na(myu)){
        average_sr[average_sr < 0] <- 0
      }
      dt_average_sr_all_mk[[index_mk]] <- data.table(Allele = alname,
                                                     AlCor = alcor,
                                                     AveSR = average_sr)
    }
  }

  dt_result_est_sr <- data.table(Marker = markers_auto,
                                 Model = meth_model_sr,
                                 Myu = myu_all_mk,
                                 Slope = slope_all_mk,
                                 Intercept = intercept_all_mk,
                                 Sigma = sigma_all_mk)

  return(list("dt_result_est_sr" = dt_result_est_sr,
              "dt_average_sr_all_mk" = dt_average_sr_all_mk))
}

#' estimate_average_sr
#'
#' @description The function to estimate average stutter ratios
estimate_average_sr <- function(dt_exp, dt_gt, dt_kit, dt_setting_est_sr, model_sr_gen, model_sr_uni, n_conv, n_max_init){
  # Extract information on setting
  markers_auto <- dt_setting_est_sr$Marker
  n_mk_auto <- length(markers_auto)
  at_auto <- dt_setting_est_sr$AT
  names(at_auto) <- markers_auto
  meth_model_bsr <- dt_setting_est_sr[["BSR model"]]
  meth_model_fsr <- dt_setting_est_sr[["FSR model"]]
  meth_model_dsr <- dt_setting_est_sr[["DSR model"]]
  meth_model_m2sr <- dt_setting_est_sr[["M2SR model"]]

  # Information on considering stutter ratios
  sr_consider_auto <- matrix(TRUE, n_mk_auto, 4)
  rownames(sr_consider_auto) <- markers_auto
  colnames(sr_consider_auto) <- c("BS", "FS", "DS", "M2S")
  sr_consider_auto[, 1] <- meth_model_bsr != "Not consider"
  sr_consider_auto[, 2] <- meth_model_fsr != "Not consider"
  sr_consider_auto[, 3] <- meth_model_dsr != "Not consider"
  sr_consider_auto[, 4] <- meth_model_m2sr != "Not consider"

  # Arrange experimental data
  tmp <- arrange_exp_data(dt_exp, dt_gt,
                          markers_auto, sr_consider_auto, at_auto)
  peak_list <- tmp$peak_list
  height_list <- tmp$height_list
  gt_list <- tmp$gt_list

  # Create a data.table for the stutter ratios
  dt_sr_data <- create_dt_sr_data(peak_list, height_list, gt_list, sr_consider_auto, dt_kit)

  # Perform MLE for BSR
  tmp <- perform_mle_sr_all_mk(dt_sr_data, dt_kit, "BS", markers_auto, meth_model_bsr, model_sr_gen, model_sr_uni, n_conv, n_max_init)
  dt_result_est_bsr <- tmp$dt_result_est_sr
  list_dt_average_bsr <- tmp$dt_average_sr_all_mk

  # Perform MLE for FSR
  tmp <- perform_mle_sr_all_mk(dt_sr_data, dt_kit, "FS", markers_auto, meth_model_fsr, model_sr_gen, model_sr_uni, n_conv, n_max_init)
  dt_result_est_fsr <- tmp$dt_result_est_sr
  list_dt_average_fsr <- tmp$dt_average_sr_all_mk

  # Perform MLE for DSR
  tmp <- perform_mle_sr_all_mk(dt_sr_data, dt_kit, "DS", markers_auto, meth_model_dsr, n_conv, n_max_init)
  dt_result_est_dsr <- tmp$dt_result_est_sr
  list_dt_average_dsr <- tmp$dt_average_sr_all_mk

  # Perform MLE for M2SR
  tmp <- perform_mle_sr_all_mk(dt_sr_data, dt_kit, "M2S", markers_auto, meth_model_m2sr, n_conv, n_max_init)
  dt_result_est_m2sr <- tmp$dt_result_est_sr
  list_dt_average_m2sr <- tmp$dt_average_sr_all_mk

  return(list("dt_result_est_bsr" = dt_result_est_bsr,
              "list_dt_average_bsr" = list_dt_average_bsr,
              "dt_result_est_fsr" = dt_result_est_fsr,
              "list_dt_average_fsr" = list_dt_average_fsr,
              "dt_result_est_dsr" = dt_result_est_dsr,
              "list_dt_average_dsr" = list_dt_average_dsr,
              "dt_result_est_m2sr" = dt_result_est_m2sr,
              "list_dt_average_m2sr" = list_dt_average_m2sr))
}
