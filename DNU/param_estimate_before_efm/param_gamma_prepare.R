#' prepare_eph_oph
#'
#' @description The function to prepare the expected and observed peak heights
#' @param peak_list The list of peaks
#' @param height_list The list of peak heights
#' @param gt_list The list of genotypes
#' @param at_auto Analytical thresholds for autosomal markers
prepare_eph_oph <- function(peak_list, height_list, gt_list, at_auto){
  n_sample <- length(peak_list)
  n_mk <- length(peak_list[[1]])
  eph_all <- oph_all <- at_all <- rep(NA, 2 * n_sample * n_mk)

  for(i in 1:n_sample){
    peak_list_one <- peak_list[[i]]
    height_list_one <- height_list[[i]]
    gt_list_one <- gt_list[[i]]

    for(j in 1:n_mk){
      gt <- gt_list_one[j, ]
      if(gt[1] != gt[2]){
        peaks <- peak_list_one[[j]]
        heights <- height_list_one[[j]]

        index_al <- which(is.element(peaks, gt))
        n_detect_al <- length(index_al)
        if(n_detect_al == 1 || n_detect_al == 2){
          oph <- heights[index_al]
          if(n_detect_al == 1){
            oph <- c(heights[index_al], 0)
          }
          current_index <- 2 * n_mk * (i - 1) + c(2 * j - 1, 2 * j)
          eph_all[current_index] <- mean(oph)
          oph_all[current_index] <- oph
          at_all[current_index] <- at_auto[j]
        }
      }
    }
  }
  eph_all <- eph_all[!is.na(eph_all)]
  oph_all <- oph_all[!is.na(oph_all)]
  at_all <- at_all[!is.na(at_all)]

  return(list("eph_all" = eph_all,
              "oph_all" = oph_all,
              "at_all" = at_all))
}

#' calc_eph_base
#'
#' @description The function to calculate the based expected peak heights
#' @param gt_comb A genotype combination
#' @param peaks A vector of peaks in one marker
#' @param tph The average total peak height of one marker
#' @param myu_ae The parameter 'myu' for AE
#' @param model_bsr The model name of BSR
#' @param myu_bsr The parameter 'myu' for BSR
#' @param slope_bsr The parameter 'slope' for BSR
#' @param intercept_bsr The parameter 'intercept' for BSR
#' @param model_fsr The model name of FSR
#' @param myu_fsr The parameter 'myu' for FSR
#' @param slope_fsr The parameter 'slope' for FSR
#' @param intercept_fsr The parameter 'intercept' for FSR
#' @param model_dsr The model name of DSR
#' @param myu_dsr The parameter 'myu' for DSR
#' @param slope_dsr The parameter 'slope' for DSR
#' @param intercept_dsr The parameter 'intercept' for DSR
#' @param model_m2sr The model name of M2SR
#' @param myu_m2sr The parameter 'myu' for M2SR
calc_eph_base <- function(gt_comb, peaks, tph,
                          myu_ae = 0,
                          model_bsr = "Not consider", myu_bsr = NA, slope_bsr = NA, intercept_bsr = NA,
                          model_fsr = "Not consider", myu_fsr = NA, slope_fsr = NA, intercept_fsr = NA,
                          model_dsr = "Not consider", myu_dsr = NA, slope_dsr = NA, intercept_dsr = NA,
                          model_m2sr = "Not consider", myu_m2sr = NA,
                          info_al_alcor = NULL, alcor_bsr = NULL, alcor_fsr = NULL, alcor_dsr = NULL){
  # Change data type
  gt_comb <- as.numeric(gt_comb)
  peaks <- as.numeric(peaks)

  # Define an output object
  n_peak <- length(peaks)
  hnc <- length(gt_comb) / 2
  eph_base <- matrix(0, hnc, n_peak)
  rownames(eph_base) <- paste0("Contributor", 1:hnc)
  colnames(eph_base) <- peaks
  for(i in 1:hnc){
    # Extract the genotype of one contributor
    gt <- gt_comb[c(2 * i - 1, 2 * i)]
    al <- unique(gt)
    n_al <- length(al)

    # Indices of allele repeat correction
    if(!is.null(info_al_alcor)){
      index_alcor <- match(al, info_al_alcor)
      bool_not_na_alcor <- !is.na(index_alcor)
    }else{
      bool_not_na_alcor <- FALSE
      index_alcor <- NA
    }

    # Back stutter
    if(model_bsr == "Not consider"){
      bs <- bsr <- rep(0, n_al)
    }else if(model_bsr == "Multiple loci together"){
      bs <- round(al - 1, 1)
      bsr <- rep(myu_bsr, n_al)
    }else{
      bs <- round(al - 1, 1)
      bsr <- rep(NA, n_al)
      if(any(bool_not_na_alcor)){
        bsr[bool_not_na_alcor] <- slope_bsr * alcor_bsr[index_alcor[bool_not_na_alcor]] + intercept_bsr
      }
    }

    # Forward stutter
    if(model_fsr == "Not consider"){
      fs <- fsr <- rep(0, n_al)
    }else if(model_fsr == "Multiple loci together"){
      fs <- round(al + 1, 1)
      fsr <- rep(myu_fsr, n_al)
    }else{
      fs <- round(al + 1, 1)
      fsr <- rep(NA, n_al)
      if(any(bool_not_na_alcor)){
        fsr[bool_not_na_alcor] <- slope_fsr * alcor_fsr[index_alcor[bool_not_na_alcor]] + intercept_fsr
      }
    }

    # Double-back stutter
    if(model_dsr == "Not consider"){
      ds <- dsr <- rep(0, n_al)
    }else if(model_dsr == "Multiple loci together"){
      ds <- round(al - 2, 1)
      dsr <- rep(myu_dsr, n_al)
    }else{
      ds <- round(al - 2, 1)
      dsr <- rep(NA, n_al)
      if(any(bool_not_na_alcor)){
        dsr[bool_not_na_alcor] <- slope_dsr * alcor_dsr[index_alcor[bool_not_na_alcor]] + intercept_dsr
      }
    }

    # Minus 2-nt stutter
    if(model_m2sr == "Not consider"){
      m2s <- m2sr <- rep(0, n_al)
    }else{
      m2s <- rep(0, n_al)
      for(j in 1:n_al){
        if(al[j] %% 1 < 0.15){
          m2s[j] <- round(al[j] - 0.8, 1)
        }else{
          m2s[j] <- round(al[j] - 0.2, 1)
        }
      }
      m2sr <- rep(myu_m2sr, n_al)
    }

    # Indices of each allele and stutter
    index_al <- match(al, peaks)
    index_bs <- match(bs, peaks)
    index_fs <- match(fs, peaks)
    index_ds <- match(ds, peaks)
    index_m2s <- match(m2s, peaks)

    # Base peak height
    if(n_al == 1){
      ph_base <- tph * myu_ae
    }else{
      ph_base <- 0.5 * tph * myu_ae
    }

    # Calculate expected peak heights
    ephs <- rep(0, n_peak)
    for(j in 1:n_al){
      if(!is.na(index_al[j])){
        ephs[index_al[j]] <- ephs[index_al[j]] + ph_base / (1 + bsr[j] + fsr[j] + dsr[j] + m2sr[j])
      }
      if(!is.na(index_bs[j])){
        ephs[index_bs[j]] <- ephs[index_bs[j]] + ph_base * bsr[j] / (1 + bsr[j] + fsr[j] + dsr[j] + m2sr[j])
      }
      if(!is.na(index_fs[j])){
        ephs[index_fs[j]] <- ephs[index_fs[j]] + ph_base * fsr[j] / (1 + bsr[j] + fsr[j] + dsr[j] + m2sr[j])
      }
      if(!is.na(index_ds[j])){
        ephs[index_ds[j]] <- ephs[index_ds[j]] + ph_base * dsr[j] / (1 + bsr[j] + fsr[j] + dsr[j] + m2sr[j])
      }
      if(!is.na(index_m2s[j])){
        ephs[index_m2s[j]] <- ephs[index_m2s[j]] + ph_base * m2sr[j] / (1 + bsr[j] + fsr[j] + dsr[j] + m2sr[j])
      }
    }
    eph_base[i, ] <- ephs
  }
  return(eph_base)
}

#' calc_eph_all_data
#'
#' @description The function to calculate the expected peak heights of all data
#' @param markers_auto The vector of autosomal markers
#' @param tph_all_data The average total peak heights of all data
calc_eph_all_data <- function(markers_auto, n_oph, tph_all_data, gt_list, peak_list,
                              dt_result_ae, dt_result_bsr, dt_result_fsr, dt_result_dsr, dt_result_m2sr,
                              dt_alcor_bsr, dt_alcor_fsr, dt_alcor_dsr){
  n_gt <- length(gt_list)
  n_mk <- length(markers_auto)
  eph_all_data <- rep(0, n_oph)

  setkey(dt_alcor_bsr, Marker)
  setkey(dt_alcor_fsr, Marker)
  setkey(dt_alcor_dsr, Marker)

  ae_myu_all_mk <- n_mk * dt_result_ae[, Alpha]

  bsr_model_all_mk <- dt_result_bsr[, Model]
  bsr_myu_all_mk <- dt_result_bsr[, Myu]
  bsr_slope_all_mk <- dt_result_bsr[, Slope]
  bsr_intercept_all_mk <- dt_result_bsr[, Intercept]
  bsr_min_us_all_mk <- dt_result_bsr[, MinUS]

  fsr_model_all_mk <- dt_result_fsr[, Model]
  fsr_myu_all_mk <- dt_result_fsr[, Myu]
  fsr_slope_all_mk <- dt_result_fsr[, Slope]
  fsr_intercept_all_mk <- dt_result_fsr[, Intercept]
  fsr_min_us_all_mk <- dt_result_fsr[, MinUS]

  dsr_model_all_mk <- dt_result_dsr[, Model]
  dsr_myu_all_mk <- dt_result_dsr[, Myu]
  dsr_slope_all_mk <- dt_result_dsr[, Slope]
  dsr_intercept_all_mk <- dt_result_dsr[, Intercept]
  dsr_min_us_all_mk <- dt_result_dsr[, MinUS]

  m2sr_model_all_mk <- dt_result_m2sr[, Model]
  m2sr_myu_all_mk <- dt_result_m2sr[, Myu]

  # Define the initial index of peaks
  index_current <- 1
  for(i in 1:n_gt){
    # Extract one data
    gt_one_data <- gt_list[[i]]
    peaks_one_data <- peak_list[[i]]
    tph <- tph_all_data[i]

    for(j in 1:n_mk){
      mk <- markers_auto[j]
      gt <- gt_one_data[j, ]
      peaks <- peaks_one_data[[j]]
      n_peak <- length(peaks)

      if(n_peak > 0){
        # Extract parameters of AE
        myu_ae <- ae_myu_all_mk[j]

        # Extract parameters of BSR
        adopt_model_bsr <- bsr_model_all_mk[j]
        myu_bsr <- bsr_myu_all_mk[j]
        slope_bsr <- bsr_slope_all_mk[j]
        intercept_bsr <- bsr_intercept_all_mk[j]
        min_us_bsr <- bsr_min_us_all_mk[j]

        # Extract parameters of FSR
        adopt_model_fsr <- fsr_model_all_mk[j]
        myu_fsr <- fsr_myu_all_mk[j]
        slope_fsr <- fsr_slope_all_mk[j]
        intercept_fsr <- fsr_intercept_all_mk[j]
        min_us_fsr <- fsr_min_us_all_mk[j]

        # Extract parameters of DSR
        adopt_model_dsr <- dsr_model_all_mk[j]
        myu_dsr <- dsr_myu_all_mk[j]
        slope_dsr <- dsr_slope_all_mk[j]
        intercept_dsr <- dsr_intercept_all_mk[j]
        min_us_dsr <- dsr_min_us_all_mk[j]

        # Extract parameters of M2SR
        adopt_model_m2sr <- m2sr_model_all_mk[j]
        myu_m2sr <- m2sr_myu_all_mk[j]

        # Assign information on allele repeat correction for BSR
        dt_alcor_bsr_mk <- dt_alcor_bsr[mk]
        if(adopt_model_bsr == "Allele"){
          alcor_bsr <- dt_alcor_bsr_mk[, Allele]
        }else if(is.element(adopt_model_bsr, c("LUS", "Multi-sequence"))){
          alcor_bsr <- dt_alcor_bsr_mk[, Corrected]
        }else{
          alcor_bsr <- NULL
        }

        # Assign information on allele repeat correction for FSR
        dt_alcor_fsr_mk <- dt_alcor_fsr[mk]
        if(adopt_model_fsr == "Allele"){
          alcor_fsr <- dt_alcor_fsr_mk[, Allele]
        }else if(is.element(adopt_model_fsr, c("LUS", "Multi-sequence"))){
          alcor_fsr <- dt_alcor_fsr_mk[, Corrected]
        }else{
          alcor_fsr <- NULL
        }

        # Assign information on allele repeat correction for DSR
        dt_alcor_dsr_mk <- dt_alcor_dsr[mk]
        if(adopt_model_dsr == "Allele"){
          alcor_dsr <- dt_alcor_dsr_mk[, Allele]
        }else if(is.element(adopt_model_dsr, c("LUS", "Multi-sequence"))){
          alcor_dsr <- dt_alcor_dsr_mk[, Corrected]
        }else{
          alcor_dsr <- NULL
        }

        # Information on alleles for allele repeat corrections
        if(is.element(adopt_model_bsr, c("Allele", "LUS", "Multi-sequence"))){
          info_al_one_mk <- dt_alcor_bsr_mk[, Allele]
        }else{
          info_al_one_mk <- NULL
        }

        # Calculate expected peak heights
        eph_base <- calc_eph_base(gt, peaks, tph,
                                  myu_ae,
                                  adopt_model_bsr, myu_bsr, slope_bsr, intercept_bsr,
                                  adopt_model_fsr, myu_fsr, slope_fsr, intercept_fsr,
                                  adopt_model_dsr, myu_dsr, slope_dsr, intercept_dsr,
                                  adopt_model_m2sr, myu_m2sr,
                                  info_al_one_mk, alcor_bsr, alcor_fsr, alcor_dsr)
        eph_base <- apply(eph_base, 2, sum)
        eph_base <- eph_base[match(peaks, as.numeric(names(eph_base)))]
        eph_all_data[index_current:(index_current + n_peak - 1)] <- eph_base

        # Update indices of peaks
        index_current <- index_current + n_peak
      }
    }
  }
  return(eph_all_data)
}
