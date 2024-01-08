#' extract_auto_mk
#'
#' @description The function to extract autosomal markers from the data.table of the information on the kit
#' @param dt_kit The data.table of the information on the kit
extract_auto_mk <- function(dt_kit){
  setkey(dt_kit, Sex_Chr)
  dt_kit_auto <- dt_kit[""]
  return(dt_kit_auto[, Marker])
}

#' extract_peak_data
#'
#' @description The function to extract peak data
#' @param dt_csp The data.table of the crime stain profiles
#' @param dt_at The data.table of the analytical thresholds
#' @param markers_auto A vector of the autosomal markers
extract_peak_data <- function(dt_csp, dt_at, markers_auto){
  # Define output objects
  peaks_list <- list()
  sizes_list <- list()
  heights_list <- list()

  # Column indices for alleles, sizes, and heights
  names_csp <- names(dt_csp)
  index_allele <- grep("Allele", names_csp)
  index_size <- grep("Size", names_csp)
  index_height <- grep("Height", names_csp)

  sn_all_csp <- unique(dt_csp[, SampleName])
  setkey(dt_csp, SampleName)
  setkey(dt_at, Marker)
  n_mk <- length(markers_auto)
  for(i in 1:length(sn_all_csp)){
    peaks_list_s1 <- list()
    sizes_list_s1 <- list()
    heights_list_s1 <- list()

    dt_csp_s1 <- dt_csp[sn_all_csp[i]]
    setkey(dt_csp_s1, Marker)
    for(j in 1:n_mk){
      mk <- markers_auto[j]
      dt_csp_s1_mk1 <- dt_csp_s1[mk]
      peaks <- as.numeric(dt_csp_s1_mk1[, index_allele, with = FALSE])
      peaks <- peaks[!is.na(peaks)]
      sizes <- as.numeric(dt_csp_s1_mk1[, index_size, with = FALSE])
      sizes <- sizes[!is.na(sizes)]
      heights <- as.numeric(dt_csp_s1_mk1[, index_height, with = FALSE])
      heights <- heights[!is.na(heights)]

      dt_at_mk1 <- dt_at[mk]
      at <- as.numeric(dt_at_mk1[, AT])
      index_above_at <- which(heights >= at)

      peaks_list_s1[[j]] <- peaks[index_above_at]
      sizes_list_s1[[j]] <- sizes[index_above_at]
      heights_list_s1[[j]] <- heights[[index_above_at]]
    }

    names(peaks_list_s1) <- markers_auto
    names(sizes_list_s1) <- markers_auto
    names(heights_list_s1) <- markers_auto

    peaks_list[[i]] <- peaks_list_s1
    sizes_list[[i]] <- sizes_list_s1
    heights_list[[i]] <- heights_list_s1
  }

  names(peaks_list) <- sn_all_csp
  names(sizes_list) <- sn_all_csp
  names(heights_list) <- sn_all_csp
  return("peaks_list" = peaks_list,
         "sizes_list" = sizes_list,
         "heights_list" = heights_list)
}

#' make_sr_consider_all_mk
#'
#' @description The function to make information on consideration of each stutter
#' @param markers_auto A vector of the autosomal markers
#' @param dt_parameter The data.table of the parameters
make_sr_consider_all_mk <- function(markers_auto, dt_parameter){
  setkey(dt_parameter, Parameter)
  dt_parameter_bsr <- dt_parameter["BSR"]
  dt_parameter_fsr <- dt_parameter["FSR"]
  dt_parameter_dsr <- dt_parameter["DSR"]
  dt_parameter_m2sr <- dt_parameter["M2SR"]

  index_mk_bsr <- match(markers_auto, dt_parameter_bsr[, Marker])
  model_bsr <- dt_parameter_bsr[, Model]
  model_bsr <- model_bsr[index_mk_bsr]
  sr_consider_all_mk_bsr <- model_bsr != "Not consider"

  index_mk_fsr <- match(markers_auto, dt_parameter_fsr[, Marker])
  model_fsr <- dt_parameter_fsr[, Model]
  model_fsr <- model_fsr[index_mk_fsr]
  sr_consider_all_mk_fsr <- model_fsr != "Not consider"

  index_mk_dsr <- match(markers_auto, dt_parameter_dsr[, Marker])
  model_dsr <- dt_parameter_dsr[, Model]
  model_dsr <- model_dsr[index_mk_dsr]
  sr_consider_all_mk_dsr <- model_dsr != "Not consider"

  index_mk_m2sr <- match(markers_auto, dt_parameter_m2sr[, Marker])
  model_m2sr <- dt_parameter_m2sr[, Model]
  model_m2sr <- model_m2sr[index_mk_m2sr]
  sr_consider_all_mk_m2sr <- model_m2sr != "Not consider"

  sr_consider_all_mk <- cbind(sr_consider_all_mk_bsr, sr_consider_all_mk_fsr, sr_consider_all_mk_dsr, sr_consider_all_mk_m2sr)
  rownames(sr_consider_all_mk) <- markers_auto
  colnames(sr_consider_all_mk) <- c("BS", "FS", "DS", "M2S")

  return(sr_consider_all_mk)
}

#' judge_gt_comb
#'
#' @description The function to judge a genotype combination
#' @param gt_comb A genotype combination
#' @param peaks A vector of peaks
#' @param sr_consider Information on consideration of back, forward, double-back, and minus 2-nt stutters
judge_gt_comb <- function(gt_comb, peaks, sr_consider){
  cand_peak <- make_cand_peak(gt_comb, sr_consider)
  return(all(is.element(peaks[peaks != 99], cand_peak)))
}

#' make_gt_comb_mk
#'
#' @description The function to make genotype combinations of one marker
#' @param hnc The hypothesized number of contributors
#' @param peaks A vector of peaks in one marker
#' @param sr_consider Information on consideration of back, forward, double-back, and minus 2-nt stutters in one marker
make_gt_comb_mk <- function(hnc, peaks, sr_consider){
  # Add drop-out
  peaks <- c(peaks, 99)

  # Genotypes of one contributor
  gts <- combinations(n = length(peaks), r = 2, v = peaks, repeats.allowed = TRUE)

  # Indices of genotypes of each contributor
  id_gt_comb <- permutations(n = nrow(gts), r = hnc, repeats.allowed = TRUE)

  # Initial gt_comb_mk
  gt_comb_mk <- matrix(0, nrow(id_gt_comb), 2 * hnc)
  for(i in 1:hnc){
    gt_comb_mk[, c(2 * (i - 1) + 1, 2 * i)] <- gts[id_gt_comb[, i], ]
  }

  # Judge each genotype combination
  bool_gt_comb_mk <- apply(gt_comb_mk, 1, judge_gt_comb, peaks, sr_consider)
  gt_comb_mk <- gt_comb_mk[bool_gt_comb_mk, , drop = FALSE]

  return(gt_comb_mk)
}

#' make_gt_comb_all
#'
#' @description The function to make all genotype combinations
#' @param hnc The hypothesized number of contributors
#' @param peaks_list The list of peaks of all crime stain profiles
#' @param sr_consider_all_mk Information on consideration of back, forward, double-back, and minus 2-nt stutters in all markers
make_gt_comb_all <- function(hnc, peaks_list, sr_consider_all_mk){
  gt_comb_all <- list()
  for(i in 1:length(peaks_list)){
    peaks_list_s1 <- peaks_list[[i]]
    gt_comb_s1 <- list()
    for(j in 1:length(peaks_list_s1)){
      gt_comb_s1[[j]] <- make_gt_comb_mk(hnc, peaks_list_s1[[j]], sr_consider_all_mk[j, ])
    }
    names(gt_comb_s1) <- names(peaks_list_s1)
    gt_comb_all[[i]] <- gt_comb_s1
  }
  names(gt_comb_all) <- names(peaks_list)
  return(gt_comb_all)
}

#' calc_m2s
#'
#' @description The function to calculate the allele call of the minus 2-nt stutter
#' @param al The parent allele of the minus 2-nt stutter
calc_m2s <- function(al){
  if(al %% 1 < 0.15){
    m2s <- round(al - 0.8, 1)
  }else{
    m2s <- round(al - 0.2, 1)
  }
  return(m2s)
}

calc_eph_base_tmp <- function(gt_comb, peaks, dt_alcor_mk1, model_m2sr, myu_m2sr, tph){
  n_peak <- length(peaks)
  n_al <- length(gt_comb)
  eph_base <- matrix(0, n_al, n_peak + 1)
  colnames(eph_base) <- c(peaks, 99)

  alleles <- dt_alcor_mk1[, Allele]
  model_bsr <- unique(dt_alcor_mk1[, Model_BSR])
  model_fsr <- unique(dt_alcor_mk1[, Model_FSR])
  model_dsr <- unique(dt_alcor_mk1[, Model_DSR])
  ave_bsrs <- dt_alcor_mk1[, Average_BSR]
  ave_fsrs <- dt_alcor_mk1[, Average_FSR]
  ave_dsrs <- dt_alcor_mk1[, Average_DSR]

  ph_base <- 0.5 * tph
  for(i in 1:n_al){
    al <- gt_comb[i]

    if(al != 99){
      index_al <- which(alleles == al)

      # No sequence data (only for parameter estimation)
      if(length(index_al) == 0){
        eph_base[i, n_peak + 1] <- NA
      }else{
        # Back stutter
        if(model_bsr == "Not consider"){
          bs <- bsr <- 0
        }else{
          bs <- round(al - 1, 1)
          bsr <- ave_bsrs[index_al]
        }

        # Forward stutter
        if(model_fsr == "Not consider"){
          fs <- fsr <- 0
        }else{
          fs <- round(al + 1, 1)
          fsr <- ave_fsrs[index_al]
        }

        # Double-back stutter
        if(model_dsr == "Not consider"){
          ds <- dsr <- 0
        }else{
          ds <- round(al - 2, 1)
          dsr <- ave_dsrs[index_al]
        }

        # Minus 2-nt stutter
        if(model_m2sr == "Not consider"){
          m2s <- m2sr <- 0
        }else{
          m2s <- calc_m2s(al)
          m2sr <- myu_m2sr
        }

        # Indices of each allele and stutter
        index_al <- match(al, peaks)
        index_bs <- match(bs, peaks)
        index_fs <- match(fs, peaks)
        index_ds <- match(ds, peaks)
        index_m2s <- match(m2s, peaks)

        # Calculate expected peak heights
        ephs <- rep(0, n_peak + 1)
        if(!is.na(index_al)){
          ephs[index_al] <- ph_base / (1 + bsr + fsr + dsr + m2sr)
        }
        if(!is.na(index_bs)){
          ephs[index_bs] <- ph_base * bsr / (1 + bsr + fsr + dsr + m2sr)
        }
        if(!is.na(index_fs)){
          ephs[index_fs] <- ph_base * fsr / (1 + bsr + fsr + dsr + m2sr)
        }
        if(!is.na(index_ds)){
          ephs[index_ds] <- ph_base * dsr / (1 + bsr + fsr + dsr + m2sr)
        }
        if(!is.na(index_m2s)){
          ephs[index_m2s] <- ph_base * m2sr / (1 + bsr + fsr + dsr + m2sr)
        }
        eph_base[i, ] <- ephs
      }
    }else{
      eph_base[i, n_peak + 1] <- ph_base
    }
  }
  return(eph_base)
}
