#' make_peak_data
#'
#' @description The function to make peak data
#' @param dt_csp The data.table of the crime stain profiles
#' @param dt_at The data.table of the analytical thresholds
#' @param dt_kit The data.table of the kit information
#' @param markers_auto A vector of the autosomal markers
make_peak_data <- function(dt_csp, dt_at, dt_kit, markers_auto){
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
  setkey(dt_kit, Marker)
  n_mk <- length(markers_auto)
  for(i in 1:length(sn_all_csp)){
    # Define objects for one sample
    peaks_list_sn1 <- list()
    sizes_list_sn1 <- list()
    heights_list_sn1 <- list()

    dt_csp_sn1 <- dt_csp[sn_all_csp[i]]
    setkey(dt_csp_sn1, Marker)
    for(j in 1:n_mk){
      # Extract data for one marker
      mk <- markers_auto[j]
      dt_csp_sn1_mk1 <- dt_csp_sn1[mk]
      peaks <- as.numeric(dt_csp_sn1_mk1[, index_allele, with = FALSE])
      peaks <- peaks[!is.na(peaks)]
      sizes <- as.numeric(dt_csp_sn1_mk1[, index_size, with = FALSE])
      sizes <- sizes[!is.na(sizes)]
      heights <- as.numeric(dt_csp_sn1_mk1[, index_height, with = FALSE])
      heights <- heights[!is.na(heights)]

      # Keep peaks above AT
      dt_at_mk1 <- dt_at[mk]
      at <- as.numeric(dt_at_mk1[, AT])
      index_above_at <- which(heights >= at)
      peaks <- peaks[index_above_at]
      sizes <- sizes[index_above_at]
      heights <- heights[index_above_at]

      # Add information on the drop-out peak
      dt_kit_mk1 <- dt_kit[mk]
      peaks <- c(peaks, 99)
      sizes <- c(sizes, median(dt_kit_mk1[, Size]))
      heights <- c(heights, 0)

      peaks_list_sn1[[j]] <- peaks
      sizes_list_sn1[[j]] <- sizes
      heights_list_sn1[[j]] <- heights
    }

    names(peaks_list_sn1) <- markers_auto
    names(sizes_list_sn1) <- markers_auto
    names(heights_list_sn1) <- markers_auto

    peaks_list[[i]] <- peaks_list_sn1
    sizes_list[[i]] <- sizes_list_sn1
    heights_list[[i]] <- heights_list_sn1
  }

  names(peaks_list) <- sn_all_csp
  names(sizes_list) <- sn_all_csp
  names(heights_list) <- sn_all_csp
  return(list("peaks_list" = peaks_list,
              "sizes_list" = sizes_list,
              "heights_list" = heights_list))
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

#' extract_m2sr_all_mk
#'
#' @description The function to extract minus 2-nt stutter ratios of all markers
#' @param markers_auto A vector of the autosomal markers
#' @param dt_parameter The data.table of the parameters
extract_m2sr_all_mk <- function(markers_auto, dt_parameter){
  setkey(dt_parameter, Parameter)
  dt_parameter_m2sr <- dt_parameter["M2SR"]
  m2sr_all_mk <- dt_parameter_m2sr[, Myu]
  m2sr_all_mk <- m2sr_all_mk[match(markers_auto, dt_parameter_m2sr[, Marker])]
  m2sr_all_mk[is.na(m2sr_all_mk)] <- 0
  return(m2sr_all_mk)
}

#' indicate_gt_al
#'
#' @description The function to make indicators of alleles for genotype combinations
#' @param hnc The hypothesized number of contributors
#' @param peaks A vector of peaks
indicate_gt_al <- function(hnc, peaks){
  # Indices of peaks for candidate genotypes of one contributor
  index_peak <- combinations(n = length(peaks), r = 2, repeats.allowed = TRUE)

  # Indices of genotypes for candidate genotype combinations of each contributor
  index_gt <- permutations(n = nrow(index_peak), r = hnc, repeats.allowed = TRUE)
  index_gt <- as.vector(t(index_gt))

  # Indices of peaks for candidate genotype combinations of each contributor
  index_peak_gt <- index_peak[index_gt , ]
  index_peak_gt <- as.vector(t(index_peak_gt))

  # Indicators of alleles for genotype combinations
  n_index <- length(index_peak_gt)
  ind_gt_al <- matrix(0, n_index, length(peaks))
  colnames(ind_gt_al) <- peaks
  ind_gt_al[cbind(1:n_index, index_peak_gt)] <- 1

  return(ind_gt_al)
}

#' indicate_gt_bs
#'
#' @description The function to make indicators of back stutters for genotype combinations
#' @param ind_gt_al The indicators of alleles for genotype combinations
#' @param peaks A vector of peaks
indicate_gt_bs <- function(ind_gt_al, peaks){
  ind_gt_bs <- ind_gt_al * 0
  bs_all <- calc_bs(peaks)
  index_bs_all <- match(bs_all, peaks)
  for(i in 1:length(peaks)){
    index_bs <- index_bs_all[i]
    if(!is.na(index_bs)){
      ind_gt_bs[, index_bs] <- ind_gt_al[, i]
    }
  }
  return(ind_gt_bs)
}

#' indicate_gt_fs
#'
#' @description The function to make indicators of forward stutters for genotype combinations
#' @param ind_gt_al The indicators of alleles for genotype combinations
#' @param peaks A vector of peaks
indicate_gt_fs <- function(ind_gt_al, peaks){
  ind_gt_fs <- ind_gt_al * 0
  fs_all <- calc_fs(peaks)
  index_fs_all <- match(fs_all, peaks)
  for(i in 1:length(peaks)){
    index_fs <- index_fs_all[i]
    if(!is.na(index_fs)){
      ind_gt_fs[, index_fs] <- ind_gt_al[, i]
    }
  }
  return(ind_gt_fs)
}

#' indicate_gt_ds
#'
#' @description The function to make indicators of double-back stutters for genotype combinations
#' @param ind_gt_al The indicators of alleles for genotype combinations
#' @param peaks A vector of peaks
indicate_gt_ds <- function(ind_gt_al, peaks){
  ind_gt_ds <- ind_gt_al * 0
  ds_all <- calc_ds(peaks)
  index_ds_all <- match(ds_all, peaks)
  for(i in 1:length(peaks)){
    index_ds <- index_ds_all[i]
    if(!is.na(index_ds)){
      ind_gt_ds[, index_ds] <- ind_gt_al[, i]
    }
  }
  return(ind_gt_ds)
}

#' indicate_gt_m2s
#'
#' @description The function to make indicators of minus 2-nt stutters for genotype combinations
#' @param ind_gt_al The indicators of alleles for genotype combinations
#' @param peaks A vector of peaks
indicate_gt_m2s <- function(ind_gt_al, peaks){
  ind_gt_m2s <- ind_gt_al * 0
  m2s_all <- calc_m2s(peaks)
  index_m2s_all <- match(m2s_all, peaks)
  for(i in 1:length(peaks)){
    index_m2s <- index_m2s_all[i]
    if(!is.na(index_m2s)){
      ind_gt_m2s[, index_m2s] <- ind_gt_al[, i]
    }
  }
  return(ind_gt_m2s)
}

#' judge_ind_gt
#'
#' @description The function to judge genotype combinations
#' @param hnc The hypothesized number of contributors
#' @param ind_gt_al The indicators of alleles for genotype combinations
#' @param ind_gt_bs The indicators of back stutters for genotype combinations
#' @param ind_gt_fs The indicators of forward stutters for genotype combinations
#' @param ind_gt_ds The indicators of double-back stutters for genotype combinations
#' @param ind_gt_m2s The indicators of minus 2-nt stutters for genotype combinations
judge_ind_gt <- function(hnc, ind_gt_al, ind_gt_bs, ind_gt_fs, ind_gt_ds, ind_gt_m2s){
  ind_gt_sum <- ind_gt_al + ind_gt_bs + ind_gt_fs + ind_gt_ds + ind_gt_m2s
  ind_gt_sum <- ind_gt_sum[, - ncol(ind_gt_sum), drop = FALSE]

  n_gt_comb <- nrow(ind_gt_sum) / (2 * hnc)
  adopt_index <- integer(0)
  for(i in 1:n_gt_comb){
    current_index <- (2 * hnc * (i - 1) + 1):(2 * hnc * i)
    ind_gt_sum_gtc1 <- ind_gt_sum[current_index, ]
    if(all(apply(ind_gt_sum_gtc1, 2, sum) > 0)){
      adopt_index <- c(adopt_index, current_index)
    }
  }
  return(adopt_index)
}

#' make_indicators_mk
#'
#' @description The function to make indicators in one marker
#' @param hnc The hypothesized number of contributors
#' @param peaks A vector of peaks
#' @param sr_consider The information on consideration of each stutter
make_indicators_mk <- function(hnc, peaks, sr_consider){
  ind_gt_al <- indicate_gt_al(hnc, peaks)

  if(sr_consider[1]){
    ind_gt_bs <- indicate_gt_bs(ind_gt_al, peaks)
  }else{
    ind_gt_bs <- ind_gt_al * 0
  }

  if(sr_consider[2]){
    ind_gt_fs <- indicate_gt_fs(ind_gt_al, peaks)
  }else{
    ind_gt_fs <- ind_gt_al * 0
  }

  if(sr_consider[3]){
    ind_gt_ds <- indicate_gt_ds(ind_gt_al, peaks)
  }else{
    ind_gt_ds <- ind_gt_al * 0
  }

  if(sr_consider[4]){
    ind_gt_m2s <- indicate_gt_m2s(ind_gt_al, peaks)
  }else{
    ind_gt_m2s <- ind_gt_al * 0
  }

  adopt_index <- judge_ind_gt(hnc, ind_gt_al, ind_gt_bs, ind_gt_fs, ind_gt_ds, ind_gt_m2s)
  ind_gt_al <- ind_gt_al[adopt_index, , drop = FALSE]
  ind_gt_bs <- ind_gt_bs[adopt_index, , drop = FALSE]
  ind_gt_fs <- ind_gt_fs[adopt_index, , drop = FALSE]
  ind_gt_ds <- ind_gt_ds[adopt_index, , drop = FALSE]
  ind_gt_m2s <- ind_gt_m2s[adopt_index, , drop = FALSE]

  return(list("ind_gt_al" = ind_gt_al,
              "ind_gt_bs" = ind_gt_bs,
              "ind_gt_fs" = ind_gt_fs,
              "ind_gt_ds" = ind_gt_ds,
              "ind_gt_m2s" = ind_gt_m2s))
}

#' calc_eph_base_mk
#'
#' @description The function to calculate the based expected peak heights in one marker
#' @param ind_gt_al The indicators of alleles for genotype combinations
#' @param ind_gt_bs The indicators of back stutters for genotype combinations
#' @param ind_gt_fs The indicators of forward stutters for genotype combinations
#' @param ind_gt_ds The indicators of double-back stutters for genotype combinations
#' @param ind_gt_m2s The indicators of minus 2-nt stutters for genotype combinations
#' @param peaks A vector of peaks
#' @param dt_alcor_mk1 The data.table of corrected alleles in one marker
#' @param m2sr The average minus 2-nt stutter ratio
calc_eph_base_mk <- function(ind_gt_al, ind_gt_bs, ind_gt_fs, ind_gt_ds, ind_gt_m2s,
                             peaks, dt_alcor_mk1, m2sr){

  alleles <- dt_alcor_mk1[, Allele]
  ave_bsrs <- dt_alcor_mk1[, Average_BSR]
  ave_fsrs <- dt_alcor_mk1[, Average_FSR]
  ave_dsrs <- dt_alcor_mk1[, Average_DSR]

  bsrs <- c(ave_bsrs[match(peaks[-length(peaks)], alleles)], 0)
  fsrs <- c(ave_fsrs[match(peaks[-length(peaks)], alleles)], 0)
  dsrs <- c(ave_dsrs[match(peaks[-length(peaks)], alleles)], 0)

  # For 'Not consider' of each stutter
  bsrs[is.na(bsrs)] <- 0
  fsrs[is.na(fsrs)] <- 0
  dsrs[is.na(dsrs)] <- 0

  bs_all <- calc_bs(peaks)
  fs_all <- calc_fs(peaks)
  ds_all <- calc_ds(peaks)
  m2s_all <- calc_m2s(peaks)

  eph_al <- eph_bs <- eph_fs <- eph_ds <- eph_m2s <- ind_gt_al * 0

  for(i in 1:length(peaks)){
    bsr <- bsrs[i]
    fsr <- fsrs[i]
    dsr <- dsrs[i]

    eph_al[, i] <- ind_gt_al[, i] / (1 + bsr + fsr + dsr + m2sr)

    bs <- bs_all[i]
    index_bs <- match(bs, peaks)
    if(!is.na(index_bs)){
      eph_bs[, index_bs] <- ind_gt_bs[, index_bs] * bsr / (1 + bsr + fsr + dsr + m2sr)
    }

    fs <- fs_all[i]
    index_fs <- match(fs, peaks)
    if(!is.na(index_fs)){
      eph_fs[, index_fs] <- ind_gt_fs[, index_fs] * fsr / (1 + bsr + fsr + dsr + m2sr)
    }

    ds <- ds_all[i]
    index_ds <- match(ds, peaks)
    if(!is.na(index_ds)){
      eph_ds[, index_ds] <- ind_gt_ds[, index_ds] * dsr / (1 + bsr + fsr + dsr + m2sr)
    }

    if(m2sr > 0){
      m2s <- m2s_all[i]
      index_m2s <- match(m2s, peaks)
      if(!is.na(index_m2s)){
        eph_m2s[, index_m2s] <- ind_gt_m2s[, index_m2s] * m2sr / (1 + bsr + fsr + dsr + m2sr)
      }
    }
  }

  eph_base_mk <- eph_al + eph_bs + eph_fs + eph_ds + eph_m2s
  return(eph_base_mk)
}

#' calc_eph_base_all
#'
#' @description The function to calculate the based expected peak heights in all markers
#' @param hnc The hypothesized number of contributors
#' @param markers_auto A vector of the autosomal markers
#' @param peaks_list The list of peaks of all crime stain profiles
#' @param sr_consider_all_mk The information on consideration of each stutter in all markers
#' @param dt_alcor The data.table of corrected alleles
#' @param m2sr_all_mk The average minus 2-nt stutter ratios of all markers
calc_eph_base_all <- function(hnc, markers_auto, peaks_list, sr_consider_all_mk, dt_alcor, m2sr_all_mk){
  setkey(dt_alcor, Marker)
  eph_base_all <- list()
  for(i in 1:length(peaks_list)){
    peaks_list_sn1 <- peaks_list[[i]]
    eph_base_sn1 <- list()
    for(j in 1:length(peaks_list_sn1)){
      peaks <- peaks_list_sn1[[j]]
      tmp <- make_indicators_mk(hnc, peaks, sr_consider_all_mk[j, ])
      ind_gt_al <- tmp$ind_gt_al
      ind_gt_bs <- tmp$ind_gt_bs
      ind_gt_fs <- tmp$ind_gt_fs
      ind_gt_ds <- tmp$ind_gt_ds
      ind_gt_m2s <- tmp$ind_gt_m2s

      dt_alcor_mk1 <- dt_alcor[markers_auto[j]]
      eph_base_sn1[[j]] <- calc_eph_base_mk(ind_gt_al, ind_gt_bs, ind_gt_fs, ind_gt_ds, ind_gt_m2s,
                                            peaks, dt_alcor_mk1, m2sr_all_mk[j])
    }
    eph_base_all[[i]] <- eph_base_sn1
  }
  return(eph_base_all)
}
