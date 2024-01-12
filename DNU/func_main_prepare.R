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
