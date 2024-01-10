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
    peaks_list_sn1 <- peaks_list[[i]]
    gt_comb_sn1 <- list()
    for(j in 1:length(peaks_list_sn1)){
      gt_comb_sn1[[j]] <- make_gt_comb_mk(hnc, peaks_list_sn1[[j]], sr_consider_all_mk[j, ])
    }
    names(gt_comb_sn1) <- names(peaks_list_sn1)
    gt_comb_all[[i]] <- gt_comb_sn1
  }
  names(gt_comb_all) <- names(peaks_list)
  return(gt_comb_all)
}

calc_eph_base_main <- function(gt_comb, peaks, dt_alcor_mk1, model_m2sr, myu_m2sr, tph){
  n_peak <- length(peaks)
  n_al <- length(gt_comb)
  eph_base <- matrix(0, n_al, n_peak)
  colnames(eph_base) <- peaks

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
      eph_base[i, n_peak] <- ph_base
    }
  }
  return(eph_base)
}

calc_eph_base_all <- function(gt_comb_all, peaks_list, dt_alcor, model_m2sr_all_mk, myu_m2sr_all_mk, tph_all_sn){

}
