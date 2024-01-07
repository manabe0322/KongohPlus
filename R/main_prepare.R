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
#' @param markers_auto A vector of the autosomal markers
extract_peak_data <- function(dt_csp, markers_auto){
  peaks_list <- list()
  sizes_list <- list()
  heights_list <- list()

  n_mk <- length(markers_auto)

  names_csp <- names(dt_csp)
  index_allele <- grep("Allele", names_csp)
  index_size <- grep("Size", names_csp)
  index_height <- grep("Height", names_csp)

  sn_all_csp <- unique(dt_csp[, SampleName])
  setkey(dt_csp, SampleName)
  for(i in 1:length(sn_all_csp)){
    peaks_list_s1 <- list()
    sizes_list_s1 <- list()
    heights_list_s1 <- list()

    dt_csp_s1 <- dt_csp[sn_all_csp[i]]
    setkey(dt_csp_s1, Marker)
    for(j in 1:n_mk){
      dt_csp_s1_mk1 <- dt_csp_s1[markers_auto[j]]
      peaks_list_s1[[j]] <- as.numeric(dt_csp_s1_mk1[, index_allele, with = FALSE])
      sizes_list_s1[[j]] <- as.numeric(dt_csp_s1_mk1[, index_size, with = FALSE])
      heights_list_s1[[j]] <- as.numeric(dt_csp_s1_mk1[, index_height, with = FALSE])
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
#' @param sr_consider Information on consideration of back, forward, double-back, and minus 2-nt stutters
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
