#' make_at_cspmk
#'
#' @description The function to make a vector of ATs of all samples and all markers
#' @param n_csp The number of crime stain profiles
#' @param markers_auto A vector of the autosomal markers
#' @param dt_at The data.table of the analytical thresholds
make_at_cspmk <- function(n_csp, markers_auto, dt_at){
  markers_at <- dt_at[, Marker]
  ats <- dt_at[, AT]
  at_mk <- ats[match(markers_auto, markers_at)]
  at_cspmk <- rep(at_mk, n_csp)
  return(at_cspmk)
}

#' make_peak_data
#'
#' @description The function to make peak data
#' @param dt_csp The data.table of the crime stain profiles
#' @param dt_kit The data.table of the kit information
#' @param markers_auto A vector of the autosomal markers
#' @param at_cspmk A vector of the analytical thresholds for each marker of all crime stain profiles
make_peak_data <- function(dt_csp, dt_kit, markers_auto, at_cspmk){
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
  markers_dt_kit <- dt_kit[, Marker]
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
      at <- at_cspmk[n_mk * (i - 1) + j]
      index_above_at <- which(heights >= at)
      peaks <- peaks[index_above_at]
      sizes <- sizes[index_above_at]
      heights <- heights[index_above_at]

      # Add information on the drop-out peak
      dt_kit_mk1 <- dt_kit[which(markers_dt_kit == mk)]
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

#' make_ref_data
#'
#' @description The function to make reference data
#' @param dt_ref The data.table of the reference profiles
#' @param markers_auto A vector of the autosomal markers
make_ref_data <- function(dt_ref, markers_auto){
  ref_list <- list()
  index_allele <- grep("Allele", names(dt_ref))
  sn_all_ref <- unique(dt_ref[, SampleName])
  setkey(dt_ref, SampleName)
  n_mk <- length(markers_auto)
  for(i in 1:length(sn_all_ref)){
    ref_list_sn1 <- list()
    dt_ref_sn1 <- dt_ref[sn_all_ref[i]]
    setkey(dt_ref_sn1, Marker)
    for(j in 1:n_mk){
      mk <- markers_auto[j]
      dt_ref_sn1_mk1 <- dt_ref_sn1[mk]
      ref_list_sn1[[j]] <- as.numeric(dt_ref_sn1_mk1[, index_allele, with = FALSE])
    }
    names(ref_list_sn1) <- markers_auto
    ref_list[[i]] <- ref_list_sn1
  }
  names(ref_list) <- sn_all_ref
  return(ref_list)
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

#' extract_myu_all_mk
#'
#' @description The function to extract minus 2-nt stutter ratios of all markers
#' @param param_name A parameter name
#' @param markers_auto A vector of the autosomal markers
#' @param dt_parameter The data.table of the parameters
extract_myu_all_mk <- function(param_name, markers_auto, dt_parameter){
  setkey(dt_parameter, Parameter)
  dt_parameter_target <- dt_parameter[param_name]
  myu_all_mk <- dt_parameter_target[, Myu]
  myu_all_mk <- myu_all_mk[match(markers_auto, dt_parameter_target[, Marker])]
  myu_all_mk[is.na(myu_all_mk)] <- 0
  names(myu_all_mk) <- markers_auto
  return(myu_all_mk)
}

#' extract_beta
#'
#' @description The function to extract the beta
#' @param dt_parameter The data.table of the parameters
extract_beta <- function(dt_parameter){
  setkey(dt_parameter, Parameter)
  dt_gamma <- dt_parameter["Gamma"]
  return(as.numeric(dt_gamma[, Beta]))
}

#' make_index_dedup_2_mk
#'
#' @description The function to make 'index_dedup_2_mk'
#' @param index_eph_gtc Information on the connection between eph_base and genotype combinations
#' @param n_al The number of alleles (i.e., the number of contributors * 2)
#' @param index_dedup_1_mk Index for deduplication (pattern 1)
make_index_dedup_2_mk <- function(index_eph_gtc, n_al, index_dedup_1_mk){
  n_index_eph_gtc <- nrow(index_eph_gtc)
  index_dedup_2_mk <- rep(0, n_index_eph_gtc)
  bool_mat <- matrix(FALSE, nrow(index_dedup_1_mk), n_al)
  for(i in 1:n_index_eph_gtc){
    ind_eph_gtc1 <- index_eph_gtc[i, ]
    for(j in 1:n_al){
      bool_mat[, j] <- index_dedup_1_mk[, j] == ind_eph_gtc1[j]
    }
    index_dedup_2_mk[i] <- which(apply(bool_mat, 1, all))
  }
  return(index_dedup_2_mk)
}
