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
  names(m2sr_all_mk) <- markers_auto
  return(m2sr_all_mk)
}

#' extract_gamma_beta
#'
#' @description The function to extract the 'beta' parameter of the gamma distribution
#' @param dt_parameter The data.table of the parameters
extract_gamma_beta <- function(dt_parameter){
  setkey(dt_parameter, Parameter)
  dt_parameter_gamma <- dt_parameter["Gamma"]
  return(as.numeric(dt_parameter_gamma[, Beta]))
}

#' set_mp
#'
#' @description The function to set mixture proportions
#' @param hnc The hypothesized number of contributors
#' @param mp_cand_c1 A vector of mixture proportions per contributor
set_mp <- function(hnc, mp_cand_c1){
  if(hnc == 1){
    return(matrix(1, 1, 1))
  }else{
    mp_cand_c1 <- round(mp_cand_c1, 5)
    n_mp_cand_c1 <- length(mp_cand_c1)

    # Mixture proportions of the contributors 1 to N - 1
    mp_other_cn <- combinations(n = n_mp_cand_c1, r = hnc - 1, v = mp_cand_c1, repeats.allowed = TRUE)
    mp_other_cn <- mp_other_cn[which(apply(mp_other_cn, 1, sum) < 1), , drop = FALSE]

    # Mixture proportions of the contributor N
    mp_cn <- round(1 - apply(mp_other_cn, 1, sum), 5)

    # Initial mixture proportions of all contributors
    mp <- cbind(mp_other_cn, mp_cn)

    # Cut overlapped mixture proportions
    mp <- t(apply(mp, 1, sort, decreasing = TRUE))
    mp_cut_count <- 1
    repeat{
      index_mp_cut <- setdiff(which(apply(mp, 1, setequal, mp[mp_cut_count, ])), mp_cut_count)
      if(length(index_mp_cut) != 0){
        mp <- mp[-index_mp_cut, , drop = FALSE]
      }
      mp_cut_count <- mp_cut_count + 1
      if(nrow(mp) < mp_cut_count){
        break
      }
    }

    # Cut mixture proportions of the contributors 1 to N - 1 other than components of 'mp_cand_c1'
    for(i in 2:hnc){
      mp <- mp[is.element(mp[, i], mp_cand_c1), , drop = FALSE]
    }

    return(mp)
  }
}
