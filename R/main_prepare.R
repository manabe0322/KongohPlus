#' extract_at_all_mk
#'
#' @description The function to extract ATs of all markers
#' @param markers_auto A vector of the autosomal markers
#' @param dt_at The data.table of the analytical thresholds
extract_at_all_mk <- function(markers_auto, dt_at){
  markers_dt <- dt_at[, Marker]
  ats_dt <- dt_at[, AT]
  at_all_mk <- ats_dt[match(markers_auto, markers_dt)]
  names(at_all_mk) <- markers_auto
  return(at_all_mk)
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
#' @param dt_param The data.table of the parameters
make_sr_consider_all_mk <- function(markers_auto, dt_param){
  setkey(dt_param, Parameter)
  dt_param_bsr <- dt_param["BSR"]
  dt_param_fsr <- dt_param["FSR"]
  dt_param_dsr <- dt_param["DSR"]
  dt_param_m2sr <- dt_param["M2SR"]

  index_mk_bsr <- match(markers_auto, dt_param_bsr[, Marker])
  model_bsr <- dt_param_bsr[, Model]
  model_bsr <- model_bsr[index_mk_bsr]
  sr_consider_all_mk_bsr <- model_bsr != "Not consider"

  index_mk_fsr <- match(markers_auto, dt_param_fsr[, Marker])
  model_fsr <- dt_param_fsr[, Model]
  model_fsr <- model_fsr[index_mk_fsr]
  sr_consider_all_mk_fsr <- model_fsr != "Not consider"

  index_mk_dsr <- match(markers_auto, dt_param_dsr[, Marker])
  model_dsr <- dt_param_dsr[, Model]
  model_dsr <- model_dsr[index_mk_dsr]
  sr_consider_all_mk_dsr <- model_dsr != "Not consider"

  index_mk_m2sr <- match(markers_auto, dt_param_m2sr[, Marker])
  model_m2sr <- dt_param_m2sr[, Model]
  model_m2sr <- model_m2sr[index_mk_m2sr]
  sr_consider_all_mk_m2sr <- model_m2sr != "Not consider"

  sr_consider_all_mk <- cbind(sr_consider_all_mk_bsr, sr_consider_all_mk_fsr, sr_consider_all_mk_dsr, sr_consider_all_mk_m2sr)
  rownames(sr_consider_all_mk) <- markers_auto
  colnames(sr_consider_all_mk) <- c("BS", "FS", "DS", "M2S")

  return(sr_consider_all_mk)
}

#' make_index_dedup2_mk
#'
#' @description The function to make 'index_dedup2_mk'
#' @param index_rph_gtc Information on the connection between rph_base and genotype combinations
#' @param index_dedup1_all_mk Index for deduplication (pattern 1)
make_index_dedup2_mk <- function(index_rph_gtc, index_dedup1_all_mk){
  paste_index_rph_gtc <- apply(index_rph_gtc, 2, paste, collapse = "-")
  paste_index_dedup1_all_mk <- apply(index_dedup1_all_mk, 2, paste, collapse = "-")
  return(match(paste_index_rph_gtc, paste_index_dedup1_all_mk))
}
