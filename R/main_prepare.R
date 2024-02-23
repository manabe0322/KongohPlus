#' extract_mp_cand
#'
#' @description The function to extract the candidate mixture proportions
#' @param hnc The hypothesized number of contributors
extract_mp_cand <- function(hnc){
  if(hnc == 1){
    mp_cand <- mp_cand_1p
  }else if(hnc == 2){
    mp_cand <- mp_cand_2p
  }else if(hnc == 3){
    mp_cand <- mp_cand_3p
  }else if(hnc == 4){
    mp_cand <- mp_cand_4p
  }else{
    mp_cand <- NULL
  }
  return(mp_cand)
}

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

#' extract_param_all_mk
#'
#' @description The function to extract minus 2-nt stutter ratios of all markers
#' @param biparam_name A biological parameter name
#' @param extparam_name The extracted parameter name
#' @param markers_auto A vector of the autosomal markers
#' @param dt_param The data.table of the parameters
extract_param_all_mk <- function(biparam_name, extparam_name, markers_auto, dt_param){
  setkey(dt_param, Parameter)
  dt_param_target <- dt_param[biparam_name]
  val_all_mk <- unlist(dt_param_target[, extparam_name, with = FALSE])
  val_all_mk <- val_all_mk[match(markers_auto, dt_param_target[, Marker])]
  val_all_mk[is.na(val_all_mk)] <- 0
  names(val_all_mk) <- markers_auto
  return(val_all_mk)
}

#' extract_beta
#'
#' @description The function to extract the beta
#' @param dt_param The data.table of the parameters
extract_beta <- function(dt_param){
  setkey(dt_param, Parameter)
  dt_gamma <- dt_param["Gamma"]
  return(as.numeric(dt_gamma[, Beta]))
}

#' make_index_dedup2_mk
#'
#' @description The function to make 'index_dedup2_mk'
#' @param index_eph_gtc Information on the connection between eph_base and genotype combinations
#' @param index_dedup1_mk Index for deduplication (pattern 1)
make_index_dedup2_mk <- function(index_eph_gtc, index_dedup1_mk){
  paste_index_eph_gtc <- apply(index_eph_gtc, 1, paste, collapse = "-")
  paste_index_dedup1_mk <- apply(index_dedup1_mk, 1, paste, collapse = "-")
  return(match(paste_index_eph_gtc, paste_index_dedup1_mk))
}
