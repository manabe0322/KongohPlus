#' make_ref_data
#'
#' @description The function to make reference data
#' @param dt_ref The data.table of the reference profiles
#' @param markers_auto A vector of the autosomal markers
make_ref_data <- function(dt_ref, markers_auto){
  ref_list <- list()
  index_allele <- grep("Allele", names(dt_ref))
  sn_all_ref <- unique(dt_ref[, SampleName])
  setindex(dt_ref, SampleName)
  n_mk <- length(markers_auto)
  for(i in 1:length(sn_all_ref)){
    ref_list_sn1 <- list()
    dt_ref_sn1 <- dt_ref[SampleName == sn_all_ref[i]]
    setindex(dt_ref_sn1, Marker)
    for(j in 1:n_mk){
      mk <- markers_auto[j]
      dt_ref_sn1_mk1 <- dt_ref_sn1[Marker == mk]
      ref_list_sn1[[j]] <- as.numeric(dt_ref_sn1_mk1[, index_allele, with = FALSE])
    }
    names(ref_list_sn1) <- markers_auto
    ref_list[[i]] <- ref_list_sn1
  }
  names(ref_list) <- sn_all_ref
  return(ref_list)
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
