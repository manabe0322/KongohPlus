#' make_prof_ref_all
#'
#' @description The function to make a matrix for reference profiles of all autosomal markers
#' @param dt_ref The data.table of the reference profile
#' @param markers_auto A vector of the autosomal markers
make_prof_ref_all <- function(dt_ref, markers_auto){
  n_mk <- length(markers_auto)
  refs <- unique(dt_ref[, SampleName])
  n_ref <- length(refs)
  setkey(dt_ref, Marker)

  prof_ref_all <- matrix(0, n_mk, 2 * n_ref)
  index_al1 <- 2 * 1:n_ref - 1
  index_al2 <- 2 * 1:n_ref
  for(i in 1:n_mk){
    mk <- markers_auto[i]
    dt_ref_one_mk <- dt_ref[mk]
    prof_ref_all[i, index_al1] <- as.numeric(dt_ref_one_mk[, Allele1])
    prof_ref_all[i, index_al2] <- as.numeric(dt_ref_one_mk[, Allele2])
  }
  rownames(prof_ref_all) <- markers_auto
  colnames(prof_ref_all) <- as.vector(sapply(refs, rep, 2))
  return(prof_ref_all)
}
