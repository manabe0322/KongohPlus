#' calc_loglike_gtc
#'
#' @description The function to calculate the log likelihood of each genotype combination in one marker
#' @param gtc_logprob_mk A vector of the genotype probabilities in one marker
#' @param lpd_gtc_mk_1st A vector of the log likelihood of each genotype in one marker (replicate 1)
#' @param lpd_gtc_mk_2nd A vector of the log likelihood of each genotype in one marker (replicate 2)
calc_loglike_gtc <- function(gtc_logprob_mk, lpd_gtc_mk_1st, lpd_gtc_mk_2nd = NULL){
  if(is.null(lpd_gtc_mk_2nd)){
    loglike_gtc <- gtc_logprob_mk + lpd_gtc_mk_1st
  }else{
    loglike_gtc <- gtc_logprob_mk + lpd_gtc_mk_1st + lpd_gtc_mk_2nd
  }
  return(loglike_gtc)
}

#' calc_log10like_use
#'
#' @description The function to calculate the log10 likelihood of the used markers
#' @param gtc_logprob_use A list of the genotype probabilities in each used marker
#' @param lpd_gtc_use_1st A list of the log likelihood of each genotype in each used marker (replicate 1)
#' @param lpd_gtc_use_2nd A list of the log likelihood of each genotype in each used marker (replicate 2)
calc_log10like_use <- function(gtc_logprob_use, lpd_gtc_use_1st, lpd_gtc_use_2nd = NULL){
  n_mk <- length(lpd_gtc_use_1st)
  log10like_use <- rep(0, n_mk)

  if(is.null(lpd_gtc_use_2nd)){
    lpd_gtc_use_2nd <- vector(mode = "list", length = n_mk)
  }

  for(i in 1:n_mk){
    loglike_gtc <- calc_loglike_gtc(gtc_logprob_use[[i]], lpd_gtc_use_1st[[i]], lpd_gtc_use_2nd[[i]])
    log10like_use[i] <- log10(sum(exp(loglike_gtc)))
  }
  return(log10like_use)
}

#' make_dt_pg
#'
#' @description The function to make a data.table for the result of probabilistic genotyping in one marker
#' @param label_contrib A character vector of the contributor names
#' @param gtc_mk Genotype combinations in one marker
#' @param gtc_logprob_mk A vector of the genotype probabilities in one marker
#' @param lpd_gtc_mk_1st A vector of the log likelihood of each genotype in one marker (replicate 1)
#' @param lpd_gtc_mk_2nd A vector of the log likelihood of each genotype in one marker (replicate 2)
make_dt_pg <- function(label_contrib, gtc_mk, gtc_logprob_mk, lpd_gtc_mk_1st, lpd_gtc_mk_2nd = NULL){
  hnc <- length(gtc_mk[[1]]) / 2
  gtc_mk <- matrix(unlist(gtc_mk), ncol = 2 * hnc, byrow = TRUE)
  gtc_mk[which(gtc_mk == 99, arr.ind = TRUE)] <- "Q"
  gtc_mk_new <- matrix("", nrow(gtc_mk), hnc)

  # Paste two alleles reversing sequence (minor to major -> major to minor)
  for(i in 1:hnc){
    gtc_mk_new[, i] <- apply(gtc_mk[, c(2 * (hnc - i) + 1, 2 * (hnc - i) + 2)], 1, paste, collapse = ", ")
  }

  # Assign NA to lpd_gtc_mk_2nd when performing one PCR
  if(is.null(lpd_gtc_mk_2nd)){
    lpd_gtc_mk_2nd <- rep(NA, length(lpd_gtc_mk_1st))
  }

  # Cut data with genotype probability = 0
  index_keep <- which(gtc_logprob_mk != -Inf)
  gtc_mk_new <- gtc_mk_new[index_keep, , drop = FALSE]
  gtc_logprob_mk <- gtc_logprob_mk[index_keep]
  lpd_gtc_mk_1st <- lpd_gtc_mk_1st[index_keep]
  lpd_gtc_mk_2nd <- lpd_gtc_mk_2nd[index_keep]

  # Calculate posterior probabilities of each genotype combination
  if(any(is.na(lpd_gtc_mk_2nd))){
    loglike_gtc <- calc_loglike_gtc(gtc_logprob_mk, lpd_gtc_mk_1st, NULL)
  }else{
    loglike_gtc <- calc_loglike_gtc(gtc_logprob_mk, lpd_gtc_mk_1st, lpd_gtc_mk_2nd)
  }
  tmp <- exp(loglike_gtc)
  post_prob_gtc <- tmp / sum(tmp)

  # Make a data.table
  mat_pg <- cbind(gtc_mk_new, gtc_logprob_mk, lpd_gtc_mk_1st, lpd_gtc_mk_2nd, post_prob_gtc)
  mat_pg <- mat_pg[gtc_logprob_mk != 0, , drop = FALSE]
  dt_pg <- as.data.table(mat_pg)
  names(dt_pg) <- c(paste0("Gt_", label_contrib), "gtc_prob", "loglike_gtc_1st", "loglike_gtc_2nd", "post_prob")

  # Change data type
  col_numeric <- c("gtc_prob", "loglike_gtc_1st", "loglike_gtc_2nd", "post_prob")
  dt_pg[, (col_numeric) := lapply(.SD, as.numeric), .SDcols = col_numeric]

  # Descending order for weights
  dt_pg <- dt_pg[order(dt_pg[, post_prob], decreasing = TRUE)]

  return(dt_pg)
}

#' make_dt_pg_list
#'
#' @description The function to make a list of data.tables for the result of probabilistic genotyping
#' @param markers_use A character vector of the used markers
#' @param label_contrib A character vector of the contributor names
#' @param gtc_use Genotype combinations in the used markers
#' @param gtc_logprob_use A list of the genotype probabilities in each used marker
#' @param lpd_gtc_use_1st A list of the log likelihood of each genotype in each used marker (replicate 1)
#' @param lpd_gtc_use_2nd A list of the log likelihood of each genotype in each used marker (replicate 2)
make_dt_pg_list <- function(markers_use, label_contrib, gtc_use, gtc_logprob_use, lpd_gtc_use_1st, lpd_gtc_use_2nd = NULL){
  n_mk <- length(gtc_use)
  dt_pg_list <- vector(mode = "list", length = n_mk)

  if(is.null(lpd_gtc_use_2nd)){
    lpd_gtc_use_2nd <- vector(mode = "list", length = n_mk)
  }

  for(i in 1:n_mk){
    dt_pg_list[[i]] <- make_dt_pg(label_contrib, gtc_use[[i]], gtc_logprob_use[[i]], lpd_gtc_use_1st[[i]], lpd_gtc_use_2nd[[i]])
  }
  names(dt_pg_list) <- markers_use
  return(dt_pg_list)
}
