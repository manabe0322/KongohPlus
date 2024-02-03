calc_weight <- function(lpd_gtc_mk){
  pd_gtc_mk <- exp(lpd_gtc_mk)
  return(pd_gtc_mk / sum(pd_gtc_mk))
}

make_dt_pg <- function(gtc_mk, lpd_gtc_mk, gtc_logprob_mk){
  hnc <- ncol(gtc_mk) / 2
  gtc_mk[which(gtc_mk == 99, arr.ind = TRUE)] <- "Q"
  gtc_mk_new <- matrix("", nrow(gtc_mk), hnc)

  # Paste two alleles
  for(i in 1:hnc){
    gtc_mk_new[, i] <- apply(gtc_mk[, c(2 * i - 1, 2 * i)], 1, paste, collapse = ", ")
  }

  # Calculate weights
  weight_mk <- calc_weight(lpd_gtc_mk)

  # Make a data.table
  mat_pg <- cbind(gtc_mk_new, weight_mk, lpd_gtc_mk, gtc_logprob_mk)
  dt_pg <- as.data.table(mat_pg)
  names(dt_pg) <- c(paste0("Genotype (U", 1:hnc, ")"), "Weight", "lpd_gtc_mk", "gtc_logprob_mk")

  # Change data type
  col_numeric <- c("Weight", "lpd_gtc_mk", "gtc_logprob_mk")
  dt_pg[, (col_numeric) := lapply(.SD, as.numeric), .SDcols = col_numeric]

  # Descending order for weights
  dt_pg <- dt_pg[order(dt_pg[, Weight], decreasing = TRUE)]

  return(dt_pg)
}

make_dt_pg_list <- function(gtc_all, lpd_gtc_all, gtc_logprob_all){
  n_cspmk <- length(gtc_all)
  dt_pg_list <- vector(mode = "list", length = n_cspmk)
  for(i in 1:n_cspmk){
    dt_pg_list[[i]] <- make_dt_pg(gtc_all[[i]], lpd_gtc_all[[i]], gtc_logprob_all[[i]])
  }
  return(dt_pg_list)
}
