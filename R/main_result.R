calc_weight <- function(lpd_gtc_mk){
  pd_gtc_mk <- exp(lpd_gtc_mk)
  return(pd_gtc_mk / sum(pd_gtc_mk))
}

make_dt_pg <- function(gtc_mk, lpd_gtc_mk, gtc_logprob_mk){
  hnc <- length(gtc_mk[[1]]) / 2
  gtc_mk <- matrix(unlist(gtc_mk), ncol = 2 * hnc, byrow = TRUE)
  gtc_mk[which(gtc_mk == 99, arr.ind = TRUE)] <- "Q"
  gtc_mk_new <- matrix("", nrow(gtc_mk), hnc)

  # Paste two alleles
  for(i in 1:hnc){
    gtc_mk_new[, i] <- apply(gtc_mk[, c(2 * (hnc - i) + 1, 2 * (hnc - i) + 2)], 1, paste, collapse = ", ")
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

make_dt_pg_list <- function(markers_auto, gtc_all_mk, lpd_gtc_all_mk, gtc_logprob_all_mk){
  n_mk <- length(gtc_all_mk)
  dt_pg_list <- vector(mode = "list", length = n_mk)
  for(i in 1:n_mk){
    dt_pg_list[[i]] <- make_dt_pg(gtc_all_mk[[i]], lpd_gtc_all_mk[[i]], gtc_logprob_all_mk[[i]])
  }
  names(dt_pg_list) <- markers_auto
  return(dt_pg_list)
}

calc_log10like_mk <- function(lpd_gtc_mk, gtc_logprob_mk){
  loglike_mk <- lpd_gtc_mk + gtc_logprob_mk
  return(log10(sum(exp(loglike_mk))))
}

calc_log10like_all_mk <- function(lpd_gtc_all_mk, gtc_logprob_all_mk){
  n_mk <- length(lpd_gtc_all_mk)
  log10like_all_mk <- rep(0, n_mk)
  for(i in 1:n_mk){
    log10like_all_mk[i] <- calc_log10like_mk(lpd_gtc_all_mk[[i]], gtc_logprob_all_mk[[i]])
  }
  return(log10like_all_mk)
}

calc_combined_log10like <- function(log10like_all_mk){
  return(sum(log10like_all_mk))
}

make_dt_log10like <- function(markers_auto, log10like_all_mk){
  combined_log10like <- calc_combined_log10like(log10like_all_mk)
  dt_log10like <- data.table(Marker = c(markers_auto, "Total"),
                             Log10_Likelihood = c(log10like_all_mk, combined_log10like))
  return(dt_log10like)
}

make_dt_lr <- function(markers_auto, log10like_hp, log10like_hd){
  dt_lr <- data.table(Marker = c(markers_auto, "Total"),
                      Log10_Likelihood_Hp = log10like_hp,
                      Log10_Likelihood_Hd = log10like_hd,
                      LR = 10^(log10like_hp - log10like_hd))
  return(dt_lr)
}
