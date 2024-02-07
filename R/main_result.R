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

make_dt_pg_list <- function(markers_auto, gtc_all, lpd_gtc_all, gtc_logprob_all){
  n_cspmk <- length(gtc_all)
  dt_pg_list <- vector(mode = "list", length = n_cspmk)
  for(i in 1:n_cspmk){
    dt_pg_list[[i]] <- make_dt_pg(gtc_all[[i]], lpd_gtc_all[[i]], gtc_logprob_all[[i]])
  }
  names(dt_pg_list) <- markers_auto
  return(dt_pg_list)
}

calc_likelihood <- function(lpd_gtc_mk, gtc_logprob_mk){
  like_gtc <- lpd_gtc_mk + gtc_logprob_mk
  return(sum(exp(like_gtc)))
}

calc_likelihood_all <- function(lpd_gtc_all, gtc_logprob_all){
  n_mk <- length(lpd_gtc_all)
  likelihood_all <- rep(0, n_mk)
  for(i in 1:n_mk){
    likelihood_all[i] <- calc_likelihood(lpd_gtc_all[[i]], gtc_logprob_all[[i]])
  }
  return(likelihood_all)
}

calc_combined_likelihood <- function(likelihood_all){
  return(exp(sum(log(likelihood_all))))
}

make_dt_likelihood <- function(markers_auto, likelihood_all){
  combined_likelihood <- calc_combined_likelihood(likelihood_all)
  dt_likelihood <- data.table(Marker = c(markers_auto, "Total"),
                              Likelihood = c(likelihood_all, combined_likelihood))
  return(dt_likelihood)
}
