#' search_min_rhat_mpdeg
#'
#' @description The function to search the minimum rhat value to determine the mixture proportions and the degradation parameters
#' @param sampling_data A matrix of the sampling data obtained by the function 'as_draws_matrix'
#' @param n_sampling The number of sampling for MCMC
#' @param n_chain The number of chains for MCMC
#' @param ind_mp The column indices of the mixture proportions in the sampling_data
#' @param ind_deg The column indices of the degradation parameters in the sampling_data
#' @param hnc The hypothesized number of contributors
search_min_rhat_mpdeg <- function(sampling_data, n_sampling, n_chain, ind_mp, ind_deg, hnc){
  # Calculate rhat values of all orders of the mixture proportions
  cand_col_index_per_chain <- permutations(n = hnc, r = hnc, repeats.allowed = FALSE)
  index_access_all <- permutations(n = nrow(cand_col_index_per_chain), r = n_chain, repeats.allowed = TRUE)
  n_access <- nrow(index_access_all)
  rhat_mp_all <- rhat_deg_all <- matrix(0, n_access, hnc)
  for(i in 1:n_access){
    index_access_1 <- index_access_all[i, ]
    sampling_data_mp <- sampling_data_deg <- matrix(0, n_sampling * n_chain, hnc)
    for(j in 1:hnc){
      sampling_data_mp_1p <- sampling_data_deg_1p <- matrix(0, n_sampling, n_chain)
      for(k in 1:n_chain){
        ind <- cand_col_index_per_chain[index_access_1[k], j]
        row_index <- (n_sampling * (k - 1) + 1):(n_sampling * k)
        sampling_data_mp_1p[, k] <- sampling_data[row_index, ind_mp[ind]]
        sampling_data_deg_1p[, k] <- sampling_data[row_index, ind_deg[ind]]
      }
      rhat_mp_all[i, j] <- rhat(sampling_data_mp_1p)
      rhat_deg_all[i, j] <- rhat(sampling_data_deg_1p)
    }
  }

  # Search index of the minimum rhat value
  index_min <- which.min(apply(rhat_mp_all, 1, sum) + apply(rhat_deg_all, 1, sum))
  index_access_best <- index_access_all[index_min, ]

  # Extract sampling data with the minimum rhat value
  sampling_data_mp <- sampling_data_deg <- matrix(0, n_sampling * n_chain, hnc)
  for(i in 1:n_chain){
    col_index_chain1 <- cand_col_index_per_chain[index_access_best[i], ]
    row_index <- (n_sampling * (i - 1) + 1):(n_sampling * i)
    sampling_data_mp[row_index, ] <- sampling_data[row_index, ind_mp[col_index_chain1]]
    sampling_data_deg[row_index, ] <- sampling_data[row_index, ind_deg[col_index_chain1]]
  }

  # Calculate the median mixture proportions and the degradation parameters
  mp_final <- apply(sampling_data_mp, 2, median)
  deg_final <- apply(sampling_data_deg, 2, median)
  index_order <- order(mp_final, decreasing = TRUE)
  mp_final <- mp_final[index_order]
  deg_final <- deg_final[index_order]

  return(list("mp_final" = mp_final,
              "deg_final" = deg_final))
}

#' determine_mpdeg
#'
#' @description The function to determine the mixture proportions and the degradation parameters
#' @param sampling_data A matrix of the sampling data obtained by the function 'as_draws_matrix'
#' @param n_sampling The number of sampling for MCMC
#' @param n_chain The number of chains for MCMC
#' @param upper_rhat The upper limit of the rhat value
determine_mpdeg <- function(sampling_data, n_sampling, n_chain, upper_rhat = 1.1){
  cols <- colnames(sampling_data)
  ind_mp <- grep("mp", cols)
  ind_deg <- grep("deg", cols)
  hnc <- length(ind_mp)

  # Calculate the medians of the mixture proportions in each chain
  median_mp <- matrix(0, hnc, n_chain)
  for(i in 1:n_chain){
    sampling_data_chain1 <- sampling_data[(n_sampling * (i - 1) + 1):(n_sampling * i), ]
    for(j in 1:hnc){
      median_mp[j, i] <- median(sampling_data_chain1[, ind_mp[j]])
    }
  }

  # Order sampling data according to the mixture proportions
  order_mp <- apply(median_mp, 2, order, decreasing = TRUE)
  rhat_mp <- rhat_deg <- rep(0, hnc)
  sampling_data_mp_list <- sampling_data_deg_list <- list()
  for(i in 1:hnc){
    sampling_data_mp_1p <- sampling_data_deg_1p <- matrix(0, n_sampling, n_chain)
    for(j in 1:n_chain){
      row_index <- (n_sampling * (j - 1) + 1):(n_sampling * j)
      sampling_data_mp_1p[, j] <- sampling_data[row_index, ind_mp[order_mp[i, j]]]
      sampling_data_deg_1p[, j] <- sampling_data[row_index, ind_deg[order_mp[i, j]]]
    }
    sampling_data_mp_list[[i]] <- sampling_data_mp_1p
    sampling_data_deg_list[[i]] <- sampling_data_deg_1p

    rhat_mp[i] <- rhat(sampling_data_mp_1p)
    rhat_deg[i] <- rhat(sampling_data_deg_1p)
  }

  if(all(rhat_mp < upper_rhat) && all(rhat_deg < upper_rhat)){
    mp_final <- sapply(sampling_data_mp_list, median)
    deg_final <- sapply(sampling_data_deg_list, median)
  }else{
    tmp <- search_min_rhat_mpdeg(sampling_data, n_sampling, n_chain, ind_mp, ind_deg, hnc)
    mp_final <- tmp$mp_final
    deg_final <- tmp$deg_final
  }

  return(list("mp_final" = mp_final,
              "deg_final" = deg_final))
}
