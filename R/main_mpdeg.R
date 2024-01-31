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
  sampling_data_mp <- sampling_data_deg <- matrix(0, n_sampling * n_chain, hnc)
  for(i in 1:n_chain){
    row_index <- (n_sampling * (i - 1) + 1):(n_sampling * i)
    sampling_data_mp[row_index, ] <- sampling_data[row_index, ind_mp[order_mp[, i]]]
    sampling_data_deg[row_index, ] <- sampling_data[row_index, ind_deg[order_mp[, i]]]
  }

  # Calculate rhat values
  rhat_mp <- rhat_deg <- rep(10, hnc)
  for(i in 1:hnc){
    rhat_mp[i] <- rhat(sampling_data_mp[, i])
    rhat_deg[i] <- rhat(sampling_data_deg[, i])
  }

  if(all(rhat_mp < upper_rhat) && all(rhat_deg < upper_rhat)){
    mp_final <- apply(sampling_data_mp, 2, median)
    deg_final <- apply(sampling_data_deg, 2, median)
  }else{
    # search_min_rhat_mpdeg
    mp_final <- NULL
    deg_final <- NULL
  }

  return(list("mp_final" = mp_final,
              "deg_final" = deg_final))
}
