divide_mp <- function(mp_select, n_divide){
  return((1 - sum(mp_select)) / n_divide)
}

judge_mp_greater0 <- function(mp_cand_1){
  return(all(mp_cand_1 > 0))
}

#' set_mp
#'
#' @description The function to set mixture proportions
#' @param hnc The hypothesized number of contributors
#' @param mp_cand_c1 A vector of mixture proportions per contributor
set_mp <- function(hnc, mp_cand_c1){
  if(hnc == 1){
    return(matrix(1, 1, 1))
  }else{
    mp_cand_c1 <- round(mp_cand_c1, 5)
    n_mp_cand_c1 <- length(mp_cand_c1)

    n_mp_per_choose <- rep(0, hnc - 1)
    for(i in 1:(hnc - 1)){
      n_mp_per_choose[i] <- choose(n_mp_cand_c1 + i - 1, i)
    }

    mp_cand <- matrix(0, sum(n_mp_per_choose) + 1, hnc)
    row_start <- 1
    for(i in 1:(hnc - 1)){
      row_indices <- row_start:(row_start + n_mp_per_choose[i] - 1)
      mp_cand[row_indices, 1:i] <- mp_select <- combinations(n = n_mp_cand_c1, r = i, v = mp_cand_c1, repeats.allowed = TRUE)
      mp_divide <- apply(mp_select, 1, divide_mp, hnc - i)
      mp_divide <- rep(mp_divide, hnc - i)
      mp_cand[row_indices, (i + 1):hnc] <- mp_divide
      row_start <- row_start + n_mp_per_choose[i]
    }
    mp_cand[nrow(mp_cand), ] <- rep(1 / hnc, hnc)
    mp_cand <- round(mp_cand, 5)

    # Cut mixture proportions which includes zero
    mp_cand <- mp_cand[apply(mp_cand, 1, judge_mp_greater0), , drop = FALSE]

    # Cut overlapped mixture proportions
    mp_cand <- t(apply(mp_cand, 1, sort, decreasing = TRUE))
    mp_cut_count <- 1
    repeat{
      index_mp_cut <- setdiff(which(apply(mp_cand, 1, setequal, mp_cand[mp_cut_count, ])), mp_cut_count)
      if(length(index_mp_cut) != 0){
        mp_cand <- mp_cand[-index_mp_cut, , drop = FALSE]
      }
      mp_cut_count <- mp_cut_count + 1
      if(nrow(mp_cand) < mp_cut_count){
        break
      }
    }

    return(mp_cand)
  }
}

library(gtools)
mp_cand_c1 <- c(0.0025, 0.005, 0.0075, 0.01, 0.025, 0.05, 0.075, 0.1, 0.2, 0.3, 0.4, 0.5)
mp_cand_2p_pre <- set_mp(2, mp_cand_c1)
mp_cand_3p_pre <- set_mp(3, mp_cand_c1)
mp_cand_4p_pre <- set_mp(4, mp_cand_c1)

mp_cand_2p <- mp_cand_2p_pre * 0
hnc <- 2
for(i in 1:hnc){
  mp_cand_2p[, i] <- mp_cand_2p_pre[, hnc - i + 1]
}

mp_cand_3p <- mp_cand_3p_pre * 0
hnc <- 3
for(i in 1:hnc){
  mp_cand_3p[, i] <- mp_cand_3p_pre[, hnc - i + 1]
}

mp_cand_4p <- mp_cand_4p_pre * 0
hnc <- 4
for(i in 1:hnc){
  mp_cand_4p[, i] <- mp_cand_4p_pre[, hnc - i + 1]
}

