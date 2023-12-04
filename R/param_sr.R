#' calc_sr
#'
#' @description The function to calculate stutter ratios
#' @param peaks A character vector of peaks
#' @param heights A numeric vector of heights
#' @param gt A character vector of the genotype in one marker
calc_sr <- function(peaks, heights, gt, sr_consider){

  # bsr: back stutter ratio
  # fsr: forward stutter ratio
  # dsr: double-back stutter ratio
  # m2sr: minus 2-nt stutter ratio
  # allele_sr: allele data for stutter ratios
  # height_sr: height data for stutter ratios (i.e., total allelic product)
  bsr <- fsr <- dsr <- m2sr <- allele_sr <- height_sr <- rep(0, 2)

  peaks <- as.numeric(peaks)

  gt <- as.numeric(gt)
  gt <- sort(unique(gt[!is.na(gt)]))

  ##############
  # Homozygote #
  ##############

  if(length(gt) == 1){
    al_1 <- gt

    # No allelic drop-out
    if(is.element(al_1, peaks)){
      allele_sr[1] <- al_1
      height_sr[1] <- sum(heights)
      h_al <- heights[peaks == al_1]

      # Back stutter
      if(sr_consider["BS"]){
        bs <- round(al_1 - 1, 1)
        if(is.element(bs, peaks)){
          bsr[1] <- heights[peaks == bs] / h_al
        }
      }

      # Forward stutter
      if(sr_consider["FS"]){
        fs <- round(al_1 + 1, 1)
        if(is.element(fs, peaks)){
          fsr[1] <- heights[peaks == fs] / h_al
        }
      }

      # Double-back stutter
      if(sr_consider["DS"]){
        ds <- round(al_1 - 2, 1)
        if(is.element(ds, peaks)){
          dsr[1] <- heights[peaks == ds] / h_al
        }
      }

      # Minus 2-nt stutter
      if(sr_consider["M2S"]){
        if(al_1 %% 1 < 0.15){
          m2s <- round(al_1 - 0.8, 1)  # Allele x or x.1
        }else{
          m2s <- round(al_1 - 0.2, 1)  # Allele x.2, x.3, ...
        }
        if(is.element(m2s, peaks)){
          m2sr[1] <- heights[peaks == m2s] / h_al
        }
      }
    }

    ################
    # Heterozygote #
    ################

  }else{

    # Other than one repeat difference
    if(abs(round(gt[1] - gt[2], 1)) != 1){

      for(i in 1:2){
        al_1 <- gt[i]

        # No allelic drop-out
        if(is.element(al_1, peaks)){
          allele_sr[i] <- al_1
          h_sum <- h_al <- heights[peaks == al_1]

          # Candidates of shared peaks
          al_2 <- gt[3 - i]
          cand_share <- al_2
          if(sr_consider["BS"]){
            cand_share <- c(cand_share, al_2 - 1)
          }
          if(sr_consider["FS"]){
            cand_share <- c(cand_share, al_2 + 1)
          }
          if(sr_consider["DS"]){
            cand_share <- c(cand_share, al_2 - 2)
          }
          if(sr_consider["M2S"]){
            if(al_2 %% 1 < 0.15){
              cand_share <- c(cand_share, al_2 - 0.8)  # Allele x or x.1
            }else{
              cand_share <- c(cand_share, al_2 - 0.2)  # Allele x.2, x.3, ...
            }
          }
          cand_share <- round(sort(unique(cand_share)), 1)

          # Back stutter
          if(sr_consider["BS"]){
            bs <- round(al_1 - 1, 1)
            if(is.element(bs, peaks) && !is.element(bs, cand_share)){
              h_bs <- heights[peaks == bs]
              bsr[i] <- h_bs / h_al
              h_sum <- h_sum + h_bs
            }
          }

          # Forward stutter
          if(sr_consider["FS"]){
            fs <- round(al_1 + 1, 1)
            if(is.element(fs, peaks) && !is.element(fs, cand_share)){
              h_fs <- heights[peaks == fs]
              fsr[i] <- h_fs / h_al
              h_sum <- h_sum + h_fs
            }
          }

          # Double-back stutter
          if(sr_consider["DS"]){
            ds <- round(al_1 - 2, 1)
            if(is.element(ds, peaks) && !is.element(ds, cand_share)){
              h_ds <- heights[peaks == ds]
              dsr[i] <- h_ds / h_al
              h_sum <- h_sum + h_ds
            }
          }

          # Minus 2-nt stutter
          if(sr_consider["M2S"]){
            if(al_1 %% 1 < 0.15){
              m2s <- round(al_1 - 0.8, 1)  # Allele x or x.1
            }else{
              m2s <- round(al_1 - 0.2, 1)  # Allele x.2, x.3, ...
            }
            if(is.element(m2s, peaks) && !is.element(m2s, cand_share)){
              h_m2s <- heights[peaks == m2s]
              m2sr[i] <- h_m2s / h_al
              h_sum <- h_sum + h_m2s
            }
          }
          height_sr[i] <- h_sum
        }
      }
    }
  }
  return(list(bsr, fsr, dsr, m2sr, allele_sr, height_sr))
}

#' get_sr_data
#'
#' @description The function to get data of stutter ratios
#' @param peak_list A list of peaks
#' @param height_list A list of heights
#' @param sr_consider_all_mk Information on consideration of back, forward, double-back, and minus 2-nt stutters in all markers
get_sr_data <- function(peak_list, height_list, gt_list, sr_consider_all_mk){
  n_sample <- length(peak_list)
  markers <- names(peak_list[[1]])
  n_mk <- length(markers)

  bsr_data <- fsr_data <- dsr_data <- m2sr_data <- allele_sr_data <- height_sr_data <- matrix(0, 2 * n_sample, n_mk)

  # Define row and column names
  tmp <- sapply(names(height_list), rep, 2)
  tmp[1, ] <- paste0(tmp[1, ], "_data1")
  tmp[2, ] <- paste0(tmp[2, ], "_data2")
  rownames(bsr_data) <- rownames(fsr_data) <- rownames(dsr_data) <- rownames(m2sr_data) <- rownames(allele_sr_data) <- rownames(height_sr_data) <- as.vector(tmp)
  colnames(bsr_data) <- colnames(fsr_data) <- colnames(dsr_data) <- colnames(m2sr_data) <- colnames(allele_sr_data) <- colnames(height_sr_data) <- markers

  for(i in 1:n_sample){
    peak_list_one <- peak_list[[i]]
    height_list_one <- height_list[[i]]
    gt_list_one <- gt_list[[i]]

    for(j in 1:n_mk){
      peaks <- peak_list_one[[j]]
      heights <- height_list_one[[j]]
      gt <- gt_list_one[j, ]

      if(length(peaks) > 0){
        sr_consider <- sr_consider_all_mk[j, ]

        sr_one_mk <- calc_sr(peaks, heights, gt, sr_consider)
        bsr_data[c(2 * i - 1, 2 * i), j] <- sr_one_mk[[1]]
        fsr_data[c(2 * i - 1, 2 * i), j] <- sr_one_mk[[2]]
        dsr_data[c(2 * i - 1, 2 * i), j] <- sr_one_mk[[3]]
        m2sr_data[c(2 * i - 1, 2 * i), j] <- sr_one_mk[[4]]
        allele_sr_data[c(2 * i - 1, 2 * i), j] <- sr_one_mk[[5]]
        height_sr_data[c(2 * i - 1, 2 * i), j] <- sr_one_mk[[6]]
      }
    }
  }

  sr_list <- list()
  sr_list[[1]] <- bsr_data
  sr_list[[2]] <- fsr_data
  sr_list[[3]] <- dsr_data
  sr_list[[4]] <- m2sr_data
  sr_list[[5]] <- allele_sr_data
  sr_list[[6]] <- height_sr_data
  names(sr_list) <- c("BSR", "FSR", "DSR", "M2SR", "Allele_SR", "Height_SR")
  return(sr_list)
}

#' extract_sr_use
#'
#' @description The function to extract SR data used
#' @param sr_one_mk A vector of the stutter ratios in one marker
#' @param sr_al_one_mk A vector of the parent alleles of each stutter ratio
#' @param sr_h_one_mk A vector of the total allelic products of each stutter ratio
#' @param lus_one_mk A vector of the LUS of each allele in one marker
extract_sr_use <- function(sr_one_mk, sr_al_one_mk, sr_h_one_mk, lus_one_mk){
  # Remove data without stutter ratio
  index_non_zero <- which(sr_one_mk != 0)
  sr_one_mk <- sr_one_mk[index_non_zero]
  sr_al_one_mk <- sr_al_one_mk[index_non_zero]
  sr_h_one_mk <- sr_h_one_mk[index_non_zero]

  # Assign the LUS value to each allele
  al_names <- names(lus_one_mk)
  sr_lus_one_mk <- rep(-1, length(sr_one_mk))
  for(i in 1:length(lus_one_mk)){
    index_al <- which(sr_al_one_mk == al_names[i])
    if(length(index_al) > 0){
      sr_lus_one_mk[index_al] <- lus_one_mk[i]
    }
  }

  # Remove data without sequence information
  index_no_seq <- which(sr_lus_one_mk == -1)
  if(length(index_no_seq) > 0){
    sr_one_mk <- sr_one_mk[- index_no_seq]
    sr_al_one_mk <- sr_al_one_mk[- index_no_seq]
    sr_h_one_mk <- sr_h_one_mk[- index_no_seq]
    sr_lus_one_mk <- sr_lus_one_mk[- index_no_seq]
  }

  # Make an output object
  sr_use_list <- list(sr_one_mk, sr_al_one_mk, sr_h_one_mk, sr_lus_one_mk)
  names(sr_use_list) <- c("stutter ratio", "parent allele", "total allelic product", "LUS of parent allele")
  return(sr_use_list)
}
