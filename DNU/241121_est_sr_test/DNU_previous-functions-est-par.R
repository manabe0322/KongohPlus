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
  # parent_al: parent allele data for stutter ratios
  # parent_ph: peak height data for parent alleles
  # total_al_prod: total allelic product data for stutter ratios
  bsr <- fsr <- dsr <- m2sr <- parent_al <- parent_ph <- total_al_prod <- rep(0, 2)

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
      parent_al[1] <- al_1
      parent_ph[1] <- h_al <- heights[peaks == al_1]
      total_al_prod[1] <- sum(heights)

      # Back stutter
      if(sr_consider["BS"]){
        bs <- calc_bs_one(al_1)
        if(is.element(bs, peaks)){
          bsr[1] <- heights[peaks == bs] / h_al
        }
      }

      # Forward stutter
      if(sr_consider["FS"]){
        fs <- calc_fs_one(al_1)
        if(is.element(fs, peaks)){
          fsr[1] <- heights[peaks == fs] / h_al
        }
      }

      # Double-back stutter
      if(sr_consider["DS"]){
        ds <- calc_ds_one(al_1)
        if(is.element(ds, peaks)){
          dsr[1] <- heights[peaks == ds] / h_al
        }
      }

      # Minus 2-nt stutter
      if(sr_consider["M2S"]){
        m2s <- calc_m2s_one(al_1)
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
          parent_al[i] <- al_1
          h_sum <- h_al <- heights[peaks == al_1]
          parent_ph[i] <- h_al

          # Candidates of shared peaks
          al_2 <- gt[3 - i]
          cand_share <- al_2
          if(sr_consider["BS"]){
            cand_share <- c(cand_share, calc_bs_one(al_2))
          }
          if(sr_consider["FS"]){
            cand_share <- c(cand_share, calc_fs_one(al_2))
          }
          if(sr_consider["DS"]){
            cand_share <- c(cand_share, calc_ds_one(al_2))
          }
          if(sr_consider["M2S"]){
            cand_share <- c(cand_share, calc_m2s_one(al_2))
          }
          cand_share <- round(sort(unique(cand_share)), 1)

          # Back stutter
          if(sr_consider["BS"]){
            bs <- calc_bs_one(al_1)
            if(is.element(bs, peaks) && !is.element(bs, cand_share)){
              h_bs <- heights[peaks == bs]
              bsr[i] <- h_bs / h_al
              h_sum <- h_sum + h_bs
            }
          }

          # Forward stutter
          if(sr_consider["FS"]){
            fs <- calc_fs_one(al_1)
            if(is.element(fs, peaks) && !is.element(fs, cand_share)){
              h_fs <- heights[peaks == fs]
              fsr[i] <- h_fs / h_al
              h_sum <- h_sum + h_fs
            }
          }

          # Double-back stutter
          if(sr_consider["DS"]){
            ds <- calc_ds_one(al_1)
            if(is.element(ds, peaks) && !is.element(ds, cand_share)){
              h_ds <- heights[peaks == ds]
              dsr[i] <- h_ds / h_al
              h_sum <- h_sum + h_ds
            }
          }

          # Minus 2-nt stutter
          if(sr_consider["M2S"]){
            m2s <- calc_m2s_one(al_1)
            if(is.element(m2s, peaks) && !is.element(m2s, cand_share)){
              h_m2s <- heights[peaks == m2s]
              m2sr[i] <- h_m2s / h_al
              h_sum <- h_sum + h_m2s
            }
          }
          total_al_prod[i] <- h_sum
        }
      }
    }
  }
  return(list(bsr, fsr, dsr, m2sr, parent_al, parent_ph, total_al_prod))
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

  bsr_data <- fsr_data <- dsr_data <- m2sr_data <- parent_al_data <- parent_ph_data <- total_al_prod_data <- matrix(0, 2 * n_sample, n_mk)

  # Define row and column names
  tmp <- sapply(names(height_list), rep, 2)
  tmp[1, ] <- paste0(tmp[1, ], "_data1")
  tmp[2, ] <- paste0(tmp[2, ], "_data2")
  rownames(bsr_data) <- rownames(fsr_data) <- rownames(dsr_data) <- rownames(m2sr_data) <- rownames(parent_al_data) <- rownames(parent_ph_data) <- rownames(total_al_prod_data) <- as.vector(tmp)
  colnames(bsr_data) <- colnames(fsr_data) <- colnames(dsr_data) <- colnames(m2sr_data) <- colnames(parent_al_data) <- colnames(parent_ph_data) <- colnames(total_al_prod_data) <- markers

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
        parent_al_data[c(2 * i - 1, 2 * i), j] <- sr_one_mk[[5]]
        parent_ph_data[c(2 * i - 1, 2 * i), j] <- sr_one_mk[[6]]
        total_al_prod_data[c(2 * i - 1, 2 * i), j] <- sr_one_mk[[7]]
      }
    }
  }

  return(list("bsr_data" = bsr_data,
              "fsr_data" = fsr_data,
              "dsr_data" = dsr_data,
              "m2sr_data" = m2sr_data,
              "parent_al_data" = parent_al_data,
              "parent_ph_data" = parent_ph_data,
              "total_al_prod_data" = total_al_prod_data))
}

#' analyze_seq
#'
#' @description The function to analyze sequence data
#' @param dt_seq The data.table of sequence data
#' @param markers_auto The vector of autosomal markers
analyze_seq <- function(dt_seq, markers_auto){
  lus_auto <- list()
  reps_auto <- list()
  seq_al_auto <- list()
  seq_count_auto <- list()

  for(i in 1:length(markers_auto)){
    # Extract data of one marker
    dt_seq_one_mk <- dt_seq[Marker == markers_auto[i]]
    alleles <- as.numeric(dt_seq_one_mk[, Allele])
    al_names <- sort(unique(alleles))
    counts <- dt_seq_one_mk[, Count]
    sequences <- dt_seq_one_mk[, RepeatRegion]

    # Investigate the number of repeats in each motif
    reps_one_mk <- lapply(sequences, pick_rep)

    # Calculate weighted LUS values
    lus_all_reps <- sapply(reps_one_mk, max)
    lus_one_mk <- rep(0, length(al_names))
    for(j in 1:length(al_names)){
      index_al <- which(alleles == al_names[j])
      lus_one_mk[j] <- calc_weighted_lus(lus_all_reps[index_al], counts[index_al])
    }
    names(lus_one_mk) <- al_names

    # Assign data of one marker
    lus_auto[[i]] <- lus_one_mk
    reps_auto[[i]] <- reps_one_mk
    seq_al_auto[[i]] <- alleles
    seq_count_auto[[i]] <- counts
  }

  names(lus_auto) <- markers_auto
  names(reps_auto) <- markers_auto
  names(seq_al_auto) <- markers_auto
  names(seq_count_auto) <- markers_auto

  return(list("lus_auto" = lus_auto,
              "reps_auto" = reps_auto,
              "seq_al_auto" = seq_al_auto,
              "seq_count_auto" = seq_count_auto))
}

#' make_sr_ca_one_mk
#'
#' @description The function to make data of corrected allele numbers in one marker
#' @param min_us Minimum uninterrupted stretch
#' @param reps_one_mk Repeat motifs of one marker
#' @param seq_al_one_mk Alleles for sequence data in one marker
#' @param seq_count_one_mk Counts for sequence data in one marker
#' @param sr_al_one_mk A vector of the parent alleles of each stutter ratio
make_sr_ca_one_mk <- function(min_us, reps_one_mk, seq_al_one_mk, seq_count_one_mk, sr_al_one_mk){
  sr_ca_one_mk <- rep(0, length(sr_al_one_mk))
  al_names <- sort(unique(seq_al_one_mk))
  for(i in 1:length(al_names)){
    index_seq_al <- which(seq_al_one_mk == al_names[i])
    sr_ca_one_mk[which(sr_al_one_mk == al_names[i])] <- calc_weighted_ca(min_us, reps_one_mk[index_seq_al], seq_count_one_mk[index_seq_al])
  }
  return(sr_ca_one_mk)
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
  sr_one_mk <- as.numeric(sr_one_mk[index_non_zero])
  sr_al_one_mk <- as.numeric(sr_al_one_mk[index_non_zero])
  sr_h_one_mk <- as.numeric(sr_h_one_mk[index_non_zero])

  # Assign the LUS value to each allele
  al_names <- as.numeric(names(lus_one_mk))
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
