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

#' pick_rep
#'
#' @description The function to investigate the number of repeats in each motif
#' @param seq A character of the sequence of the repeat region
pick_rep <- function(seq){
  seq_sep <- strsplit(seq, " ")[[1]]
  reps <- integer(0)
  for(i in 1:length(seq_sep)){
    motif <- seq_sep[i]
    pos_bracket <- regexpr("]", motif)[1]
    if(pos_bracket != -1){
      reps <- c(reps, as.numeric(strsplit(motif, ']')[[1]][2]))
    }
  }
  return(reps)
}

#' calc_weighted_lus
#'
#' @description The function to calculate the weighted LUS of an allele
#' @param lus_one_al A vector of the LUS values of an allele
#' @param count_one_al A vector of the number of counts of an allele
calc_weighted_lus <- function(lus_one_al, count_one_al){
  return(sum(lus_one_al * count_one_al / sum(count_one_al)))
}

#' calc_weighted_ca
#'
#' @description The function to calculate the weighted corrected allele number
#' @param min_us Minimum uninterrupted stretch
#' @param reps_one_al Repeat motifs of one allele
#' @param seq_count_one_al Counts for sequence data in one allele
calc_weighted_ca <- function(min_us, reps_one_al, seq_count_one_al){
  # Calculate corrected allele numbers in each sequence
  n_seq_al <- length(reps_one_al)
  ca_one_al <- rep(0, n_seq_al)
  for(i in 1:n_seq_al){
    tmp <- reps_one_al[[i]] - min_us
    ca_one_al[i] <- sum(tmp[tmp > 0])
  }

  # Calculate the weighted corrected allele number
  weighted_ca <- sum(ca_one_al * seq_count_one_al / sum(seq_count_one_al))
  return(weighted_ca)
}

#' analyze_seq
#'
#' @description The function to analyze sequence data
#' @param dt_seq_init The initial data.table of sequence data
#' @param markers_auto The vector of autosomal markers
analyze_seq <- function(dt_seq_init, markers_auto){
  lus_auto <- list()
  reps_auto <- list()
  seq_al_auto <- list()
  seq_count_auto <- list()

  for(i in 1:length(markers_auto)){
    # Extract data of one marker
    dt_seq_one_mk <- dt_seq_init[Marker == markers_auto[i]]
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

#' estimate_average_sr
#'
#' @description The function to estimate average stutter ratios
estimate_average_sr <- function(dt_exp, dt_gt, dt_seq, dt_kit, dt_setting_est_sr, kit_name,
                                model_sr_gen, model_sr_uni,
                                warmup, sampling, chains, pu_remove, pu_dist, at_sex,
                                n_parallel_chains = detectCores()){
  # Extract information on the kit
  markers_all <- unique(dt_kit$Marker)
  n_mk_all <- length(markers_all)
  dye_all_mk <- extract_dye(dt_kit, markers_all)

  # Extract information on setting
  markers_auto <- dt_setting_est_sr$Marker
  n_mk_auto <- length(markers_auto)
  at_auto <- dt_setting_est_sr$AT
  names(at_auto) <- markers_auto
  meth_model_bsr_auto <- dt_setting_est_sr$Model_BSR
  meth_model_fsr_auto <- dt_setting_est_sr$Model_FSR
  meth_model_dsr_auto <- dt_setting_est_sr$Model_DSR
  meth_model_m2sr_auto <- dt_setting_est_sr$Model_M2SR

  # Make an object for ATs of all markers
  at_all_mk <- rep(NA, n_mk_all)
  names(at_all_mk) <- markers_all
  at_all_mk[match(markers_auto, markers_all)] <- at_auto
  at_all_mk[is.na(at_all_mk)] <- at_sex

  # Information on considering stutter ratios
  sr_consider_auto <- matrix(TRUE, n_mk_auto, 4)
  rownames(sr_consider_auto) <- markers_auto
  colnames(sr_consider_auto) <- c("BS", "FS", "DS", "M2S")
  sr_consider_auto[, 1] <- meth_model_bsr_auto != "Not consider"
  sr_consider_auto[, 2] <- meth_model_fsr_auto != "Not consider"
  sr_consider_auto[, 3] <- meth_model_dsr_auto != "Not consider"
  sr_consider_auto[, 4] <- meth_model_m2sr_auto != "Not consider"

  # Arrange experimental data
  tmp <- arrange_exp_data(dt_exp, dt_gt,
                          markers_auto, sr_consider_auto,
                          markers_all, at_all_mk, dye_all_mk,
                          pu_remove, pu_dist)
  peak_list <- tmp$peak_list
  size_list <- tmp$size_list
  height_list <- tmp$height_list
  gt_list <- tmp$gt_list

  # Prepare data for the stutter ratios
  tmp <- get_sr_data(peak_list, height_list, gt_list, sr_consider_auto)
  bsr_data <- tmp$bsr_data
  fsr_data <- tmp$fsr_data
  dsr_data <- tmp$dsr_data
  m2sr_data <- tmp$m2sr_data
  parent_al_data <- tmp$parent_al_data
  parent_ph_data <- tmp$parent_ph_data
  total_al_prod_data <- tmp$total_al_prod_data

  # Prepare data related to the sequence data
  tmp <- analyze_seq(dt_seq, markers_auto)
  lus_auto <- tmp$lus_auto
  reps_auto <- tmp$reps_auto
  seq_al_auto <- tmp$seq_al_auto
  seq_count_auto <- tmp$seq_count_auto
}
