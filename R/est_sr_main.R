#' extract_sr_data
#'
#' @description The function to extract stutter ratios
#' @param peaks A character vector of peaks
#' @param heights A numeric vector of heights
#' @param gt A character vector of the genotype in one marker
#' @param sr_consider Information on consideration of back, forward, double-back, and minus 2-nt stutters in one marker
extract_sr_data <- function(peaks, heights, gt, sr_consider){
  peaks <- as.numeric(peaks)
  n_peak <- length(peaks)
  gt <- as.numeric(gt)
  gt <- sort(unique(gt[!is.na(gt)]))
  n_al <- length(gt)

  parent_al <- stutter <- rep("", n_peak)
  stutter_ph <- parent_ph <- total_al_prod <- stutter_ratio <- rep(NA, n_peak)
  bool_al <- bool_bs <- bool_fs <- bool_ds <- bool_m2s <- bool_shared_allele <- bool_shared_stutter <- rep(FALSE, n_peak)

  # At least one peak
  if(n_peak > 0){
    bool_al <- is.element(peaks, gt)
    if(sr_consider["BS"]){
      bs <- calc_bs(gt)
      bool_bs <- is.element(peaks, bs)
    }
    if(sr_consider["FS"]){
      fs <- calc_fs(gt)
      bool_fs <- is.element(peaks, fs)
    }
    if(sr_consider["DS"]){
      ds <- calc_ds(gt)
      bool_ds <- is.element(peaks, ds)
    }
    if(sr_consider["M2S"]){
      m2s <- calc_m2s(gt)
      bool_m2s <- is.element(peaks, m2s)
    }

    bool_shared_allele_base <- bool_al & bool_bs + bool_fs + bool_ds + bool_m2s > 0
    bool_shared_stutter <- bool_bs + bool_fs + bool_ds + bool_m2s > 1

    for(i in 1:n_peak){
      if(!bool_al[i]){
        peak <- peaks[i]
        if(!bool_shared_stutter[i]){
          if(bool_bs[i]){
            al <- calc_parent_al_bs_one(peak)
          }else if(bool_fs[i]){
            al <- calc_parent_al_fs_one(peak)
          }else if(bool_ds[i]){
            al <- calc_parent_al_ds_one(peak)
          }else if(bool_m2s[i]){
            al <- calc_parent_al_m2s_one(peak)
          }else{
            al <- NA
          }

          index_al <- which(is.element(peaks, al))

          parent_al[i] <- as.character(al)
          stutter[i] <- as.character(peak)
          stutter_ph[i] <- heights[i]

          # The parent allele is above AT
          if(length(index_al) == 1){
            parent_ph[i] <- heights[index_al]
            total_al_prod[i] <- stutter_ph[i] + parent_ph[i]
            stutter_ratio[i] <- stutter_ph[i] / parent_ph[i]
            bool_shared_allele[i] <- bool_shared_allele_base[index_al]
          }
        }

        # Shared stutters
        else{
          parent_al[i] <- "Multiple"
          stutter[i] <- as.character(peak)
          stutter_ph[i] <- heights[i]
        }
      }
    }

    parent_al <- parent_al[!bool_al]
    stutter <- stutter[!bool_al]
    stutter_ph <- stutter_ph[!bool_al]
    parent_ph <- parent_ph[!bool_al]
    total_al_prod <- total_al_prod[!bool_al]
    stutter_ratio <- stutter_ratio[!bool_al]
    bool_bs <- bool_bs[!bool_al]
    bool_fs <- bool_fs[!bool_al]
    bool_ds <- bool_ds[!bool_al]
    bool_m2s <- bool_m2s[!bool_al]
    bool_shared_allele <- bool_shared_allele[!bool_al]
  }

  return(list("parent_al" = parent_al,
              "stutter" = stutter,
              "stutter_ph" = stutter_ph,
              "parent_ph" = parent_ph,
              "total_al_prod" = total_al_prod,
              "stutter_ratio" = stutter_ratio,
              "bool_bs" = bool_bs,
              "bool_fs" = bool_fs,
              "bool_ds" = bool_ds,
              "bool_m2s" = bool_m2s,
              "bool_shared_allele" = bool_shared_allele))
}

#' create_dt_sr_data
#'
#' @description The function to create a data.table of the stutter ratios
#' @param peak_list A list of peaks
#' @param height_list A list of heights
#' @param gt_list A list of genotypes
#' @param sr_consider_auto Information on consideration of back, forward, double-back, and minus 2-nt stutters in autosomal markers
#' @param wlus_auto A list of weight LUS values in autosomal markers
create_dt_sr_data <- function(peak_list, height_list, gt_list, sr_consider_auto, wlus_auto){
  sample_names <- names(peak_list)
  n_sample <- length(sample_names)
  markers <- names(peak_list[[1]])
  n_mk <- length(markers)

  n_max_data <- 8 * n_sample * n_mk
  sample_name_all <- rep("", n_max_data)
  marker_all <- rep("", n_max_data)
  parent_al_all <- rep("", n_max_data)
  stutter_all <- rep("", n_max_data)
  stutter_ratio_all <- rep(NA, n_max_data)
  stutter_ph_all <- rep(NA, n_max_data)
  parent_ph_all <- rep(NA, n_max_data)
  total_al_prod_all <- rep(NA, n_max_data)
  bool_bs_all <- rep(FALSE, n_max_data)
  bool_fs_all <- rep(FALSE, n_max_data)
  bool_ds_all <- rep(FALSE, n_max_data)
  bool_m2s_all <- rep(FALSE, n_max_data)
  bool_shared_allele_all <- rep(FALSE, n_max_data)
  bool_seq_all <- rep(FALSE, n_max_data)

  current_index <- 1
  for(i in 1:n_sample){
    peak_list_one <- peak_list[[i]]
    height_list_one <- height_list[[i]]
    gt_list_one <- gt_list[[i]]

    for(j in 1:n_mk){
      peaks <- peak_list_one[[j]]
      if(length(peaks) > 0){
        heights <- height_list_one[[j]]
        gt <- gt_list_one[j, ]
        sr_consider <- sr_consider_auto[j, ]
        al_names_with_seq <- names(wlus_auto[[j]])

        tmp <- extract_sr_data(peaks, heights, gt, sr_consider)
        parent_al <- tmp$parent_al
        n_data <- length(parent_al)

        if(n_data > 0){
          index_assign <- current_index:(current_index + n_data - 1)
          sample_name_all[index_assign] <- sample_names[i]
          marker_all[index_assign] <- markers[j]
          parent_al_all[index_assign] <- parent_al
          stutter_all[index_assign] <- tmp$stutter
          stutter_ratio_all[index_assign] <- tmp$stutter_ratio
          stutter_ph_all[index_assign] <- tmp$stutter_ph
          parent_ph_all[index_assign] <- tmp$parent_ph
          total_al_prod_all[index_assign] <- tmp$total_al_prod
          bool_bs_all[index_assign] <- tmp$bool_bs
          bool_fs_all[index_assign] <- tmp$bool_fs
          bool_ds_all[index_assign] <- tmp$bool_ds
          bool_m2s_all[index_assign] <- tmp$bool_m2s
          bool_shared_allele_all[index_assign] <- tmp$bool_shared_allele
          bool_seq_all[index_assign] <- is.element(parent_al, al_names_with_seq)

          current_index <- current_index + n_data
        }
      }
    }
  }

  sample_name_all <- sample_name_all[1:(current_index - 1)]
  marker_all <- marker_all[1:(current_index - 1)]
  parent_al_all <- parent_al_all[1:(current_index - 1)]
  stutter_all <- stutter_all[1:(current_index - 1)]
  stutter_ratio_all <- stutter_ratio_all[1:(current_index - 1)]
  stutter_ph_all <- stutter_ph_all[1:(current_index - 1)]
  parent_ph_all <- parent_ph_all[1:(current_index - 1)]
  total_al_prod_all <- total_al_prod_all[1:(current_index - 1)]
  bool_bs_all <- bool_bs_all[1:(current_index - 1)]
  bool_fs_all <- bool_fs_all[1:(current_index - 1)]
  bool_ds_all <- bool_ds_all[1:(current_index - 1)]
  bool_m2s_all <- bool_m2s_all[1:(current_index - 1)]
  bool_shared_allele_all <- bool_shared_allele_all[1:(current_index - 1)]
  bool_seq_all <- bool_seq_all[1:(current_index - 1)]

  dt_sr_data <- data.table(SampleName = sample_name_all,
                           Marker = marker_all,
                           Allele = parent_al_all,
                           Stutter = stutter_all,
                           SR = stutter_ratio_all,
                           PH_Stutter = stutter_ph_all,
                           PH_Allele = parent_ph_all,
                           TAP = total_al_prod_all,
                           Bool_BS = bool_bs_all,
                           Bool_FS = bool_fs_all,
                           Bool_DS = bool_ds_all,
                           Bool_M2S = bool_m2s_all,
                           Bool_SharedAllele = bool_shared_allele_all,
                           Bool_SeqData = bool_seq_all)

  return(dt_sr_data)
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
#' @param count_one_al Counts for sequence data in one allele
calc_weighted_ca <- function(min_us, reps_one_al, count_one_al){
  # Calculate corrected allele numbers in each sequence
  n_seq_al <- length(reps_one_al)
  ca_one_al <- rep(0, n_seq_al)
  for(i in 1:n_seq_al){
    tmp <- reps_one_al[[i]] - min_us
    ca_one_al[i] <- sum(tmp[tmp > 0])
  }

  # Calculate the weighted corrected allele number
  weighted_ca <- sum(ca_one_al * count_one_al / sum(count_one_al))
  return(weighted_ca)
}

#' correct_al_rep
#'
#' @description The function to correct allele repeat numbers
#' @param dt_seq The data.table of sequence data
#' @param markers_auto A vector of the autosomal markers
#' @param min_us_auto A vector of the minimum uninterrupted stretches of the autosomal markers
correct_al_rep <- function(dt_seq, markers_auto, min_us_auto){
  setindex(dt_seq, Marker)
  n_mk <- length(markers_auto)
  wlus_auto <- vector(mode = "list", length = n_mk)
  wca_auto <- vector(mode = "list", length = n_mk)

  for(i in 1:n_mk){
    # Extract data of one marker
    dt_seq_one_mk <- dt_seq[Marker == markers_auto[i]]
    alleles <- as.numeric(dt_seq_one_mk$Allele)
    al_names <- sort(unique(alleles))
    n_al <- length(al_names)
    counts <- dt_seq_one_mk$Count
    sequences <- dt_seq_one_mk$RepeatRegion
    min_us <- min_us_auto[i]

    # Investigate the number of repeats in each motif
    reps_one_mk <- lapply(sequences, pick_rep)

    # Calculate weighted LUS values
    lus_all_reps <- sapply(reps_one_mk, max)
    lus_one_mk <- rep(0, n_al)
    for(j in 1:n_al){
      index_al <- which(alleles == al_names[j])
      lus_one_mk[j] <- calc_weighted_lus(lus_all_reps[index_al], counts[index_al])
    }
    names(lus_one_mk) <- al_names

    # Calculate weighted CA values
    ca_one_mk <- rep(0, n_al)
    for(j in 1:n_al){
      index_al <- which(alleles == al_names[j])
      ca_one_mk[j] <- calc_weighted_ca(min_us, reps_one_mk[index_al], counts[index_al])
    }
    names(ca_one_mk) <- al_names

    # Assign data of one marker
    wlus_auto[[i]] <- lus_one_mk
    wca_auto[[i]] <- ca_one_mk
  }

  names(wlus_auto) <- markers_auto
  names(wca_auto) <- markers_auto

  return(list("wlus_auto" = wlus_auto,
              "wca_auto" = wca_auto))
}

#' run_mcmc_sr_gen
#'
#' @description The function to perform MCMC for the stutter ratio (locus-specific with slope)
#' @param model_sr_gen The stan model of the stutter ratio
#' @param data_list_sr The data list of the stutter ratio
#' @param warmup The number of warm-up iterations for each chain
#' @param sampling The number of sampling iterations for each chain
#' @param chains The number of chains
#' @param n_parallel_chains The number of parallel chains
run_mcmc_sr_gen <- function(model_sr_gen, data_list_sr, warmup, sampling, chains, n_parallel_chains){
  # Perform MCMC
  fit_sr <- model_sr_gen$sample(
    data_list_sr,
    iter_warmup = warmup,
    iter_sampling = sampling,
    chains = chains,
    parallel_chains = n_parallel_chains
  )

  # Summary data
  fit_sr_summary <- fit_sr$summary()
  medians <- fit_sr_summary$median
  variables <- fit_sr_summary$variable
  rhats <- fit_sr_summary$rhat
  median_lp <- medians[which(variables == "lp__")]
  median_slope <- medians[which(variables == "slope")]
  median_intercept <- medians[which(variables == "intercept")]
  median_sigma <- medians[which(variables == "sigma")]
  rhat_slope <- rhats[which(variables == "slope")]
  rhat_intercept <- rhats[which(variables == "intercept")]
  rhat_sigma <- rhats[which(variables == "sigma")]
  rhat_lp <- rhats[which(variables == "lp__")]

  # Trace plot
  fit_sr$draws() %>% mcmc_trace(pars = c("slope", "intercept", "sigma")) -> traceplot

  return(list("median_lp" = median_lp,
              "median_slope" = median_slope,
              "median_intercept" = median_intercept,
              "median_sigma" = median_sigma,
              "rhat_slope" = rhat_slope,
              "rhat_intercept" = rhat_intercept,
              "rhat_sigma" = rhat_sigma,
              "rhat_lp" = rhat_lp,
              "traceplot" = traceplot))
}

#' calc_average_sr
#'
#' @description The function to calculate the average stutter ratio
#' @param allele A vector of alleles or corrected alleles
#' @param slope The slope parameter
#' @param intercept The intercept parameter
calc_average_sr <- function(allele, slope, intercept){
  return(allele * slope + intercept)
}

#' run_mcmc_sr_uni
#'
#' @description The function to perform MCMC for the stutter ratio (uniform)
#' @param model_sr_uni The stan model of the stutter ratio
#' @param data_list_sr The data list of the stutter ratio
#' @param warmup The number of warm-up iterations for each chain
#' @param sampling The number of sampling iterations for each chain
#' @param chains The number of chains
#' @param n_parallel_chains The number of parallel chains
run_mcmc_sr_uni <- function(model_sr_uni, data_list_sr, warmup, sampling, chains, n_parallel_chains){
  # Perform MCMC
  fit_sr <- model_sr_uni$sample(
    data_list_sr,
    iter_warmup = warmup,
    iter_sampling = sampling,
    chains = chains,
    parallel_chains = n_parallel_chains
  )

  # Summary data
  fit_sr_summary <- fit_sr$summary()
  medians <- fit_sr_summary$median
  variables <- fit_sr_summary$variable
  rhats <- fit_sr_summary$rhat
  median_lp <- medians[which(variables == "lp__")]
  median_myu <- medians[which(variables == "myu")]
  median_sigma <- medians[which(variables == "sigma")]
  rhat_myu <- rhats[which(variables == "myu")]
  rhat_sigma <- rhats[which(variables == "sigma")]
  rhat_lp <- rhats[which(variables == "lp__")]

  # Trace plot
  fit_sr$draws() %>% mcmc_trace(pars = c("myu", "sigma")) -> traceplot

  return(list("median_lp" = median_lp,
              "median_myu" = median_myu,
              "median_sigma" = median_sigma,
              "rhat_myu" = rhat_myu,
              "rhat_sigma" = rhat_sigma,
              "rhat_lp" = rhat_lp,
              "traceplot" = traceplot))
}

#' perform_mcmc_sr_all_mk
#'
#' @description The function to perform MCMC for SR of all markers
perform_mcmc_sr_all_mk <- function(dt_sr_data, stutter_type, markers_auto, meth_model_sr, wlus_auto, wca_auto,
                                   model_sr_gen, model_sr_uni,
                                   warmup, sampling, chains, n_parallel_chains){
  if(stutter_type == "BS"){
    dt_st_type <- dt_sr_data[dt_sr_data$Bool_BS]
  }else if(stutter_type == "FS"){
    dt_st_type <- dt_sr_data[dt_sr_data$Bool_FS]
  }else if(stutter_type == "DS"){
    dt_st_type <- dt_sr_data[dt_sr_data$Bool_DS]
  }else if(stutter_type == "M2S"){
    dt_st_type <- dt_sr_data[dt_sr_data$Bool_M2S]
  }else{
    dt_st_type <- dt_sr_data[rep(FALSE, nrow(dt_sr_data))]
  }

  n_mk <- length(markers_auto)
  setindex(dt_st_type, Marker, Allele)
  median_myu_all_mk <- rep(NA, n_mk)
  median_slope_all_mk <- rep(NA, n_mk)
  median_intercept_all_mk <- rep(NA, n_mk)
  median_sigma_all_mk <- rep(NA, n_mk)
  median_lp_all_mk <- rep(NA, n_mk)
  rhat_myu_all_mk <- rep(NA, n_mk)
  rhat_slope_all_mk <- rep(NA, n_mk)
  rhat_intercept_all_mk <- rep(NA, n_mk)
  rhat_sigma_all_mk <- rep(NA, n_mk)
  rhat_lp_all_mk <- rep(NA, n_mk)
  trace_all_mk <- vector(mode = "list", length = n_mk)
  names(trace_all_mk) <- markers_auto
  dt_average_sr_all_mk <- vector(mode = "list", length = n_mk)
  names(dt_average_sr_all_mk) <- markers_auto

  index_multi_mk <- integer(0)
  y_sr_multi_mk <- numeric(0)
  tap_multi_mk <- numeric(0)
  tmp_alname_multi_mk <- vector(mode = "list", length = n_mk)
  tmp_alcor_multi_mk <- vector(mode = "list", length = n_mk)

  for(i in 1:n_mk){
    mk <- markers_auto[i]
    meth <- meth_model_sr[i]
    if(meth != "Not consider"){
      dt_st_mk <- dt_st_type[Marker == mk]
      dt_st_mk <- dt_st_mk[Allele != "Multiple"]
      dt_st_mk <- dt_st_mk[!dt_st_mk$Bool_SharedAllele]

      if(nrow(dt_st_mk) != 0){
        if(meth == "LUS"){
          dt_st_mk <- dt_st_mk[dt_st_mk$Bool_SeqData]
          st_alleles <- dt_st_mk$Allele

          alcor <- wlus_auto[[mk]]
          alname <- names(alcor)
          index_adopt <- which(is.element(st_alleles, alname))
          st_alleles <- st_alleles[index_adopt]

          index_alcor <- match(st_alleles, alname)
          x_al <- alcor[index_alcor]
          y_sr <- dt_st_mk$SR[index_adopt]
          tap <- dt_st_mk$TAP[index_adopt]
        }else if(meth == "Multi-seq"){
          dt_st_mk <- dt_st_mk[dt_st_mk$Bool_SeqData]
          st_alleles <- dt_st_mk$Allele

          alcor <- wca_auto[[mk]]
          alname <- names(alcor)
          index_adopt <- which(is.element(st_alleles, alname))
          st_alleles <- st_alleles[index_adopt]

          index_alcor <- match(st_alleles, alname)
          x_al <- alcor[index_alcor]
          y_sr <- dt_st_mk$SR[index_adopt]
          tap <- dt_st_mk$TAP[index_adopt]
        }else{
          st_alleles <- dt_st_mk$Allele
          x_al <- as.numeric(st_alleles)
          y_sr <- dt_st_mk$SR
          tap <- dt_st_mk$TAP

          tmp_alname <- c(names(wlus_auto[[mk]]), st_alleles)
          alcor <- sort(unique(as.numeric(tmp_alname)))
          alname <- as.character(alcor)
        }

        if(meth == "Uniform (Multiple loci)"){
          index_multi_mk <- c(index_multi_mk, i)
          y_sr_multi_mk <- c(y_sr_multi_mk, y_sr)
          tap_multi_mk <- c(tap_multi_mk, tap)
          tmp_alname_multi_mk[[i]] <- alname
          tmp_alcor_multi_mk[[i]] <- alcor
        }else if(meth == "Uniform"){
          data_list_sr <- list(N = length(x_al), SR = y_sr, Height = tap)
          tmp <- run_mcmc_sr_uni(model_sr_uni, data_list_sr, warmup, sampling, chains, n_parallel_chains)
          median_myu_all_mk[i] <- tmp$median_myu
          median_sigma_all_mk[i] <- tmp$median_sigma
          median_lp_all_mk[i] <- tmp$median_lp
          rhat_myu_all_mk[i] <- tmp$rhat_myu
          rhat_sigma_all_mk[i] <- tmp$rhat_sigma
          rhat_lp_all_mk[i] <- tmp$rhat_lp
          trace_all_mk[[i]] <- tmp$traceplot

          average_sr <- rep(median_myu_all_mk[i], length(alcor))
          average_sr[average_sr < 0] <- 0
          dt_average_sr_all_mk[[i]] <- data.table(Allele = alname,
                                                  AlCor = alcor,
                                                  AveSR = average_sr)
        }else{
          data_list_sr <- list(N = length(x_al), SR = y_sr, XVal = x_al, Height = tap)
          tmp <- run_mcmc_sr_gen(model_sr_gen, data_list_sr, warmup, sampling, chains, n_parallel_chains)
          median_slope_all_mk[i] <- tmp$median_slope
          median_intercept_all_mk[i] <- tmp$median_intercept
          median_sigma_all_mk[i] <- tmp$median_sigma
          median_lp_all_mk[i] <- tmp$median_lp
          rhat_slope_all_mk[i] <- tmp$rhat_slope
          rhat_intercept_all_mk[i] <- tmp$rhat_intercept
          rhat_sigma_all_mk[i] <- tmp$rhat_sigma
          rhat_lp_all_mk[i] <- tmp$rhat_lp
          trace_all_mk[[i]] <- tmp$traceplot

          average_sr <- calc_average_sr(alcor, median_slope_all_mk[i], median_intercept_all_mk[i])
          average_sr[average_sr < 0] <- 0
          dt_average_sr_all_mk[[i]] <- data.table(Allele = alname,
                                                  AlCor = alcor,
                                                  AveSR = average_sr)
        }
      }else{
        if(meth == "Uniform (Multiple loci)"){
          index_multi_mk <- c(index_multi_mk, i)
        }
      }
    }else{
      dt_average_sr_all_mk[[i]] <- data.table(Allele = character(0),
                                              AlCor = numeric(0),
                                              AveSR = numeric(0))
    }
  }

  # Uniform (Multiple loci)
  if(length(y_sr_multi_mk) != 0){
    data_list_sr <- list(N = length(y_sr_multi_mk), SR = y_sr_multi_mk, Height = tap_multi_mk)
    tmp <- run_mcmc_sr_uni(model_sr_uni, data_list_sr, warmup, sampling, chains, n_parallel_chains)
    median_myu_all_mk[index_multi_mk] <- tmp$median_myu
    median_sigma_all_mk[index_multi_mk] <- tmp$median_sigma
    median_lp_all_mk[index_multi_mk] <- tmp$median_lp
    rhat_myu_all_mk[index_multi_mk] <- tmp$rhat_myu
    rhat_sigma_all_mk[index_multi_mk] <- tmp$rhat_sigma
    rhat_lp_all_mk[index_multi_mk] <- tmp$rhat_lp

    for(i in 1:length(index_multi_mk)){
      index_mk <- index_multi_mk[i]

      trace_all_mk[[index_mk]] <- tmp$traceplot

      alname <- tmp_alname_multi_mk[[index_mk]]
      alcor <- tmp_alcor_multi_mk[[index_mk]]
      average_sr <- rep(median_myu_all_mk[index_mk], length(alcor))
      average_sr[average_sr < 0] <- 0
      dt_average_sr_all_mk[[index_mk]] <- data.table(Allele = alname,
                                                     AlCor = alcor,
                                                     AveSR = average_sr)
    }
  }

  dt_result_est_sr <- data.table(Marker = markers_auto,
                                 Model = meth_model_sr,
                                 Myu = median_myu_all_mk,
                                 Slope = median_slope_all_mk,
                                 Intercept = median_intercept_all_mk,
                                 Sigma = median_sigma_all_mk,
                                 Lp = median_lp_all_mk,
                                 Rhat_Myu = rhat_myu_all_mk,
                                 Rhat_Slope = rhat_slope_all_mk,
                                 Rhat_Intercept = rhat_intercept_all_mk,
                                 Rhat_Sigma = rhat_sigma_all_mk,
                                 Rhat_Lp = rhat_lp_all_mk)

  return(list("dt_result_est_sr" = dt_result_est_sr,
              "trace_all_mk" = trace_all_mk,
              "dt_average_sr_all_mk" = dt_average_sr_all_mk))
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
  meth_model_bsr <- dt_setting_est_sr[["BSR model"]]
  meth_model_fsr <- dt_setting_est_sr[["FSR model"]]
  meth_model_dsr <- dt_setting_est_sr[["DSR model"]]
  meth_model_m2sr <- dt_setting_est_sr[["M2SR model"]]

  # Make an object for ATs of all markers
  at_all_mk <- rep(NA, n_mk_all)
  names(at_all_mk) <- markers_all
  at_all_mk[match(markers_auto, markers_all)] <- at_auto
  at_all_mk[is.na(at_all_mk)] <- at_sex

  # Information on considering stutter ratios
  sr_consider_auto <- matrix(TRUE, n_mk_auto, 4)
  rownames(sr_consider_auto) <- markers_auto
  colnames(sr_consider_auto) <- c("BS", "FS", "DS", "M2S")
  sr_consider_auto[, 1] <- meth_model_bsr != "Not consider"
  sr_consider_auto[, 2] <- meth_model_fsr != "Not consider"
  sr_consider_auto[, 3] <- meth_model_dsr != "Not consider"
  sr_consider_auto[, 4] <- meth_model_m2sr != "Not consider"

  # Arrange experimental data
  tmp <- arrange_exp_data(dt_exp, dt_gt,
                          markers_auto, sr_consider_auto,
                          markers_all, at_all_mk, dye_all_mk,
                          pu_remove, pu_dist)
  peak_list <- tmp$peak_list
  size_list <- tmp$size_list
  height_list <- tmp$height_list
  gt_list <- tmp$gt_list

  # Prepare data related to the sequence data
  min_us_auto <- dt_setting_est_sr[["Minimum uninterrupted stretch (for Multi-seq)"]]
  tmp <- correct_al_rep(dt_seq, markers_auto, min_us_auto)
  wlus_auto <- tmp$wlus_auto
  wca_auto <- tmp$wca_auto

  # Create a data.table for the stutter ratios
  dt_sr_data <- create_dt_sr_data(peak_list, height_list, gt_list, sr_consider_auto, wlus_auto)

  # Perform MCMC for BSR
  tmp <- perform_mcmc_sr_all_mk(dt_sr_data, "BS", markers_auto, meth_model_bsr, wlus_auto, wca_auto,
                                model_sr_gen, model_sr_uni,
                                warmup, sampling, chains, n_parallel_chains)
  dt_result_est_bsr <- tmp$dt_result_est_sr
  list_trace_bsr <- tmp$trace_all_mk
  list_dt_average_bsr <- tmp$dt_average_sr_all_mk

  # Perform MCMC for FSR
  tmp <- perform_mcmc_sr_all_mk(dt_sr_data, "FS", markers_auto, meth_model_fsr, wlus_auto, wca_auto,
                                model_sr_gen, model_sr_uni,
                                warmup, sampling, chains, n_parallel_chains)
  dt_result_est_fsr <- tmp$dt_result_est_sr
  list_trace_fsr <- tmp$trace_all_mk
  list_dt_average_fsr <- tmp$dt_average_sr_all_mk

  # Perform MCMC for DSR
  tmp <- perform_mcmc_sr_all_mk(dt_sr_data, "DS", markers_auto, meth_model_dsr, wlus_auto, wca_auto,
                                model_sr_gen, model_sr_uni,
                                warmup, sampling, chains, n_parallel_chains)
  dt_result_est_dsr <- tmp$dt_result_est_sr
  list_trace_dsr <- tmp$trace_all_mk
  list_dt_average_dsr <- tmp$dt_average_sr_all_mk

  # Perform MCMC for M2SR
  tmp <- perform_mcmc_sr_all_mk(dt_sr_data, "M2S", markers_auto, meth_model_m2sr, wlus_auto, wca_auto,
                                model_sr_gen, model_sr_uni,
                                warmup, sampling, chains, n_parallel_chains)
  dt_result_est_m2sr <- tmp$dt_result_est_sr
  list_trace_m2sr <- tmp$trace_all_mk
  list_dt_average_m2sr <- tmp$dt_average_sr_all_mk

  return(list("dt_result_est_bsr" = dt_result_est_bsr,
              "list_trace_bsr" = list_trace_bsr,
              "list_dt_average_bsr" = list_dt_average_bsr,
              "dt_result_est_fsr" = dt_result_est_fsr,
              "list_trace_fsr" = list_trace_fsr,
              "list_dt_average_fsr" = list_dt_average_fsr,
              "dt_result_est_dsr" = dt_result_est_dsr,
              "list_trace_dsr" = list_trace_dsr,
              "list_dt_average_dsr" = list_dt_average_dsr,
              "dt_result_est_m2sr" = dt_result_est_m2sr,
              "list_trace_m2sr" = list_trace_m2sr,
              "list_dt_average_m2sr" = list_dt_average_m2sr))
}
