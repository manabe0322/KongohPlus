#' make_index_dedup2_mk
#'
#' @description The function to make 'index_dedup2_mk'
#' @param index_rph_gtc Information on the connection between rph_base and genotype combinations
#' @param index_dedup1_all_mk Index for deduplication (pattern 1)
make_index_dedup2_mk <- function(index_rph_gtc, index_dedup1_all_mk){
  paste_index_rph_gtc <- apply(index_rph_gtc, 2, paste, collapse = "-")
  paste_index_dedup1_all_mk <- apply(index_dedup1_all_mk, 2, paste, collapse = "-")
  return(match(paste_index_rph_gtc, paste_index_dedup1_all_mk))
}

#' perform_deconvo
#'
#' @description The function to perform mixture deconvolution
#' @param dt_csp The data.table of the crime stain profiles
#' @param dt_af The data.table of the allele frequencies
#' @param dt_kit The data.table of the information on the kit
#' @param dt_info_mk The data.table of the information on the ignored marker and the analytical thresholds for each marker
#' @param dt_par The data.table of the parameters
#' @param dt_sr The data.table of the stutter ratios
#' @param dt_al_prob The data.table of the setting of the allele probabilities
#' @param dt_gt_fltr The data.table of the genotype filters
#' @param dt_dropin The data.table of the setting of the drop-in
#' @param hnc The hypothesized number of contributors
#' @param warmup The number of warm-up iterations for each chain
#' @param sampling The number of sampling iterations for each chain
#' @param chains The number of chains
perform_deconvo <- function(dt_csp, dt_af, dt_kit, dt_info_mk, dt_par, dt_sr, dt_al_prob, dt_gt_fltr, dt_dropin,
                            hnc, warmup, sampling, chains){
  # Start computation
  start_time <- proc.time()

  # Number of replicates
  n_rep <- length(unique(dt_csp[, SampleName]))

  # Extract used markers and ATs
  markers_use <- dt_info_mk$Marker[!dt_info_mk$Ignore]
  ats_use <- dt_info_mk$AT[!dt_info_mk$Ignore]
  n_mk <- length(markers_use)

  if(n_mk != 0){
    dt_kit_auto <- extract_kit_auto(dt_kit)

    # Extract peak data from crime stain profiles
    tmp <- make_peak_data(dt_csp, dt_kit_auto, markers_use, ats_use)
    peaks_list <- tmp$peaks_list
    sizes_list <- tmp$sizes_list
    heights_list <- tmp$heights_list
    mean_size <- mean(dt_kit_auto$Size)

    # Extract parameters for AE
    ae_alpha_all_mk <- dt_par$Alpha[dt_par$Parameter == "AE"]
    ae_alpha_use <- ae_alpha_all_mk[is.element(dt_par$Marker[dt_par$Parameter == "AE"], markers_use)]

    # Extract methods to estimate allele probabilities
    al_prob_subcategory <- dt_al_prob[, Subcategory]
    al_prob_value <- dt_al_prob[, Value]
    method_al_prob <- al_prob_value[al_prob_subcategory == "Method"]
    maf <- as.numeric(al_prob_value[al_prob_subcategory == "Maf"])
    theta <- as.numeric(al_prob_value[al_prob_subcategory == "Theta"])

    # Set allele probabilities
    tmp <- set_al_prob(markers_use, dt_af, peaks_list, method_al_prob, maf)
    pop_al_list <- tmp$pop_al_list
    pop_prob_list <- tmp$pop_prob_list

    # Extract filters to remove unrealistic genotype combinations
    gt_fltr_subcategory <- dt_gt_fltr[, Subcategory]
    gt_fltr_value <- as.numeric(dt_gt_fltr[, Value])
    bsr_fltr <- gt_fltr_value[gt_fltr_subcategory == "BSR"]
    fsr_fltr <- gt_fltr_value[gt_fltr_subcategory == "FSR"]
    dsr_fltr <- gt_fltr_value[gt_fltr_subcategory == "DSR"]
    m2sr_fltr <- gt_fltr_value[gt_fltr_subcategory == "M2SR"]
    n_dropin_fltr <- gt_fltr_value[gt_fltr_subcategory == "N_Dropin"]
    h_dropin_fltr <- gt_fltr_value[gt_fltr_subcategory == "H_Dropin"]

    # Extract drop-in parameters
    dropin_subcategory <- dt_dropin[, Subcategory]
    dropin_value <- as.numeric(dt_dropin[, Value])
    logprob_di <- log(dropin_value[dropin_subcategory == "Probability"])
    shape_di <- dropin_value[dropin_subcategory == "Shape"]
    rate_di <- dropin_value[dropin_subcategory == "Rate"]
    mt_di <- dropin_value[dropin_subcategory == "MT"]

    # setindex of data.table
    setindex(dt_sr, Marker)

    # Define objects
    cand_gt_use <- vector(mode = "list", length = n_mk)
    index_gtc_adopt_cpp_use <- vector(mode = "list", length = n_mk)
    bool_gtc_di_mk_rep_use <- vector(mode = "list", length = n_mk)
    n_peak_use <- rep(NA, n_mk)
    n_gtc_use <- rep(NA, n_mk)

    for(i in 1:n_mk){
      # Extract peak data of one marker
      mk <- markers_use[i]
      peaks <- peaks_list[[mk]]
      sizes <- sizes_list[[mk]]
      heights_all_rep <- heights_list[[mk]]
      n_peak_use[i] <- n_peak <- length(peaks)

      # Extract parameters of one marker
      at <- ats_use[i]
      dt_sr_mk <- dt_sr[Marker == mk]
      sr_alleles_mk <- as.numeric(dt_sr_mk[, Allele])
      bsrs_mk <- dt_sr_mk[, Average_BSR]
      fsrs_mk <- dt_sr_mk[, Average_FSR]
      dsrs_mk <- dt_sr_mk[, Average_DSR]
      m2srs_mk <- dt_sr_mk[, Average_M2SR]
      index_sr <- match(peaks, sr_alleles_mk)
      bsrs <- fsrs <- dsrs <- m2srs <- rep(0, n_peak)
      for(j in 1:length(index_sr)){
        if(!is.na(index_sr[j])){
          bsrs[j] <- bsrs_mk[index_sr[j]]
          fsrs[j] <- fsrs_mk[index_sr[j]]
          dsrs[j] <- dsrs_mk[index_sr[j]]
          m2srs[j] <- m2srs_mk[index_sr[j]]
        }
      }

      # All candidate genotypes of one contributor
      cand_gt <- combinations(n = n_peak, r = 2, v = peaks, repeats.allowed = TRUE)
      cand_gt_use[[i]] <- cand_gt <- asplit(cand_gt, 1)

      # Indices of genotypes of each contributor in all candidate genotype combinations
      index_gtc <- permutations(n = length(cand_gt), r = hnc, repeats.allowed = TRUE)
      index_gtc_cpp <- index_gtc - 1 # C++ index
      index_gtc_cpp <- asplit(index_gtc_cpp, 1)

      # Indices of genotypes of each contributor in each adopted genotype combination
      tmp <- judge_gtc_mk(cand_gt, index_gtc_cpp, peaks, heights_all_rep,
                          bsrs, fsrs, dsrs, m2srs,
                          bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr, n_dropin_fltr, h_dropin_fltr)
      index_gtc_adopt_cpp_use[[i]] <- tmp[[1]]
      bool_gtc_di_mk_rep_use[[i]] <- tmp[[2]]

      n_gtc_use[i] <- n_gtc <- length(index_gtc_adopt_cpp_use[[i]])
    }

    if(all(n_gtc_use != 0)){
      ##############################
      # Objects for all replicates #
      ##############################

      max_tph <- 0
      rph_base_use <- vector(mode = "list", length = n_mk)
      n_rph_base_use <- rep(NA, n_mk)
      gtc_use <- vector(mode = "list", length = n_mk)
      names(gtc_use) <- markers_use
      gtc_logprob_use <- vector(mode = "list", length = n_mk)

      ###########################
      # Objects for replicate 1 #
      ###########################

      # Objects for general peaks
      n_general_use_1st <- rep(NA, n_mk)
      index_dedup1_general_use_1st <- NULL
      size_dedup1_general_use_1st <- numeric(0)
      oph_dedup1_general_use_1st <- numeric(0)
      id_mk_dedup1_general_use_1st <- integer(0)

      # Objects for drop-out peaks
      n_dropout_use_1st <- rep(NA, n_mk)
      index_dedup1_dropout_use_1st <- NULL
      size_dedup1_dropout_use_1st <- numeric(0)
      id_mk_dedup1_dropout_use_1st <- integer(0)

      # Objects for drop-in peaks
      n_dropin_use_1st <- rep(NA, n_mk)
      index_dedup1_dropin_use_1st <- NULL
      loglike_di_use_1st <- numeric(0)

      # Objects for ignored peaks
      n_ignore_use_1st <- rep(NA, n_mk)
      index_dedup1_ignore_use_1st <- NULL

      # Other objects
      index_rph_gtc_use_1st <- vector(mode = "list", length = n_mk)
      sumup_index_dedup1_1st <- 0

      ###########################
      # Objects for replicate 2 #
      ###########################

      # Objects for general peaks
      n_general_use_2nd <- rep(NA, n_mk)
      index_dedup1_general_use_2nd <- NULL
      size_dedup1_general_use_2nd <- numeric(0)
      oph_dedup1_general_use_2nd <- numeric(0)
      id_mk_dedup1_general_use_2nd <- integer(0)

      # Objects for drop-out peaks
      n_dropout_use_2nd <- rep(NA, n_mk)
      index_dedup1_dropout_use_2nd <- NULL
      size_dedup1_dropout_use_2nd <- numeric(0)
      id_mk_dedup1_dropout_use_2nd <- integer(0)

      # Objects for drop-in peaks
      n_dropin_use_2nd <- rep(NA, n_mk)
      index_dedup1_dropin_use_2nd <- NULL
      loglike_di_use_2nd <- numeric(0)

      # Objects for ignored peaks
      n_ignore_use_2nd <- rep(NA, n_mk)
      index_dedup1_ignore_use_2nd <- NULL

      # Other objects
      index_rph_gtc_use_2nd <- vector(mode = "list", length = n_mk)
      sumup_index_dedup1_2nd <- 0

      for(i in 1:n_mk){
        # Extract peak data of one marker
        mk <- markers_use[i]
        peaks <- peaks_list[[mk]]
        sizes <- sizes_list[[mk]]
        heights_all_rep <- heights_list[[mk]]
        n_peak <- n_peak_use[i]

        # Extract parameters of one marker
        at <- ats_use[i]
        dt_sr_mk <- dt_sr[Marker == mk]
        sr_alleles_mk <- as.numeric(dt_sr_mk[, Allele])
        bsrs_mk <- dt_sr_mk[, Average_BSR]
        fsrs_mk <- dt_sr_mk[, Average_FSR]
        dsrs_mk <- dt_sr_mk[, Average_DSR]
        m2srs_mk <- dt_sr_mk[, Average_M2SR]
        index_sr <- match(peaks, sr_alleles_mk)
        bsrs <- fsrs <- dsrs <- m2srs <- rep(0, n_peak)
        for(j in 1:length(index_sr)){
          if(!is.na(index_sr[j])){
            bsrs[j] <- bsrs_mk[index_sr[j]]
            fsrs[j] <- fsrs_mk[index_sr[j]]
            dsrs[j] <- dsrs_mk[index_sr[j]]
            m2srs[j] <- m2srs_mk[index_sr[j]]
          }
        }

        # The base of the relative peak heights
        rph_base_use[[i]] <- rph_base_mk <- make_rph_base(peaks, bsrs, fsrs, dsrs, m2srs)
        n_rph_base_use[i] <- length(rph_base_mk)

        # Extract genotype combinations
        n_gtc <- n_gtc_use[i]
        cand_gt <- cand_gt_use[[i]]
        index_gtc_adopt_cpp <- index_gtc_adopt_cpp_use[[i]]
        gtc_use[[i]] <- gtc_mk <- extract_gtc_adopt(cand_gt, index_gtc_adopt_cpp)

        # Extract allele probabilities
        pop_al <- pop_al_list[[mk]]
        pop_prob <- pop_prob_list[[mk]]

        # Calculate genotype probabilities
        n_obs <- make_n_obs(pop_al, prof_ref = numeric(0))
        gtc_prob_mk <- calc_gtc_prob_mk(gtc_mk, pop_al, pop_prob, theta, prof_ref = numeric(0), hyp_id = rep(0, hnc), n_obs)
        gtc_logprob_use[[i]] <- log(gtc_prob_mk)

        ###############
        # Replicate 1 #
        ###############

        id_rep <- 1

        # Connect rph_base and genotype combinations
        index_rph_gtc_cpp <- connect_rph_gtc_mk(cand_gt, index_gtc_adopt_cpp, peaks, bsrs, fsrs, dsrs, m2srs, bool_gtc_di_mk_rep_use[[i]][[id_rep]])
        index_rph_gtc_use_1st[[i]] <- index_rph_gtc <- matrix(unlist(index_rph_gtc_cpp) + 1, nrow = 2 * hnc) # R and stan index

        # Delete duplicated indices -> index for deduplication (pattern 1: Before summation of the rph_base)
        bool_duplicated_adopt <- !duplicated(t(index_rph_gtc))
        index_dedup1_tmp <- index_rph_gtc[, bool_duplicated_adopt, drop = FALSE]

        # Sizes for de-duplicated objects (pattern 1)
        size_dedup1_tmp <- rep(sizes, n_gtc)
        size_dedup1_tmp <- size_dedup1_tmp[bool_duplicated_adopt]

        # Update the maximum of the total peak height per marker
        heights <- heights_all_rep[[id_rep]]
        oph_sum <- sum(heights)
        if(oph_sum > max_tph){
          max_tph <- oph_sum
        }

        # Observe peak heights for de-duplicated objects (pattern 1)
        oph_dedup1_tmp <- rep(heights, n_gtc)
        oph_dedup1_tmp <- oph_dedup1_tmp[bool_duplicated_adopt]

        # Indices of each peak type
        rph_dedup1_tmp <- matrix(rph_base_mk[index_dedup1_tmp], nrow = 2 * hnc)
        rph_dedup1_tmp <- apply(rph_dedup1_tmp, 2, sum)
        index_oph0 <- which(oph_dedup1_tmp == 0)
        index_ignore <- c(which(rph_dedup1_tmp == 0), intersect(which(rph_dedup1_tmp < 0), index_oph0)) # intersect(which(rph_dedup1_tmp < 0), index_oph0): for duplicate
        index_dropin <- setdiff(which(rph_dedup1_tmp < 0), index_ignore)
        index_dropout <- setdiff(index_oph0, c(index_ignore, index_dropin))
        index_general <- setdiff(1:length(rph_dedup1_tmp), c(index_ignore, index_dropin, index_dropout))

        # Assign values for general peaks
        n_general_use_1st[i] <- length(index_general)
        if(n_general_use_1st[i] != 0){
          index_dedup1_general_use_1st <- cbind(index_dedup1_general_use_1st, index_dedup1_tmp[, index_general, drop = FALSE] + sumup_index_dedup1_1st)
          size_dedup1_general_use_1st <- c(size_dedup1_general_use_1st, size_dedup1_tmp[index_general])
          oph_dedup1_general_use_1st <- c(oph_dedup1_general_use_1st, oph_dedup1_tmp[index_general])
          id_mk_dedup1_general_use_1st <- c(id_mk_dedup1_general_use_1st, rep(i, n_general_use_1st[i]))
        }

        # Assign values for drop-out peaks
        n_dropout_use_1st[i] <- length(index_dropout)
        if(n_dropout_use_1st[i] != 0){
          index_dedup1_dropout_use_1st <- cbind(index_dedup1_dropout_use_1st, index_dedup1_tmp[, index_dropout, drop = FALSE] + sumup_index_dedup1_1st)
          size_dedup1_dropout_use_1st <- c(size_dedup1_dropout_use_1st, size_dedup1_tmp[index_dropout])
          id_mk_dedup1_dropout_use_1st <- c(id_mk_dedup1_dropout_use_1st, rep(i, n_dropout_use_1st[i]))
        }

        # Assign values for drop-in peaks
        n_dropin_use_1st[i] <- length(index_dropin)
        if(n_dropin_use_1st[i] != 0){
          index_dedup1_dropin_use_1st <- cbind(index_dedup1_dropin_use_1st, index_dedup1_tmp[, index_dropin, drop = FALSE] + sumup_index_dedup1_1st)
          loglike_di_use_1st <- c(loglike_di_use_1st, log(dgamma(oph_dedup1_tmp[index_dropin] - mt_di + 1, shape_di, rate_di)) + logprob_di)
        }

        # Assign values for ignored peaks
        n_ignore_use_1st[i] <- length(index_ignore)
        if(n_ignore_use_1st[i] != 0){
          index_dedup1_ignore_use_1st <- cbind(index_dedup1_ignore_use_1st, index_dedup1_tmp[, index_ignore, drop = FALSE] + sumup_index_dedup1_1st)
        }

        sumup_index_dedup1_1st <- sumup_index_dedup1_1st + n_rph_base_use[i]

        ###############
        # Replicate 2 #
        ###############

        if(n_rep != 1){
          id_rep <- 2

          # Connect rph_base and genotype combinations
          index_rph_gtc_cpp <- connect_rph_gtc_mk(cand_gt, index_gtc_adopt_cpp, peaks, bsrs, fsrs, dsrs, m2srs, bool_gtc_di_mk_rep_use[[i]][[id_rep]])
          index_rph_gtc_use_2nd[[i]] <- index_rph_gtc <- matrix(unlist(index_rph_gtc_cpp) + 1, nrow = 2 * hnc) # R and stan index

          # Delete duplicated indices -> index for deduplication (pattern 1: Before summation of the rph_base)
          bool_duplicated_adopt <- !duplicated(t(index_rph_gtc))
          index_dedup1_tmp <- index_rph_gtc[, bool_duplicated_adopt, drop = FALSE]

          # Sizes for de-duplicated objects (pattern 1)
          size_dedup1_tmp <- rep(sizes, n_gtc)
          size_dedup1_tmp <- size_dedup1_tmp[bool_duplicated_adopt]

          # Update the maximum of the total peak height per marker
          heights <- heights_all_rep[[id_rep]]
          oph_sum <- sum(heights)
          if(oph_sum > max_tph){
            max_tph <- oph_sum
          }

          # Observe peak heights for de-duplicated objects (pattern 1)
          oph_dedup1_tmp <- rep(heights, n_gtc)
          oph_dedup1_tmp <- oph_dedup1_tmp[bool_duplicated_adopt]

          # Indices of each peak type
          rph_dedup1_tmp <- matrix(rph_base_mk[index_dedup1_tmp], nrow = 2 * hnc)
          rph_dedup1_tmp <- apply(rph_dedup1_tmp, 2, sum)
          index_oph0 <- which(oph_dedup1_tmp == 0)
          index_ignore <- c(which(rph_dedup1_tmp == 0), intersect(which(rph_dedup1_tmp < 0), index_oph0)) # intersect(which(rph_dedup1_tmp < 0), index_oph0): for duplicate
          index_dropin <- setdiff(which(rph_dedup1_tmp < 0), index_ignore)
          index_dropout <- setdiff(index_oph0, c(index_ignore, index_dropin))
          index_general <- setdiff(1:length(rph_dedup1_tmp), c(index_ignore, index_dropin, index_dropout))

          # Assign values for general peaks
          n_general_use_2nd[i] <- length(index_general)
          if(n_general_use_2nd[i] != 0){
            index_dedup1_general_use_2nd <- cbind(index_dedup1_general_use_2nd, index_dedup1_tmp[, index_general, drop = FALSE] + sumup_index_dedup1_2nd)
            size_dedup1_general_use_2nd <- c(size_dedup1_general_use_2nd, size_dedup1_tmp[index_general])
            oph_dedup1_general_use_2nd <- c(oph_dedup1_general_use_2nd, oph_dedup1_tmp[index_general])
            id_mk_dedup1_general_use_2nd <- c(id_mk_dedup1_general_use_2nd, rep(i, n_general_use_2nd[i]))
          }

          # Assign values for drop-out peaks
          n_dropout_use_2nd[i] <- length(index_dropout)
          if(n_dropout_use_2nd[i] != 0){
            index_dedup1_dropout_use_2nd <- cbind(index_dedup1_dropout_use_2nd, index_dedup1_tmp[, index_dropout, drop = FALSE] + sumup_index_dedup1_2nd)
            size_dedup1_dropout_use_2nd <- c(size_dedup1_dropout_use_2nd, size_dedup1_tmp[index_dropout])
            id_mk_dedup1_dropout_use_2nd <- c(id_mk_dedup1_dropout_use_2nd, rep(i, n_dropout_use_2nd[i]))
          }

          # Assign values for drop-in peaks
          n_dropin_use_2nd[i] <- length(index_dropin)
          if(n_dropin_use_2nd[i] != 0){
            index_dedup1_dropin_use_2nd <- cbind(index_dedup1_dropin_use_2nd, index_dedup1_tmp[, index_dropin, drop = FALSE] + sumup_index_dedup1_2nd)
            loglike_di_use_2nd <- c(loglike_di_use_2nd, log(dgamma(oph_dedup1_tmp[index_dropin] - mt_di + 1, shape_di, rate_di)) + logprob_di)
          }

          # Assign values for ignored peaks
          n_ignore_use_2nd[i] <- length(index_ignore)
          if(n_ignore_use_2nd[i] != 0){
            index_dedup1_ignore_use_2nd <- cbind(index_dedup1_ignore_use_2nd, index_dedup1_tmp[, index_ignore, drop = FALSE] + sumup_index_dedup1_2nd)
          }

          sumup_index_dedup1_2nd <- sumup_index_dedup1_2nd + n_rph_base_use[i]

        }
      }

      #####################################
      # De-duplicated objects (pattern 2) #
      #####################################

      # Objects for all replicates
      n_dedup2_use <- n_peak_use * n_gtc_use
      n_dedup2_sum <- sum(n_dedup2_use)

      # Objects for replicate 1
      index_dedup1_use_1st <- cbind(index_dedup1_general_use_1st, index_dedup1_dropout_use_1st, index_dedup1_dropin_use_1st, index_dedup1_ignore_use_1st)
      index_dedup2_use_1st <- rep(NA, n_dedup2_sum)

      # Objects for replicate 2
      index_dedup1_use_2nd <- cbind(index_dedup1_general_use_2nd, index_dedup1_dropout_use_2nd, index_dedup1_dropin_use_2nd, index_dedup1_ignore_use_2nd)
      index_dedup2_use_2nd <- rep(NA, n_dedup2_sum)

      ###############
      # Replicate 1 #
      ###############

      index_dedup2_start <- 1
      sumup_index_dedup2 <- 0
      for(i in 1:n_mk){
        index_rph_gtc <- index_rph_gtc_use_1st[[i]] + sumup_index_dedup2
        n_dedup2 <- n_dedup2_use[i]
        index_dedup2_use_1st[index_dedup2_start:(index_dedup2_start + n_dedup2 - 1)] <- make_index_dedup2_mk(index_rph_gtc, index_dedup1_use_1st)

        index_dedup2_start <- index_dedup2_start + n_dedup2
        sumup_index_dedup2 <- sumup_index_dedup2 + n_rph_base_use[i]
      }

      ###############
      # Replicate 2 #
      ###############

      if(n_rep != 1){
        index_dedup2_start <- 1
        sumup_index_dedup2 <- 0
        for(i in 1:n_mk){
          index_rph_gtc <- index_rph_gtc_use_2nd[[i]] + sumup_index_dedup2
          n_dedup2 <- n_dedup2_use[i]
          index_dedup2_use_2nd[index_dedup2_start:(index_dedup2_start + n_dedup2 - 1)] <- make_index_dedup2_mk(index_rph_gtc, index_dedup1_use_2nd)

          index_dedup2_start <- index_dedup2_start + n_dedup2
          sumup_index_dedup2 <- sumup_index_dedup2 + n_rph_base_use[i]
        }
      }

      ########
      # MCMC #
      ########

      var_allele <- 1141481300 #616764
      var_stutter <- 1731569 #5258
      beta_allele <- 0.0003141005 #0.0256
      beta_stutter <- 0.004219315 #0.0763
      cv_allele <- 10.9
      cv_stutter <- 4.16

      if(n_rep == 1){
        data_list_main <- list(hnc = hnc,
                               n_mk = n_mk,
                               ltph = max_tph / 4,
                               htph = max_tph / 2,
                               n_rph_base_sum = sum(n_rph_base_use),
                               n_rph_base_use = n_rph_base_use,
                               rph_base_use = unlist(rph_base_use),
                               ae_alpha_use = ae_alpha_use,
                               ats_use = ats_use,
                               mean_size = mean_size,
                               cv_allele = cv_allele, #var_allele = var_allele,
                               cv_stutter = cv_stutter, #var_stutter = var_stutter,
                               n_dedup1_sum = ncol(index_dedup1_use_1st),
                               n_general_sum = sum(n_general_use_1st),
                               n_general_use = n_general_use_1st,
                               index_dedup1_general = index_dedup1_general_use_1st,
                               size_dedup1_general = size_dedup1_general_use_1st,
                               oph_dedup1_general = oph_dedup1_general_use_1st,
                               id_mk_dedup1_general = id_mk_dedup1_general_use_1st,
                               n_dropout_sum = sum(n_dropout_use_1st),
                               n_dropout_use = n_dropout_use_1st,
                               index_dedup1_dropout = index_dedup1_dropout_use_1st,
                               size_dedup1_dropout = size_dedup1_dropout_use_1st,
                               id_mk_dedup1_dropout = id_mk_dedup1_dropout_use_1st,
                               n_dropin_sum = sum(n_dropin_use_1st),
                               loglike_di = loglike_di_use_1st,
                               n_dedup2_sum = n_dedup2_sum,
                               n_dedup2_use = n_dedup2_use,
                               index_dedup2_use = index_dedup2_use_1st,
                               n_peak_use = n_peak_use,
                               n_gtc_use = n_gtc_use,
                               gtc_logprob_use = unlist(gtc_logprob_use))

        # Perform MCMC
        model_main <- stan_package_model(name = "model_main", package = "KongohPlus")
        fit_main <- model_main$sample(
          data_list_main,
          iter_warmup = warmup,
          iter_sampling = sampling,
          chains = chains,
          parallel_chains = detectCores(),
          adapt_delta = 0.99,
          max_treedepth = 15,
          step_size = 1
        )

        # Trace plot
        draws <- fit_main$draws()
        draws %>% mcmc_trace(pars = c("lp__")) -> traceplot

        # Extract parameters
        fit_main_summary <- fit_main$summary()
        variables <- fit_main_summary$variable
        medians <- fit_main_summary$median
        rhats <- fit_main_summary$rhat
        lp <- medians[which(variables == "lp__")]
        tph <- medians[which(variables == "tph")]
        lambda <- medians[grep("lambda", variables)]
        mp <- medians[grep("mp", variables)]
        deg <- medians[grep("deg", variables)]
        pae <- medians[grep("pae", variables)]
        ae <- n_mk * pae
        rhat_lp <- rhats[which(variables == "lp__")]
        rhat_tph <- rhats[which(variables == "tph")]
        rhat_lambda <- rhats[grep("lambda", variables)]
        rhat_mp <- rhats[grep("mp", variables)]
        rhat_deg <- rhats[grep("deg", variables)]
        rhat_pae <- rhats[grep("pae", variables)]

        # Calculate likelihood under the estimated parameter
        index_dedup1_use_cpp <- asplit(index_dedup1_use_1st - 1, 1) # C++ index
        index_dedup2_use_cpp <- index_dedup2_use_1st - 1 # C++ index
        ae_dedup1_use <- ae[id_mk_dedup1_1st]
        lpd_gtc_use_1st <- calc_lpd_gtc_all_mk_cpp(n_mk, mp, deg, hnc, n_dedup1_sum_1st, n_general_sum_1st, n_dropout_sum_1st, n_dropin_sum_1st,
                                                   unlist(rph_base_use), index_dedup1_use_cpp, unlist(size_dedup1_use_1st), mean_size, ae_dedup1_use, oph_dedup1_use_1st,
                                                   at_dropout_1st, loglike_di_1st, cv, tph,
                                                   index_dedup2_use_cpp, n_dedup2_use, n_peak_use, n_gtc_use)
        lpd_gtc_use_2nd <- NULL
        log10like_use <- calc_log10like_use(gtc_logprob_use, lpd_gtc_use_1st)

        dt_pg_list <- make_dt_pg_list(markers_use, paste0("U", 1:hnc), gtc_use, gtc_logprob_use, lpd_gtc_use_1st)

        dt_mcmc <- data.table(Parameter = c("Log probability in stan", "Total peak height", rep("lambda", hnc), rep("Mixture proportion", hnc), rep("Degradation", hnc), rep("Amplification efficiency", n_mk)),
                              Contributor = c("-", "-", paste0("U", 1:hnc), paste0("U", 1:hnc), paste0("U", 1:hnc), rep("-", n_mk)),
                              Marker = c("-", "-", rep("-", 3 * hnc), markers_use),
                              Estimated_value = c(lp, tph, lambda[hnc:1], mp[hnc:1], deg[hnc:1], ae),
                              Rhat = c(rhat_lp, rhat_tph, rhat_lambda[hnc:1], rhat_mp[hnc:1], rhat_deg[hnc:1], rhat_pae))
      }else{
        data_list_main <- list(hnc = hnc,
                               n_mk = n_mk,
                               ltph = max_tph / 4,
                               htph = max_tph / 2,
                               n_rph_base_sum = sum(n_rph_base_use),
                               rph_base_use = unlist(rph_base_use),
                               ae_alpha_use = ae_alpha_use,
                               cv = cv,
                               mean_size = mean_size,
                               n_dedup1_sum_1st = n_dedup1_sum_1st,
                               index_dedup1_use_1st = index_dedup1_use_1st,
                               size_dedup1_use_1st = size_dedup1_use_1st,
                               oph_dedup1_use_1st = oph_dedup1_use_1st,
                               id_mk_dedup1_1st = id_mk_dedup1_1st,
                               n_general_sum_1st = n_general_sum_1st,
                               n_dropout_sum_1st = n_dropout_sum_1st,
                               at_dropout_1st = at_dropout_1st,
                               n_dropin_sum_1st = n_dropin_sum_1st,
                               loglike_di_1st = loglike_di_1st,
                               n_dedup1_sum_2nd = n_dedup1_sum_2nd,
                               index_dedup1_use_2nd = index_dedup1_use_2nd,
                               size_dedup1_use_2nd = size_dedup1_use_2nd,
                               oph_dedup1_use_2nd = oph_dedup1_use_2nd,
                               id_mk_dedup1_2nd = id_mk_dedup1_2nd,
                               n_general_sum_2nd = n_general_sum_2nd,
                               n_dropout_sum_2nd = n_dropout_sum_2nd,
                               at_dropout_2nd = at_dropout_2nd,
                               n_dropin_sum_2nd = n_dropin_sum_2nd,
                               loglike_di_2nd = loglike_di_2nd,
                               n_dedup2_sum = n_dedup2_sum,
                               n_dedup2_use = n_dedup2_use,
                               index_dedup2_use_1st = index_dedup2_use_1st,
                               index_dedup2_use_2nd = index_dedup2_use_2nd,
                               n_peak_use = n_peak_use,
                               n_gtc_use = n_gtc_use,
                               gtc_logprob_use = unlist(gtc_logprob_use))

        # Perform MCMC
        model_main_duplicate <- stan_package_model(name = "model_main_duplicate", package = "KongohPlus")
        fit_main <- model_main_duplicate$sample(
          data_list_main,
          iter_warmup = warmup,
          iter_sampling = sampling,
          chains = chains,
          parallel_chains = detectCores()
        )

        # Trace plot
        draws <- fit_main$draws()
        draws %>% mcmc_trace(pars = c("lp__")) -> traceplot

        # Extract parameters
        fit_main_summary <- fit_main$summary()
        variables <- fit_main_summary$variable
        medians <- fit_main_summary$median
        rhats <- fit_main_summary$rhat
        lp <- medians[which(variables == "lp__")]
        tph_1st <- medians[which(variables == "tph_1st")]
        tph_2nd <- medians[which(variables == "tph_2nd")]
        lambda <- medians[grep("lambda", variables)]
        mp <- medians[grep("mp", variables)]
        deg <- medians[grep("deg", variables)]
        pae <- medians[grep("pae", variables)]
        ae <- n_mk * pae
        rhat_lp <- rhats[which(variables == "lp__")]
        rhat_tph_1st <- rhats[which(variables == "tph_1st")]
        rhat_tph_2nd <- rhats[which(variables == "tph_2nd")]
        rhat_lambda <- rhats[grep("lambda", variables)]
        rhat_mp <- rhats[grep("mp", variables)]
        rhat_deg <- rhats[grep("deg", variables)]
        rhat_pae <- rhats[grep("pae", variables)]

        # Calculate likelihood under the estimated parameter
        index_dedup1_use_cpp_1st <- asplit(index_dedup1_use_1st - 1, 1) # C++ index
        index_dedup2_use_cpp_1st <- index_dedup2_use_1st - 1 # C++ index
        ae_dedup1_use_1st <- ae[id_mk_dedup1_1st]
        lpd_gtc_use_1st <- calc_lpd_gtc_all_mk_cpp(n_mk, mp, deg, hnc, n_dedup1_sum_1st, n_general_sum_1st, n_dropout_sum_1st, n_dropin_sum_1st,
                                                   unlist(rph_base_use), index_dedup1_use_cpp_1st, unlist(size_dedup1_use_1st), mean_size, ae_dedup1_use_1st, oph_dedup1_use_1st,
                                                   at_dropout_1st, loglike_di_1st, cv, tph_1st,
                                                   index_dedup2_use_cpp_1st, n_dedup2_use, n_peak_use, n_gtc_use)
        index_dedup1_use_cpp_2nd <- asplit(index_dedup1_use_2nd - 1, 1) # C++ index
        index_dedup2_use_cpp_2nd <- index_dedup2_use_2nd - 1 # C++ index
        ae_dedup1_use_2nd <- ae[id_mk_dedup1_2nd]
        lpd_gtc_use_2nd <- calc_lpd_gtc_all_mk_cpp(n_mk, mp, deg, hnc, n_dedup1_sum_2nd, n_general_sum_2nd, n_dropout_sum_2nd, n_dropin_sum_2nd,
                                                   unlist(rph_base_use), index_dedup1_use_cpp_2nd, unlist(size_dedup1_use_2nd), mean_size, ae_dedup1_use_2nd, oph_dedup1_use_2nd,
                                                   at_dropout_2nd, loglike_di_2nd, cv, tph_2nd,
                                                   index_dedup2_use_cpp_2nd, n_dedup2_use, n_peak_use, n_gtc_use)
        log10like_use <- calc_log10like_use(gtc_logprob_use, lpd_gtc_use_1st, lpd_gtc_use_2nd)

        dt_pg_list <- make_dt_pg_list(markers_use, paste0("U", 1:hnc), gtc_use, gtc_logprob_use, lpd_gtc_use_1st, lpd_gtc_use_2nd)

        dt_mcmc <- data.table(Parameter = c("Log probability in stan", "Total peak height (1st)", "Total peak height (2nd)", rep("lambda", hnc), rep("Mixture proportion", hnc), rep("Degradation", hnc), rep("Amplification efficiency", n_mk)),
                              Contributor = c("-", "-", "-", paste0("U", 1:hnc), paste0("U", 1:hnc), paste0("U", 1:hnc), rep("-", n_mk)),
                              Marker = c("-", "-", "-", rep("-", 3 * hnc), markers_use),
                              Estimated_value = c(lp, tph_1st, tph_2nd, lambda[hnc:1], mp[hnc:1], deg[hnc:1], ae),
                              Rhat = c(rhat_lp, rhat_tph_1st, rhat_tph_2nd, rhat_lambda[hnc:1], rhat_mp[hnc:1], rhat_deg[hnc:1], rhat_pae))
      }

      dt_log10like <- data.table(Marker = c(markers_use, "Total"),
                                 Log10_Likelihood = c(log10like_use, sum(log10like_use)))

      dt_mpdeg <- data.table(Contributor = paste0("U", 1:hnc),
                             MP = mp[hnc:1],
                             Deg = deg[hnc:1])

      # End computation
      end_time <- proc.time()
      comp_time <- end_time[3] - start_time[3]

      return(list("dt_csp" = dt_csp,
                  "dt_af" = dt_af,
                  "dt_kit" = dt_kit,
                  "dt_info_mk" = dt_info_mk,
                  "dt_par" = dt_par,
                  "dt_sr" = dt_sr,
                  "dt_al_prob" = dt_al_prob,
                  "dt_gt_fltr" = dt_gt_fltr,
                  "dt_dropin" = dt_dropin,
                  "warmup" = warmup,
                  "sampling" = sampling,
                  "chains" = chains,
                  "hnc" = hnc,
                  "markers_use" = markers_use,
                  "peaks_list" = peaks_list,
                  "pop_al_list" = pop_al_list,
                  "pop_prob_list" = pop_prob_list,
                  "theta" = theta,
                  "gtc_use" = gtc_use,
                  "lpd_gtc_use_1st" = lpd_gtc_use_1st,
                  "lpd_gtc_use_2nd" = lpd_gtc_use_2nd,
                  "dt_pg_list" = dt_pg_list,
                  "dt_log10like" = dt_log10like,
                  "dt_mpdeg" = dt_mpdeg,
                  "dt_mcmc" = dt_mcmc,
                  "traceplot" = traceplot,
                  "comp_time" = comp_time,
                  "error_message" = NULL))
    }else{
      markers_unexplainable <- markers_use[n_gtc_use == 0]
      if(hnc == 1){
        error_message <- paste0("The following marker(s) cannot be explained under ", hnc, " contributor.<br>")
      }else{
        error_message <- paste0("The following marker(s) cannot be explained under ", hnc, " contributors.<br>")
      }
      for(i in 1:length(markers_unexplainable)){
        error_message <- paste0(error_message, "  - ", markers_unexplainable[i], "<br>")
      }

      return(list("dt_csp" = NULL,
                  "dt_af" = NULL,
                  "dt_kit" = NULL,
                  "dt_info_mk" = NULL,
                  "dt_par" = NULL,
                  "dt_sr" = NULL,
                  "dt_al_prob" = NULL,
                  "dt_gt_fltr" = NULL,
                  "dt_dropin" = NULL,
                  "warmup" = NULL,
                  "sampling" = NULL,
                  "chains" = NULL,
                  "hnc" = NULL,
                  "markers_use" = NULL,
                  "peaks_list" = NULL,
                  "pop_al_list" = NULL,
                  "pop_prob_list" = NULL,
                  "theta" = NULL,
                  "gtc_use" = NULL,
                  "lpd_gtc_use_1st" = NULL,
                  "lpd_gtc_use_2nd" = NULL,
                  "dt_pg_list" = NULL,
                  "dt_log10like" = NULL,
                  "dt_mpdeg" = NULL,
                  "dt_mcmc" = NULL,
                  "traceplot" = NULL,
                  "comp_time" = NULL,
                  "error_message" = error_message))
    }
  }else{
    return(list("dt_csp" = NULL,
                "dt_af" = NULL,
                "dt_kit" = NULL,
                "dt_info_mk" = NULL,
                "dt_par" = NULL,
                "dt_sr" = NULL,
                "dt_al_prob" = NULL,
                "dt_gt_fltr" = NULL,
                "dt_dropin" = NULL,
                "warmup" = NULL,
                "sampling" = NULL,
                "chains" = NULL,
                "hnc" = NULL,
                "markers_use" = NULL,
                "peaks_list" = NULL,
                "pop_al_list" = NULL,
                "pop_prob_list" = NULL,
                "theta" = NULL,
                "gtc_use" = NULL,
                "lpd_gtc_use_1st" = NULL,
                "lpd_gtc_use_2nd" = NULL,
                "dt_pg_list" = NULL,
                "dt_log10like" = NULL,
                "dt_mpdeg" = NULL,
                "dt_mcmc" = NULL,
                "traceplot" = NULL,
                "comp_time" = NULL,
                "error_message" = "All markers are ignored!"))
  }
}
