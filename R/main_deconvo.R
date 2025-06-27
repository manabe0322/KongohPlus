#' create_gt_perm_obj
#'
#' @description The function to create objects for genotype permutations
#' @param cand_gt A list of all possible candidate genotypes
#' @param index_gt_comb_adopt_cpp A list of the indices of the adopted genotype combinations (C++ index)
#' @param bool_gt_comb_di_mk_rep A list of indicators for drop-in peaks in each genotype combination
create_gt_perm_obj <- function(cand_gt, index_gt_comb_adopt_cpp, log_gt_comb_prob_mk, bool_gt_comb_di_mk_rep){
  list_gt_perm <- list()
  list_log_gt_perm_prob_mk <- list()

  n_rep <- length(bool_gt_comb_di_mk_rep)
  bool_gt_perm_di_mk_rep <- vector(mode = "list", length = n_rep)
  for(i in 1:n_rep){
    bool_gt_perm_di_mk_rep[[i]] <- list()
  }

  for(i in 1:length(index_gt_comb_adopt_cpp)){
    index_gt_perm_cpp <- permn(index_gt_comb_adopt_cpp[[i]]) # C++ index
    tmp <- extract_gts_adopt(cand_gt, index_gt_perm_cpp)
    n_perm <- length(tmp)
    list_gt_perm[[i]] <- tmp
    list_log_gt_perm_prob_mk[[i]] <- rep(log_gt_comb_prob_mk[i], n_perm)

    for(j in 1:n_rep){
      bool_gt_comb_di <- bool_gt_comb_di_mk_rep[[j]][[i]]
      count <- length(bool_gt_perm_di_mk_rep[[j]])
      for(k in 1:n_perm){
        bool_gt_perm_di_mk_rep[[j]][[count + k]] <- bool_gt_comb_di
      }
    }
  }

  gt_perm <- unlist(list_gt_perm, recursive = FALSE, use.names = FALSE)
  log_gt_perm_prob_mk <- unlist(list_log_gt_perm_prob_mk)

  return(list(gt_perm = gt_perm, log_gt_perm_prob_mk = log_gt_perm_prob_mk, bool_gt_perm_di_mk_rep = bool_gt_perm_di_mk_rep))
}

#' create_index_dedup2_mk
#'
#' @description The function to create 'index_dedup2_mk'
#' @param index_rph_gtc_cpp C++ indices of the relative peak heights for each genotype combinations
#' @param index_dedup1_rep1_mk1 C++ indices of the relative peak heights for de-duplicated objects (pattern 1)
create_index_dedup2_mk <- function(index_rph_gtc_cpp, index_dedup1_rep1_mk1){
  paste_index_rph_gtc_cpp <- sapply(index_rph_gtc_cpp, paste, collapse = "-")
  paste_index_dedup1_rep1_mk1 <- sapply(index_dedup1_rep1_mk1, paste, collapse = "-")
  return(match(paste_index_rph_gtc_cpp, paste_index_dedup1_rep1_mk1) - 1)
}

#' generate_init_deconvo_single
#'
#' @description The function to generate initial values for MLE (single-source profile)
#' @param n_rep The number of replicates
#' @param hnc The assumed number of contributors
generate_init_deconvo_single <- function(n_rep){
  tph_pre <- rbeta(n_rep, 1, 1)
  deg <- rbeta(1, 1, 1)

  return(c(tph, deg))
}

#' generate_init_deconvo_mix
#'
#' @description The function to generate initial values for MLE (mixture profile)
#' @param n_rep The number of replicates
#' @param hnc The assumed number of contributors
generate_init_deconvo_mix <- function(n_rep, hnc){
  tph_pre <- rbeta(n_rep, 1, 1)
  mp <- rdirichlet(1, rep(1, hnc))
  mp <- mp[1, 1:(hnc - 1)]
  deg <- rbeta(hnc, 1, 1)

  return(c(tph_pre, mp, deg))
}

#' make_dt_pg_list
#'
#' @description The function to make a list of data.tables for the result of probabilistic genotyping
#' @param markers_use A character vector of the used markers
#' @param label_contrib A character vector of the contributor names
#' @param gt_perm_use Genotype permutations in each used markers
#' @param gtc_logprob_use A list of the genotype probabilities in each used marker
#' @param lpd_gtc_all_rep A list of the log likelihood of obtaining observed peak heights for each genotype combination in all markers and replicates
#' @param loglike_gtc_all_mk A list of the log likelihood combining genotype probabilities and observed peak heights
make_dt_pg_list <- function(markers_use, label_contrib, gt_perm_use, gtc_logprob_use, lpd_gtc_all_rep, loglike_gtc_all_mk){
  hnc <- length(label_contrib)
  n_mk <- length(gt_perm_use)
  n_rep <- length(lpd_gtc_all_rep)
  dt_pg_list <- vector(mode = "list", length = n_mk)

  for(i in 1:n_mk){

    # Paste two alleles of each contributor
    gtc_mk <- matrix(unlist(gt_perm_use[[i]]), ncol = 2 * hnc, byrow = TRUE)
    tmp_gtc_mk <- gtc_mk
    tmp_gtc_mk[which(tmp_gtc_mk == 99, arr.ind = TRUE)] <- "Q"
    gtc_mk_bind <- matrix("", nrow(tmp_gtc_mk), hnc)
    for(j in 1:hnc){
      gtc_mk_bind[, j] <- apply(tmp_gtc_mk[, c(2 * j - 1, 2 * j), drop = FALSE], 1, paste, collapse = ", ")
    }

    # Calculate posterior probabilities of each genotype combination
    loglike_gtc <- loglike_gtc_all_mk[[i]]
    tmp <- exp(loglike_gtc)
    post_prob_gtc <- tmp / sum(tmp)

    # Make a data.table
    mat_pg <- cbind(gtc_mk, gtc_mk_bind, gtc_logprob_use[[i]])
    for(j in 1:n_rep){
      mat_pg <- cbind(mat_pg, lpd_gtc_all_rep[[j]][[i]])
    }
    mat_pg <- cbind(mat_pg, loglike_gtc, post_prob_gtc)
    mat_pg <- mat_pg[post_prob_gtc != 0, , drop = FALSE]
    dt_pg <- as.data.table(mat_pg)
    names(dt_pg) <- c(paste0("Allele", 1:(2 * hnc)), label_contrib, "Logprob_gtc", paste0("Loglike_ph_rep", 1:n_rep), "Loglike", "Probability")

    # Change data type
    col_numeric <- c(paste0("Allele", 1:(2 * hnc)), "Logprob_gtc", paste0("Loglike_ph_rep", 1:n_rep), "Loglike", "Probability")
    dt_pg[, (col_numeric) := lapply(.SD, as.numeric), .SDcols = col_numeric]

    # Descending order for posterior probabilities
    dt_pg <- dt_pg[order(dt_pg[, Probability], decreasing = TRUE)]

    dt_pg_list[[i]] <- dt_pg
  }
  names(dt_pg_list) <- markers_use
  return(dt_pg_list)
}

#' perform_deconvo
#'
#' @description The function to perform mixture deconvolution
#' @param hnc The assumed number of contributors
#' @param n_rep The number of replicates
#' @param markers_use A vector of the autosomal markers used
#' @param ats_use A vector of the analytical thresholds for each marker
#' @param peaks_list A list of the peaks
#' @param sizes_list A list of the sizes
#' @param heights_list A list of the heights
#' @param min_size The minimum size in the selected kit
#' @param pop_al_list A list of the alleles in the population database
#' @param pop_prob_list A list of the allele probabilities in the population database
#' @param theta The theta value for sub-population effect
#' @param bsr_fltr The back stutter filter to narrow down a list of the genotype combinations
#' @param fsr_fltr The forward stutter filter to narrow down a list of the genotype combinations
#' @param dsr_fltr The double-back stutter filter to narrow down a list of the genotype combinations
#' @param m2sr_fltr The minus 2-nt stutter filter to narrow down a list of the genotype combinations
#' @param n_dropin_fltr The number of drop-in filter to narrow down a list of the genotype combinations
#' @param h_dropin_fltr The drop-in peak height filter to narrow down a list of the genotype combinations
#' @param freq_di The drop-in frequency
#' @param shape_di The shape parameter of the gamma distribution for the drop-in model
#' @param rate_di The rate parameter of the gamma distribution for the drop-in model
#' @param ae_use A vector of the locus-specific amplification efficiencies for each marker
#' @param cv_allele The coefficient of variation for allelels
#' @param cv_stutter The coefficient of variation for stutters
#' @param list_dt_average_sr A list of the data.table for the average stutter ratios
#' @param min_mp Minimum mixture proportion
#' @param min_deg Minimum degradation parameter
#' @param stop_iter Number of iterations without improvement before stopping one MLE run
#' @param max_iter Maximum number of iterations in one MLE run
#' @param num_converged Number of converged MLE runs before stopping
#' @param max_mle_runs Maximum number of MLE attempts
perform_deconvo <- function(hnc, n_rep, markers_use, ats_use, peaks_list, sizes_list, heights_list, min_size,
                            pop_al_list, pop_prob_list, theta,
                            bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr, n_dropin_fltr, h_dropin_fltr,
                            freq_di, shape_di, rate_di,
                            ae_use, cv_allele, cv_stutter, list_dt_average_sr,
                            min_mp, min_deg, stop_iter, max_iter, num_converged, max_mle_runs){

  # --- Start computation --- #

  start_time <- proc.time()

  # --- Define objects --- #

  # Peak
  n_mk <- length(markers_use)
  n_peak_use <- rep(NA, n_mk)
  names(n_peak_use) <- markers_use
  median_sizes_use <- rep(NA, n_mk)
  tph_per_mk_use <- matrix(0, n_mk, n_rep)

  # Genotype combinations
  n_gt_comb_use <- rep(NA, n_mk)
  n_gt_perm_use <- rep(0, n_mk)
  gt_perm_use <- vector(mode = "list", length = n_mk)
  names(gt_perm_use) <- markers_use

  # Genotype probabilities
  gtc_logprob_use <- vector(mode = "list", length = n_mk)
  names(gtc_logprob_use) <- markers_use

  # Relative peak heights
  rph_base_use <- vector(mode = "list", length = n_mk)
  names(rph_base_use) <- markers_use

  # CV squared
  cv_square_base_use <- vector(mode = "list", length = n_mk)
  names(cv_square_base_use) <- markers_use

  # Scale of degradation corrected by sizes
  deg_size_base_use <- vector(mode = "list", length = n_mk)
  names(deg_size_base_use) <- markers_use

  # Observed peak heights
  oph_base_use <- vector(mode = "list", length = n_rep)
  for(i in 1:n_rep){
    oph_base_use[[i]] <- vector(mode = "list", length = n_mk)
  }

  # Drop-in
  logprob_di_use <- rep(NA, n_mk)
  names(logprob_di_use) <- markers_use
  loglike_di_list <- vector(mode = "list", length = n_rep)
  for(i in 1:n_rep){
    loglike_di_list[[i]] <- vector(mode = "list", length = n_mk)
  }

  # De-duplicated objects
  n_general_list <- vector(mode = "list", length = n_rep)
  n_dropout_list <- vector(mode = "list", length = n_rep)
  n_dropin_list <- vector(mode = "list", length = n_rep)
  n_ignore_list <- vector(mode = "list", length = n_rep)
  index_dedup1_list_cpp <- vector(mode = "list", length = n_rep)
  index_dedup2_list_cpp <- vector(mode = "list", length = n_rep)
  for(i in 1:n_rep){
    n_general_list[[i]] <- rep(NA, n_mk)
    n_dropout_list[[i]] <- rep(NA, n_mk)
    n_dropin_list[[i]] <- rep(NA, n_mk)
    n_ignore_list[[i]] <- rep(NA, n_mk)
    index_dedup1_list_cpp[[i]] <- vector(mode = "list", length = n_mk)
    index_dedup2_list_cpp[[i]] <- vector(mode = "list", length = n_mk)
  }

  for(i in 1:n_mk){

    # --- Extract data for one marker --- #

    # Extract peak data
    mk <- markers_use[i]
    peaks <- peaks_list[[mk]]
    sizes <- sizes_list[[mk]]
    heights_all_rep <- heights_list[[mk]]
    n_peak_use[i] <- n_peak <- length(peaks)

    # Extract parameters
    at <- ats_use[mk]
    ae <- ae_use[mk]
    dt_ave_sr <- list_dt_average_sr[[mk]]
    sr_alleles_mk <- as.numeric(dt_ave_sr$Allele)
    bsrs_mk <- dt_ave_sr[["Average BSR"]]
    fsrs_mk <- dt_ave_sr[["Average FSR"]]
    dsrs_mk <- dt_ave_sr[["Average DSR"]]
    m2srs_mk <- dt_ave_sr[["Average M2SR"]]
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

    # Extract allele probabilities
    pop_al <- pop_al_list[[mk]]
    pop_prob <- pop_prob_list[[mk]]

    # --- Calculate the median size in the current marker --- #

    median_sizes_use[i] <- median(sizes)

    # --- Calculate total peak heights in each replicate --- #

    tph_per_mk_use[i, ] <- sapply(heights_all_rep, sum) / 2

    # --- Calculate the probability of drop-in in the current marker --- #

    auc_below_at <- pgamma(q = at, shape = shape_di, rate = rate_di)
    logprob_di_use[i] <- log(freq_di / (1 - auc_below_at))

    # --- Prepare the relative peak heights --- #

    rph_base_use[[i]] <- rph_base_mk <- make_rph_base(peaks, bsrs, fsrs, dsrs, m2srs, ae)

    # --- Prepare the CV squared values for each peak --- #

    cv_square_base_use[[i]] <- cv_square_base_mk <- c(rep(cv_allele^2, n_peak), rep(cv_stutter^2, 4 * n_peak), 0, -1)

    # --- Prepare the scale of degradation corrected by sizes --- #

    deg_size_base_use[[i]] <- c(rep((sizes - min_size) / min_size, 5), 0, -1)

    # --- Generate all possible genotype combinations --- #

    # All candidate genotypes of one contributor
    cand_gt <- combinations(n = n_peak, r = 2, v = peaks, repeats.allowed = TRUE)
    cand_gt <- asplit(cand_gt, 1)

    # Indices of genotypes of each contributor in all candidate genotype combinations
    index_gt_comb <- combinations(n = length(cand_gt), r = hnc, repeats.allowed = TRUE)
    index_gt_comb_cpp <- index_gt_comb - 1 # C++ index
    index_gt_comb_cpp <- asplit(index_gt_comb_cpp, 1)

    # Indices of genotypes of each contributor in each adopted genotype combination
    tmp <- judge_gt_comb_mk(cand_gt, index_gt_comb_cpp, peaks, heights_all_rep,
                            bsrs, fsrs, dsrs, m2srs,
                            bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr, n_dropin_fltr, h_dropin_fltr)
    index_gt_comb_adopt_cpp <- tmp[[1]]
    bool_gt_comb_di_mk_rep <- tmp[[2]]

    # The number of adopted genotype combinations
    n_gt_comb_use[i] <- length(index_gt_comb_adopt_cpp)

    if(!any(is.element(n_gt_comb_use, 0))){

      # --- Calculate genotype probabilities in each genotype combination --- #

      gt_comb <- extract_gts_adopt(cand_gt, index_gt_comb_adopt_cpp)
      n_obs <- make_n_obs(pop_al, prof_ref = numeric(0))
      obs_al <- as.numeric(names(n_obs))
      gt_comb_prob_mk <- calc_gts_prob_mk(gt_comb, pop_al, pop_prob, theta, prof_ref = numeric(0), hyp_id = rep(0, hnc), n_obs, obs_al)
      log_gt_comb_prob_mk <- log(gt_comb_prob_mk)

      # --- Create objects for genotype permutations --- #

      tmp <- create_gt_perm_obj(cand_gt, index_gt_comb_adopt_cpp, log_gt_comb_prob_mk, bool_gt_comb_di_mk_rep)
      gt_perm_use[[i]] <- gt_perm <- tmp$gt_perm
      n_gt_perm_use[i] <- n_gt_perm <- length(gt_perm)
      gtc_logprob_use[[i]] <- tmp$log_gt_perm_prob_mk
      bool_gt_perm_di_mk_rep <- tmp$bool_gt_perm_di_mk_rep

      for(j in 1:n_rep){

        # --- Prepare the observed peak heights --- #

        oph_base_use[[j]][[i]] <- c(rep(heights_all_rep[[j]], 5), 0, -1)

        # --- Prepare objects for de-duplication (pattern 1) --- #

        # Search the indices of the relative peak heights for each genotype combinations
        index_rph_gtc_cpp <- connect_rph_gtc_mk(gt_perm, peaks, bsrs, fsrs, dsrs, m2srs, bool_gt_perm_di_mk_rep[[j]])
        index_rph_gtc_cpp <- unlist(index_rph_gtc_cpp, recursive = FALSE)
        index_rph_gtc <- matrix(unlist(index_rph_gtc_cpp) + 1, nrow = 2 * hnc) # R index

        # Delete duplicated indices
        bool_duplicated_adopt <- !duplicated(t(index_rph_gtc))
        index_dedup1_tmp <- index_rph_gtc[, bool_duplicated_adopt, drop = FALSE]

        # Observe peak heights for de-duplicated objects (pattern 1)
        heights <- heights_all_rep[[j]]
        oph_dedup1_tmp <- rep(heights, n_gt_perm)
        oph_dedup1_tmp <- oph_dedup1_tmp[bool_duplicated_adopt]

        # Calculate relative peak heights in each peak for de-duplicated objects
        rph_dedup1_tmp <- matrix(rph_base_mk[index_dedup1_tmp], nrow = 2 * hnc)
        rph_dedup1_tmp <- apply(rph_dedup1_tmp, 2, sum)

        # Indices of each peak type
        index_oph0 <- which(oph_dedup1_tmp == 0)
        index_ignore <- c(which(rph_dedup1_tmp == 0), intersect(which(rph_dedup1_tmp < 0), index_oph0)) # {rph = 0: allele 99 without assuming drop-out alleles} or {rph < 0 and oph = 0: drop-in peaks in other replicates}
        index_dropin <- setdiff(which(rph_dedup1_tmp < 0), index_ignore)
        index_dropout <- setdiff(index_oph0, c(index_ignore, index_dropin))
        index_general <- setdiff(1:length(rph_dedup1_tmp), c(index_ignore, index_dropin, index_dropout))

        # The number of indices for each peak type
        n_general_list[[j]][i] <- length(index_general)
        n_dropout_list[[j]][i] <- length(index_dropout)
        n_dropin_list[[j]][i] <- length(index_dropin)
        n_ignore_list[[j]][i] <- length(index_ignore)

        # Assign C++ indices of the relative peak heights for de-duplicated objects (pattern 1)
        index_dedup1_list_cpp[[j]][[i]] <- c(asplit(index_dedup1_tmp[, index_general, drop = FALSE] - 1, 2),
                                             asplit(index_dedup1_tmp[, index_dropout, drop = FALSE] - 1, 2),
                                             asplit(index_dedup1_tmp[, index_dropin, drop = FALSE] - 1, 2),
                                             asplit(index_dedup1_tmp[, index_ignore, drop = FALSE] - 1, 2))

        # Assign log-likelihood of the drop-in for de-duplicated objects (pattern 1)
        if(length(index_dropin) != 0){
          loglike_di_list[[j]][[i]] <- log(dgamma(oph_dedup1_tmp[index_dropin], shape_di, rate_di)) + logprob_di_use[i]
        }else{
          loglike_di_list[[j]][[i]] <- numeric(0)
        }

        # --- Prepare objects for de-duplication (pattern 2) --- #

        index_dedup2_list_cpp[[j]][[i]] <- create_index_dedup2_mk(index_rph_gtc_cpp, index_dedup1_list_cpp[[j]][[i]])
      }
    }else{
      n_gt_perm_use[i] <- 0
    }
  }

  if(!any(is.element(n_gt_perm_use, 0))){

    # --- Determine the lower and upper values for the parameter tph --- #

    median_sizes_use <- rep(median_sizes_use, n_rep)
    tph_per_mk_use <- as.vector(tph_per_mk_use)
    fit_exp <- lm(log(tph_per_mk_use) ~ median_sizes_use)
    mu_log <- coef(fit_exp)[1] + coef(fit_exp)[2] * min_size
    sigma_log <- sd(residuals(fit_exp))
    tph_lower <- exp(mu_log - 3 * sigma_log)
    tph_upper <- exp(mu_log + 3 * sigma_log)

    # --- Bayesian optimization --- #

    if(hnc == 1){

      # likelihood function for single-source profiles
      like_func_single <- function(par){
        # Prepare parameters
        tph <- par[1:n_rep] * (tph_upper - tph_lower) + tph_lower
        deg <- par[n_rep + 1] * (1 - min_deg) + min_deg
        mp <- 1

        if(any(deg == 0)){
          loglike <- Inf
        }else{
          # Calculate log probability densities of all genotype combinations
          lpd_gtc_all_rep <- calc_lpd_gtc_all_rep(tph, mp, deg, hnc, ats_use, n_general_list, n_dropout_list, n_dropin_list, n_ignore_list,
                                                  rph_base_use, cv_square_base_use, deg_size_base_use, oph_base_use, loglike_di_list,
                                                  index_dedup1_list_cpp, index_dedup2_list_cpp, n_gt_perm_use, n_peak_use)

          # Calculate log likelihood in each marker
          tmp <- calc_loglike_cpp(lpd_gtc_all_rep, gtc_logprob_use)
          loglike <- - sum(tmp[[1]])

          if(is.na(loglike) || is.nan(loglike) || is.infinite(loglike)){
            loglike <- Inf
          }
        }

        return(loglike)
      }

      lower_bounds <- c(rep(0, n_rep), 0)
      upper_bounds <- c(rep(1, n_rep), 1)

      # Perform simulated annealing
      current_opt <- NULL
      current_score <- Inf
      current_score_round <- Inf
      count_mle <- 0
      count_solution <- 0
      init_list <- list()
      opt_list <- list()
      repeat{
        count_mle <- count_mle + 1

        init_list[[count_mle]] <- init <- generate_init_deconvo_single(n_rep)

        opt_list[[count_mle]] <- opt <- GenSA(
          par = init,
          fn = like_func_single,
          lower = lower_bounds,
          upper = upper_bounds,
          control = list(verbose = TRUE, nb.stop.improvement = stop_iter, maxit = max_iter)
        )

        score <- opt$value
        if(score < current_score){
          current_opt <- opt
          current_score <- score
        }

        score_round <- round(score, 3)
        if(score_round < current_score_round){
          current_score_round <- score_round
          count_solution <- 1
        }else if(score_round == current_score_round){
          count_solution <- count_solution + 1
        }

        if(count_solution == num_converged){
          break
        }else if(count_mle == max_mle_runs){
          break
        }
      }

      # Extract estimated parameters
      est_params <- current_opt$par
      tph <- est_params[1:n_rep] * (tph_upper - tph_lower) + tph_lower
      mp <- 1
      deg <- est_params[n_rep + 1] * (1 - min_deg) + min_deg

    }else{

      # likelihood function for mixture profiles
      like_func_mix <- function(par){
        # Prepare parameters
        tph <- par[1:n_rep] * (tph_upper - tph_lower) + tph_lower
        mp_pre <- par[(n_rep + 1):(n_rep + hnc - 1)] * (1 - 2 * min_mp) + min_mp
        deg <- par[(n_rep + hnc):(n_rep + 2 * hnc - 1)] * (1 - min_deg) + min_deg

        mp <- c(mp_pre, 1 - sum(mp_pre))
        if(any(mp == 0)){
          loglike <- Inf
        }else if(any(deg == 0)){
          loglike <- Inf
        }else{

          # Calculate log probability densities of all genotype combinations
          lpd_gtc_all_rep <- calc_lpd_gtc_all_rep(tph, mp, deg, hnc, ats_use, n_general_list, n_dropout_list, n_dropin_list, n_ignore_list,
                                                  rph_base_use, cv_square_base_use, deg_size_base_use, oph_base_use, loglike_di_list,
                                                  index_dedup1_list_cpp, index_dedup2_list_cpp, n_gt_perm_use, n_peak_use)

          # Calculate log likelihood in each marker
          tmp <- calc_loglike_cpp(lpd_gtc_all_rep, gtc_logprob_use)
          loglike <- - sum(tmp[[1]])

          if(is.na(loglike) || is.nan(loglike) || is.infinite(loglike)){
            loglike <- Inf
          }
        }

        return(loglike)
      }

      lower_bounds <- c(rep(0, n_rep), rep(0, hnc - 1), rep(0, hnc))
      upper_bounds <- c(rep(1, n_rep), rep(1, hnc - 1), rep(1, hnc))

      # Perform simulated annealing
      current_opt <- NULL
      current_score <- Inf
      current_score_round <- Inf
      count_mle <- 0
      count_solution <- 0
      init_list <- list()
      opt_list <- list()
      repeat{
        count_mle <- count_mle + 1

        init_list[[count_mle]] <- init <- generate_init_deconvo_mix(n_rep, hnc)

        opt_list[[count_mle]] <- opt <- GenSA(
          par = init,
          fn = like_func_mix,
          lower = lower_bounds,
          upper = upper_bounds,
          control = list(verbose = TRUE, nb.stop.improvement = stop_iter, maxit = max_iter)
        )

        score <- opt$value
        if(score < current_score){
          current_opt <- opt
          current_score <- score
        }

        score_round <- round(score, 3)
        if(score_round < current_score_round){
          current_score_round <- score_round
          count_solution <- 1
        }else if(score_round == current_score_round){
          count_solution <- count_solution + 1
        }

        if(count_solution == num_converged){
          break
        }else if(count_mle == max_mle_runs){
          break
        }
      }

      # Extract estimated parameters
      est_params <- current_opt$par
      tph <- est_params[1:n_rep] * (tph_upper - tph_lower) + tph_lower
      mp_pre <- est_params[(n_rep + 1):(n_rep + hnc - 1)] * (1 - 2 * min_mp) + min_mp
      mp <- c(mp_pre, 1 - sum(mp_pre))
      deg <- est_params[(n_rep + hnc):(n_rep + 2 * hnc - 1)] * (1 - min_deg) + min_deg

      # Descending order for estimated mixture proportions
      index_order <- order(mp, decreasing = TRUE)
      mp <- mp[index_order]
      deg <- deg[index_order]
    }

    if(count_solution == num_converged){
      warning_message <- NULL
    }else{
      warning_message <- paste0("- A total of ", count_mle, " MLE runs were performed, of which ", count_solution, " converged under ", hnc, " contributor. This number is below the specified threshold of ", num_converged, " converged MLE runs , which may compromise the reliability of the results.<br>")
    }

    # --- Calculate the likelihood using the estimated parameters --- #

    # Calculate log probability densities of all genotype combinations
    lpd_gtc_all_rep <- calc_lpd_gtc_all_rep(tph, mp, deg, hnc, ats_use, n_general_list, n_dropout_list, n_dropin_list, n_ignore_list,
                                            rph_base_use, cv_square_base_use, deg_size_base_use, oph_base_use, loglike_di_list,
                                            index_dedup1_list_cpp, index_dedup2_list_cpp, n_gt_perm_use, n_peak_use)

    # Calculate log likelihood in each marker
    tmp <- calc_loglike_cpp(lpd_gtc_all_rep, gtc_logprob_use)
    loglike_all_mk <- tmp[[1]]
    loglike_gtc_all_mk <- tmp[[2]]

    dt_pg_list <- make_dt_pg_list(markers_use, paste0("U", 1:hnc), gt_perm_use, gtc_logprob_use, lpd_gtc_all_rep, loglike_gtc_all_mk)

    log10like_all_mk <- loglike_all_mk / log(10)
    dt_log10like <- data.table(Marker = c(markers_use, "Overall"),
                               Log10_Likelihood = c(log10like_all_mk, sum(log10like_all_mk)))

    dt_param <- data.table("Parameter" = c(rep("Template amount", n_rep), rep("Mixture proportion", hnc), rep("Degradation", hnc)),
                           "Contributor" = c(rep(NA, n_rep), paste0("U", 1:hnc), paste0("U", 1:hnc)),
                           "Replicate" = c(1:n_rep, rep(NA, hnc), rep(NA, hnc)),
                           "Estimated value" = c(tph, mp, deg))

    # End computation
    end_time <- proc.time()
    comp_time <- end_time[3] - start_time[3]

    return(list("dt_pg_list" = dt_pg_list,
                "dt_log10like" = dt_log10like,
                "dt_param" = dt_param,
                "logprob_di_use" = logprob_di_use,
                "comp_time" = comp_time,
                "init_list" = init_list,
                "opt_list" = opt_list,
                "warning_message" = warning_message,
                "error_message" = NULL))
  }else{
    markers_unexplainable <- markers_use[n_gt_perm_use == 0]
    if(hnc == 1){
      error_message <- paste0("The following marker(s) cannot be explained under ", hnc, " contributor.<br>")
    }else{
      error_message <- paste0("The following marker(s) cannot be explained under ", hnc, " contributors.<br>")
    }
    for(i in 1:length(markers_unexplainable)){
      error_message <- paste0(error_message, "  - ", markers_unexplainable[i], "<br>")
    }

    return(list("dt_pg_list" = NULL,
                "dt_log10like" = NULL,
                "dt_param" = NULL,
                "logprob_di_use" = NULL,
                "comp_time" = NULL,
                "init_list" = NULL,
                "opt_list" = NULL,
                "warning_message" = NULL,
                "error_message" = error_message))
  }
}

#' perform_deconvo_all_noc
#'
#' @description The function to perform mixture deconvolution for all assumed number of contributors
#' @param dt_csp The data.table of the crime stain profiles
#' @param dt_af The data.table of the allele frequencies
#' @param dt_kit_auto The data.table of the information on the autosomal markers in a kit
#' @param dt_ignore_mk The data.table of the information on the ignored markers
#' @param dt_anameth The data.table of the analysis method
#' @param list_dt_average_sr A list of the data.table for the average stutter ratios
#' @param hnc_all A vector of the assumed number of contributors
#' @param theta The theta value for sub-population effect
#' @param min_mp Minimum mixture proportion
#' @param min_deg Minimum degradation parameter
#' @param stop_iter Number of iterations without improvement before stopping one MLE run
#' @param max_iter Maximum number of iterations in one MLE run
#' @param num_converged Number of converged MLE runs before stopping
#' @param max_mle_runs Maximum number of MLE attempts
perform_deconvo_all_noc <- function(dt_csp, dt_af, dt_kit_auto, dt_ignore_mk, dt_anameth, list_dt_average_sr, hnc_all, theta,
                                    min_mp, min_deg, stop_iter, max_iter, num_converged, max_mle_runs){
  markers_use <- dt_ignore_mk$Marker[!dt_ignore_mk$Ignore]
  n_mk <- length(markers_use)
  if(n_mk != 0){
    n_rep <- length(unique(dt_csp$SampleName))

    # Extract ATs of the used markers
    markers_at <- dt_anameth$Subcategory[dt_anameth$Category == "AT"]
    ats_auto <- dt_anameth$Value[dt_anameth$Category == "AT"]
    ats_use <- as.numeric(ats_auto[match(markers_use, markers_at)])
    names(ats_use) <- markers_use
    rm(markers_at, ats_auto)

    # Extract peak data from crime stain profiles
    tmp <- make_peak_data(dt_csp, dt_kit_auto, markers_use, ats_use)
    peaks_list <- tmp$peaks_list
    sizes_list <- tmp$sizes_list
    heights_list <- tmp$heights_list
    min_size <- min(dt_kit_auto$Size)

    # Set allele probabilities
    tmp <- set_al_prob(markers_use, dt_af, peaks_list)
    pop_al_list <- tmp$pop_al_list
    pop_prob_list <- tmp$pop_prob_list

    # Extract filters to remove unrealistic genotype combinations
    bsr_fltr <- dt_anameth$Value[dt_anameth$Subcategory == "Maximum BSR"]
    fsr_fltr <- dt_anameth$Value[dt_anameth$Subcategory == "Maximum FSR"]
    dsr_fltr <- dt_anameth$Value[dt_anameth$Subcategory == "Maximum DSR"]
    m2sr_fltr <- dt_anameth$Value[dt_anameth$Subcategory == "Maximum M2SR"]
    n_dropin_fltr <- dt_anameth$Value[dt_anameth$Subcategory == "Maximum number of drop-in"]
    h_dropin_fltr <- dt_anameth$Value[dt_anameth$Subcategory == "Maximum height of drop-in"]

    # Extract expected LSAE
    expect_ae <- dt_anameth$Value[dt_anameth$Category == "LSAE"]
    names(expect_ae) <- dt_anameth$Subcategory[dt_anameth$Category == "LSAE"]
    ae_use <- expect_ae[match(markers_use, names(expect_ae))]

    # Extract peak height parameters
    cv_allele <- dt_anameth$Value[dt_anameth$Subcategory == "CV for allele"]
    cv_stutter <- dt_anameth$Value[dt_anameth$Subcategory == "CV for stutter"]

    # Extract drop-in parameters
    freq_di <- dt_anameth$Value[dt_anameth$Subcategory == "Frequency"]
    shape_di <- dt_anameth$Value[dt_anameth$Subcategory == "Shape parameter"]
    rate_di <- dt_anameth$Value[dt_anameth$Subcategory == "Rate parameter"]

    n_hnc <- length(hnc_all)
    data_hnc_list <- vector(mode = "list", length = n_hnc)
    for(i in 1:n_hnc){
      hnc <- hnc_all[i]
      data_hnc <- try(perform_deconvo(hnc, n_rep, markers_use, ats_use, peaks_list, sizes_list, heights_list, min_size,
                                      pop_al_list, pop_prob_list, theta,
                                      bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr, n_dropin_fltr, h_dropin_fltr,
                                      freq_di, shape_di, rate_di,
                                      ae_use, cv_allele, cv_stutter, list_dt_average_sr,
                                      min_mp, min_deg, stop_iter, max_iter, num_converged, max_mle_runs))
      if(!inherits(data_hnc, "try-error")){
        data_hnc_list[[i]] <- data_hnc
      }else{
        data_hnc_list[[i]] <- list("dt_pg_list" = NULL,
                                   "dt_log10like" = NULL,
                                   "dt_param" = NULL,
                                   "logprob_di_use" = NULL,
                                   "comp_time" = NULL,
                                   "init_list" = NULL,
                                   "opt_list" = NULL,
                                   "warning_message" = NULL,
                                   "error_message" = "Unexpected error")
      }
    }
  }else{
    peaks_list <- NULL
    heights_list <- NULL
    pop_al_list <- NULL
    pop_prob_list <- NULL
    ae_use <- NULL
    min_size <- NULL
    data_hnc_list <- NULL
  }

  return(list("hnc_all" = hnc_all,
              "markers_use" = markers_use,
              "peaks_list" = peaks_list,
              "heights_list" = heights_list,
              "pop_al_list" = pop_al_list,
              "pop_prob_list" = pop_prob_list,
              "ae_use" = ae_use,
              "min_size" = min_size,
              "data_hnc_list" = data_hnc_list))
}
