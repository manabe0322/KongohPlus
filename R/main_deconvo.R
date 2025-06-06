#' create_index_dedup2_mk
#'
#' @description The function to create 'index_dedup2_mk'
#' @param index_rph_gtc_cpp C++ indices of the relative peak heights for each genotype combinations
#' @param index_dedup1_rep1_mk1 C++ indices of the relative peak heights for de-duplicated objects (pattern 1)
create_index_dedup2_mk <- function(index_rph_gtc_cpp, index_dedup1_rep1_mk1){
  paste_index_rph_gtc_cpp <- sapply(index_rph_gtc_cpp, paste, collapse = "-")
  paste_index_dedup1_rep1_mk1 <- sapply(index_dedup1_rep1_mk1, paste, collapse = "-")
  return(match(paste_index_rph_gtc_cpp, paste_index_dedup1_rep1_mk1))
}

#' softmax
#'
#' @description The softmax function
#' @param x A numeric vector
softmax <- function(x){
  exp_x <- exp(x - max(x))
  return(exp_x / sum(exp_x))
}

make_bounds_list <- function(n_rep, hnc, max_tph){
  bounds <- list()

  # tph1, tph2, ...
  for(i in 1:n_rep){
    bounds[[paste0("tph", i)]] <- c(1, max_tph)
  }

  # lambda1
  bounds[["lambda1"]] <- c(-10, 10)

  # delta2, delta3, ...
  if(hnc > 1){
    for(i in 2:hnc){
      bounds[[paste0("delta", i)]] <- c(-5, 5)
    }
  }

  # deg1, deg2, ...
  for(i in 1:hnc){
    bounds[[paste0("deg", i)]] <- c(0, 1)
  }

  return(bounds)
}

#' generate_init_sobol_1p
#'
#' @description The function to generate initial values based on Sobol sequence for one contributor
#' @param x Sobol point
#' @param n_rep The number of replicates
#' @param tph_mean The mean of the total peak height
#' @param tph_sd The SD of the total peak height
generate_init_sobol_1p <- function(x, n_rep, tph_mean, tph_sd){
  # Dimensions needed in total: tph + deg
  dim_total <- n_rep + 1

  # Get only one Sobol point (x-th)
  sobol_vec <- sobol(n = x, d = dim_total, randomize = "digital.shift")[x, ]

  # Index partitioning
  tph_raw <- sobol_vec[1:n_rep]
  deg_raw <- sobol_vec[n_rep + 1]

  # Scaling
  tph <- qnorm(tph_raw, mean = tph_mean, sd = tph_sd)
  deg <- qbeta(deg_raw, shape1 = 1, shape2 = 1)

  return(list(list(tph = tph, deg = deg)))
}

#' run_mle_1p
#'
#' @description The function to perform MLE for one contributor
#' @param model_name The name of the stan model for deconvolution
#' @param data_list_main The input data for stan
#' @param n_rep The number of replicates
#' @param tph_mean The mean of the total peak height
#' @param tph_sd The SD of the total peak height
#' @param num_init The number of initial points for Pathfinder
run_mle_1p <- function(model_name, data_list_main, n_rep, tph_mean, tph_sd, num_init){
  warn_optim <- FALSE

  init_list <- lapply(1:num_init, function(x) generate_init_sobol_1p(x, n_rep, tph_mean, tph_sd))

  tmp_param <- matrix(NA, num_init, n_rep + 2)
  model_main <- stan_package_model(name = model_name, package = "KongohPlus")
  for(i in 1:num_init){
    fit <- suppressWarnings(model_main$optimize(data_list_main, init = init_list[[i]], refresh = 0, show_messages = FALSE, show_exceptions = FALSE))
    tmp <- try(fit$draws())
    if(!inherits(tmp, "try-error")){
      tmp_param[i, ] <- tmp[1, ]
    }
  }

  index_max_lp <- which.max(tmp_param[, 1])
  if(length(index_max_lp) == 1){
    est_param <- tmp_param[index_max_lp, -1]
  }else{
    warn_optim <- TRUE

    tmp_param <- matrix(NA, num_init, n_rep + 3)
    for(i in 1:num_init){
      fit <- model_main$pathfinder(data_list_main, init = init_list[[i]])
      tmp <- try(fit$draws())
      if(!inherits(tmp, "try-error")){
        tmp_param[i, ] <- tmp[which.max(tmp[, 2]), ]
      }
    }

    index_max_lp <- which.max(tmp_param[, 2])
    if(length(index_max_lp) == 1){
      est_param <- tmp_param[index_max_lp, c(-1, -2)]
    }else{
      est_param <- rep(NA, n_rep + 1)
    }
  }

  names(est_param) <- c(paste0("tph[", 1:n_rep, "]"), "deg[1]")
  return(list(est_param = est_param, warn_optim = warn_optim, init_list = init_list))
}

#' generate_init_sobol
#'
#' @description The function to generate initial values based on Sobol sequence
#' @param x Sobol point
#' @param n_rep The number of replicates
#' @param hnc The hypothesized number of contributors
#' @param tph_mean The mean of the total peak height
#' @param tph_sd The SD of the total peak height
generate_init_sobol <- function(x, n_rep, hnc, tph_mean, tph_sd){
  # Dimensions needed in total: tph + lambda + deg
  dim_total <- n_rep + hnc + hnc

  # Get only one Sobol point (x-th)
  sobol_vec <- sobol(n = x, d = dim_total, randomize = "digital.shift")[x, ]

  # Index partitioning
  tph_raw <- sobol_vec[1:n_rep]
  lambda_raw <- sobol_vec[(n_rep + 1):(n_rep + hnc)]
  deg_raw <- sobol_vec[(n_rep + hnc + 1):(dim_total)]

  # Scaling
  tph <- qnorm(tph_raw, mean = tph_mean, sd = tph_sd)
  lambda <- sort(qnorm(lambda_raw, mean = 0, sd = 1))
  deg <- qbeta(deg_raw, shape1 = 1, shape2 = 1)

  return(list(list(tph = tph, lambda = lambda, deg = deg)))
}

#' run_mle
#'
#' @description The function to perform MLE
#' @param model_name The name of the stan model for deconvolution
#' @param data_list_main The input data for stan
#' @param n_rep The number of replicates
#' @param hnc The hypothesized number of contributors
#' @param tph_mean The mean of the total peak height
#' @param tph_sd The SD of the total peak height
#' @param num_init The number of initial points for Pathfinder
run_mle <- function(model_name, data_list_main, n_rep, hnc, tph_mean, tph_sd, num_init){
  warn_optim <- FALSE

  init_list <- lapply(1:num_init, function(x) generate_init_sobol(x, n_rep, hnc, tph_mean, tph_sd))

  if(hnc == 2){
    tmp_param <- matrix(NA, num_init, 1 + n_rep + 3 * hnc)
    model_main <- stan_package_model(name = model_name, package = "KongohPlus")

    for(i in 1:num_init){
      # Create output_dir
      outdir <- file.path(tempdir(), paste0("run_", i))
      dir.create(outdir, showWarnings = FALSE)

      fit <- suppressWarnings(model_main$optimize(data_list_main, init = init_list[[i]], output_dir = outdir, refresh = 0, show_messages = FALSE, show_exceptions = FALSE))
      tmp <- try(fit$draws())
      if(!inherits(tmp, "try-error")){
        tmp_param[i, ] <- tmp[1, ]
      }
    }

    index_max_lp <- which.max(tmp_param[, 1])
    if(length(index_max_lp) == 1){
      est_param <- tmp_param[index_max_lp, -1]
    }else{
      warn_optim <- TRUE

      tmp_param <- matrix(NA, num_init, 2 + n_rep + 3 * hnc)
      for(i in 1:num_init){
        # Create output_dir
        outdir <- file.path(tempdir(), paste0("run_", i))
        dir.create(outdir, showWarnings = FALSE)

        fit <- model_main$pathfinder(data_list_main, init = init_list[[i]], output_dir = outdir)
        tmp <- try(fit$draws())
        if(!inherits(tmp, "try-error")){
          tmp_param[i, ] <- tmp[which.max(tmp[, 2]), ]
        }
      }

      index_max_lp <- which.max(tmp_param[, 2])
      if(length(index_max_lp) == 1){
        est_param <- tmp_param[index_max_lp, c(-1, -2)]
      }else{
        est_param <- rep(NA, n_rep + 3 * hnc)
      }
    }

  }else{
    cores <- min(4, detectCores() - 1)
    cl <- makeCluster(cores)
    clusterExport(cl, varlist = c("model_name", "data_list_main", "init_list", "n_rep", "hnc"), envir = environment())
    clusterEvalQ(cl, library(cmdstanr))
    clusterEvalQ(cl, library(instantiate))

    tmp_param <- parLapply(cl, 1:num_init, function(i){
      model_main <- stan_package_model(name = model_name, package = "KongohPlus")

      # Create output_dir
      outdir <- file.path(tempdir(), paste0("run_", i))
      dir.create(outdir, showWarnings = FALSE)

      fit <- suppressWarnings(model_main$optimize(data_list_main, init = init_list[[i]], output_dir = outdir, refresh = 0, show_messages = FALSE, show_exceptions = FALSE))
      tmp <- try(fit$draws())
      if(!inherits(tmp, "try-error")){
        return(as.vector(tmp[1, ]))
      }else{
        return(rep(NA, 1 + n_rep + 3 * hnc))
      }
    })
    stopCluster(cl)
    tmp_param <- matrix(unlist(tmp_param), ncol = 1 + n_rep + 3 * hnc, byrow = TRUE)

    index_max_lp <- which.max(tmp_param[, 1])
    if(length(index_max_lp) == 1){
      est_param <- tmp_param[index_max_lp, -1]
    }else{
      warn_optim <- TRUE

      cores <- min(4, detectCores() - 1)
      cl <- makeCluster(cores)
      clusterExport(cl, varlist = c("model_name", "data_list_main", "init_list", "n_rep", "hnc"), envir = environment())
      clusterEvalQ(cl, library(cmdstanr))
      clusterEvalQ(cl, library(instantiate))

      tmp_param <- parLapply(cl, 1:num_init, function(i){
        model_main <- stan_package_model(name = model_name, package = "KongohPlus")

        # Create output_dir
        outdir <- file.path(tempdir(), paste0("run_", i))
        dir.create(outdir, showWarnings = FALSE)

        fit <- model_main$pathfinder(data_list_main, init = init_list[[i]], output_dir = outdir)
        tmp <- try(fit$draws())
        if(!inherits(tmp, "try-error")){
          return(as.vector(tmp[which.max(tmp[, 2]), ]))
        }else{
          return(rep(NA, 2 + n_rep + 3 * hnc))
        }
      })
      stopCluster(cl)
      tmp_param <- matrix(unlist(tmp_param), ncol = 2 + n_rep + 3 * hnc, byrow = TRUE)

      index_max_lp <- which.max(tmp_param[, 2])
      if(length(index_max_lp) == 1){
        est_param <- tmp_param[index_max_lp, c(-1, -2)]
      }else{
        est_param <- rep(NA, n_rep + 3 * hnc)
      }
    }
  }

  names(est_param) <- c(paste0("tph[", 1:n_rep, "]"), paste0("lambda[", 1:hnc, "]"), paste0("deg[", 1:hnc, "]"), paste0("mp[", 1:hnc, "]"))
  return(list(est_param = est_param, warn_optim = warn_optim, init_list = init_list))
}

#' make_dt_pg_list
#'
#' @description The function to make a list of data.tables for the result of probabilistic genotyping
#' @param markers_use A character vector of the used markers
#' @param label_contrib A character vector of the contributor names
#' @param gtc_use Genotype combinations in the used markers
#' @param gtc_logprob_use A list of the genotype probabilities in each used marker
#' @param lpd_gtc_use A list of the log likelihood of each genotype in all markers and replicates
make_dt_pg_list <- function(markers_use, label_contrib, gtc_use, gtc_logprob_use, lpd_gtc_use){
  hnc <- length(label_contrib)
  n_mk <- length(gtc_use)
  n_rep <- length(lpd_gtc_use)
  dt_pg_list <- vector(mode = "list", length = n_mk)

  for(i in 1:n_mk){

    # Reverse order of the contributor (minor to major -> major to minor)
    gtc_mk <- matrix(unlist(gtc_use[[i]]), ncol = 2 * hnc, byrow = TRUE)
    gtc_mk_reorder <- gtc_mk * 0
    for(j in 1:hnc){
      gtc_mk_reorder[, c(2 * j - 1, 2 * j)] <- gtc_mk[, c(2 * (hnc - j) + 1, 2 * (hnc - j) + 2)]
    }

    # Paste two alleles of each contributor
    tmp_gtc_mk <- gtc_mk_reorder
    tmp_gtc_mk[which(tmp_gtc_mk == 99, arr.ind = TRUE)] <- "Q"
    gtc_mk_bind <- matrix("", nrow(gtc_mk), hnc)
    for(j in 1:hnc){
      gtc_mk_bind[, j] <- apply(tmp_gtc_mk[, c(2 * j - 1, 2 * j), drop = FALSE], 1, paste, collapse = ", ")
    }

    # Calculate posterior probabilities of each genotype combination
    gtc_logprob_mk <- gtc_logprob_use[[i]]
    loglike_gtc <- gtc_logprob_mk
    for(j in 1:n_rep){
      loglike_gtc <- loglike_gtc + lpd_gtc_use[[j]][[i]]
    }
    tmp <- exp(loglike_gtc)
    post_prob_gtc <- tmp / sum(tmp)

    # Make a data.table
    mat_pg <- cbind(gtc_mk_reorder, gtc_mk_bind, gtc_logprob_mk)
    for(j in 1:n_rep){
      mat_pg <- cbind(mat_pg, lpd_gtc_use[[j]][[i]])
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
#' @param hnc The number of contributors
#' @param cv_allele The coefficient of variation for allelels
#' @param cv_stutter The coefficient of variation for stutters
#' @param list_dt_average_sr A list of the data.table for the average stutter ratios
#' @param num_init The number of initial points for Pathfinder
perform_deconvo <- function(hnc, n_rep, markers_use, ats_use, peaks_list, sizes_list, heights_list, min_size,
                            pop_al_list, pop_prob_list, theta,
                            bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr, n_dropin_fltr, h_dropin_fltr,
                            freq_di, shape_di, rate_di,
                            ae_use, cv_allele, cv_stutter, list_dt_average_sr,
                            num_init){

  # --- Start computation --- #

  start_time <- proc.time()

  # --- Define objects --- #

  # Peak
  n_mk <- length(markers_use)
  n_peak_use <- rep(NA, n_mk)
  names(n_peak_use) <- markers_use
  tph_per_mk_use <- matrix(0, n_mk, n_rep)

  # Drop-in
  logprob_di_use <- rep(NA, n_mk)
  names(logprob_di_use) <- markers_use

  # Relative peak heights
  rph_base_use <- vector(mode = "list", length = n_mk)
  names(rph_base_use) <- markers_use

  # Genotype combinations
  n_gtc_use <- rep(NA, n_mk)
  gtc_use <- vector(mode = "list", length = n_mk)
  names(gtc_use) <- markers_use

  # Genotype probabilities
  gtc_logprob_use <- vector(mode = "list", length = n_mk)
  names(gtc_logprob_use) <- markers_use

  # De-duplicated objects
  n_general_list <- vector(mode = "list", length = n_rep)
  n_dropout_list <- vector(mode = "list", length = n_rep)
  n_dropin_list <- vector(mode = "list", length = n_rep)
  n_ignore_list <- vector(mode = "list", length = n_rep)
  index_dedup1_list_cpp <- vector(mode = "list", length = n_rep)
  cv_squared_dedup1_list <- vector(mode = "list", length = n_rep)
  deg_size_dedup1_list <- vector(mode = "list", length = n_rep)
  oph_dedup1_list <- vector(mode = "list", length = n_rep)
  loglike_di_list <- vector(mode = "list", length = n_rep)
  index_dedup2_list_cpp <- vector(mode = "list", length = n_rep)
  for(i in 1:n_rep){
    n_general_list[[i]] <- rep(NA, n_mk)
    n_dropout_list[[i]] <- rep(NA, n_mk)
    n_dropin_list[[i]] <- rep(NA, n_mk)
    n_ignore_list[[i]] <- rep(NA, n_mk)
    index_dedup1_list_cpp[[i]] <- vector(mode = "list", length = n_mk)
    cv_squared_dedup1_list[[i]] <- vector(mode = "list", length = n_mk)
    deg_size_dedup1_list[[i]] <- vector(mode = "list", length = n_mk)
    oph_dedup1_list[[i]] <- vector(mode = "list", length = n_mk)
    loglike_di_list[[i]] <- vector(mode = "list", length = n_mk)
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

    # --- Calculate total peak heights in each replicate --- #

    tph_per_mk_use[i, ] <- sapply(heights_all_rep, sum) / 2

    # --- Calculate the probability of drop-in in the current marker --- #

    auc_below_at <- pgamma(q = at, shape = shape_di, rate = rate_di)
    logprob_di_use[i] <- log(freq_di / (1 - auc_below_at))

    # --- Prepare the base of the relative peak heights --- #

    rph_base_use[[i]] <- rph_base_mk <- make_rph_base(peaks, bsrs, fsrs, dsrs, m2srs, ae)

    # --- Prepare the CV squared values for each peak --- #

    cv_square_base_mk <- c(rep(cv_allele^2, n_peak), rep(cv_stutter^2, 4 * n_peak), 0, -1)

    # --- List all possible genotype combinations --- #

    # All candidate genotypes of one contributor
    cand_gt <- combinations(n = n_peak, r = 2, v = peaks, repeats.allowed = TRUE)
    cand_gt <- asplit(cand_gt, 1)

    # Indices of genotypes of each contributor in all candidate genotype combinations
    index_gtc <- permutations(n = length(cand_gt), r = hnc, repeats.allowed = TRUE)
    index_gtc_cpp <- index_gtc - 1 # C++ index
    index_gtc_cpp <- asplit(index_gtc_cpp, 1)

    # Indices of genotypes of each contributor in each adopted genotype combination
    tmp <- judge_gtc_mk(cand_gt, index_gtc_cpp, peaks, heights_all_rep,
                        bsrs, fsrs, dsrs, m2srs,
                        bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr, n_dropin_fltr, h_dropin_fltr)
    index_gtc_adopt_cpp <- tmp[[1]]
    bool_gtc_di_mk_rep <- tmp[[2]]

    # The number of adopted genotype combinations
    n_gtc_use[i] <- n_gtc <- length(index_gtc_adopt_cpp)

    # --- Prepare objects for de-duplication (pattern 1) --- #

    if(!any(is.element(n_gtc_use, 0))){

      # Extract genotype combinations
      gtc_use[[i]] <- gtc_mk <- extract_gtc_adopt(cand_gt, index_gtc_adopt_cpp)

      # Calculate genotype probabilities
      n_obs <- make_n_obs(pop_al, prof_ref = numeric(0))
      obs_al <- as.numeric(names(n_obs))
      gtc_prob_mk <- calc_gtc_prob_mk(gtc_mk, pop_al, pop_prob, theta, prof_ref = numeric(0), hyp_id = rep(0, hnc), n_obs, obs_al)
      gtc_logprob_use[[i]] <- log(gtc_prob_mk)

      for(j in 1:n_rep){

        # Search the indices of the relative peak heights for each genotype combinations
        index_rph_gtc_cpp <- connect_rph_gtc_mk(cand_gt, index_gtc_adopt_cpp, peaks, bsrs, fsrs, dsrs, m2srs, bool_gtc_di_mk_rep[[j]])
        index_rph_gtc_cpp <- unlist(index_rph_gtc_cpp, recursive = FALSE)
        index_rph_gtc <- matrix(unlist(index_rph_gtc_cpp) + 1, nrow = 2 * hnc) # R and stan index

        # Delete duplicated indices
        bool_duplicated_adopt <- !duplicated(t(index_rph_gtc))
        index_dedup1_tmp <- index_rph_gtc[, bool_duplicated_adopt, drop = FALSE]

        # CV squared for de-duplicated objects (pattern 1)
        cv_square_dedup1_tmp <- matrix(cv_square_base_mk[index_dedup1_tmp], nrow = 2 * hnc)

        # Scale of degradation  corrected by sizes for de-duplicated objects (pattern 1)
        deg_size_dedup1_tmp <- rep((sizes - min_size) / min_size, n_gtc)
        deg_size_dedup1_tmp <- deg_size_dedup1_tmp[bool_duplicated_adopt]

        # Observe peak heights for de-duplicated objects (pattern 1)
        heights <- heights_all_rep[[j]]
        oph_dedup1_tmp <- rep(heights, n_gtc)
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

        # Assign cv squared for de-duplicated objects (pattern 1)
        cv_squared_dedup1_list[[j]][[i]] <- c(asplit(cv_square_dedup1_tmp[, index_general, drop = FALSE], 2),
                                              asplit(cv_square_dedup1_tmp[, index_dropout, drop = FALSE], 2))

        # Assign scale of degradation corrected by sizes for de-duplicated objects (pattern 1)
        deg_size_dedup1_list[[j]][[i]] <- c(deg_size_dedup1_tmp[index_general], deg_size_dedup1_tmp[index_dropout])

        # Assign observe peak heights for de-duplicated objects (pattern 1)
        oph_dedup1_list[[j]][[i]] <- c(oph_dedup1_tmp[index_general], oph_dedup1_tmp[index_dropout])

        # Assign log-likelihood of the drop-in for de-duplicated objects (pattern 1)
        if(length(index_dropin) != 0){
          loglike_di_list[[j]][[i]] <- log(dgamma(oph_dedup1_tmp[index_dropin], shape_di, rate_di)) + logprob_di_use[i]
        }else{
          loglike_di_list[[j]][[i]] <- numeric(0)
        }

        # --- Prepare objects for de-duplication (pattern 2) --- #
        index_dedup2_list_cpp[[j]][[i]] <- create_index_dedup2_mk(index_rph_gtc_cpp, index_dedup1_list_cpp[[j]][[i]])
      }
    }
  }

  if(!any(is.element(n_gtc_use, 0))){

    # --- Bayesian optimization --- #

    if(hnc == 1){

      # likelihood function
      like_func_single <- function(){

      }


    }else{

      # likelihood function
      like_func_mix <- function(...){
        # Prepare parameters
        args <- list(...)
        tph <- unlist(args[grepl("^tph", names(args))])[1:n_rep]
        lambda1 <- args$lambda1
        delta <- unlist(args[grepl("^delta", names(args))])[1:(hnc - 1)]
        deg <- unlist(args[grepl("^deg", names(args))])[1:hnc]

        # Prepare lambda
        lambda <- rep(NA, hnc)
        lambda[1] <- lambda1
        for(i in 1:(hnc - 1)){
          lambda[i + 1] <- lambda[i] + delta[i]
        }

        # Prepare mp
        mp <- softmax(lambda)

        # Calculate log probability densities of all genotype combinations
        lpd_gtc_all_rep <- vector(mode = "list", length = n_rep)
        for(i in 1:n_rep){
          lpd_gtc_all_rep[[i]] <- calc_lpd_gtc_rep1(tph, mp, deg, hnc, ats_use, n_general_list[[i]], n_dropout_list[[i]], n_dropin_list[[i]], n_ignore_list[[i]], rph_base_use,
                                                    index_dedup1_list_cpp[[i]], cv_squared_dedup1_list[[i]], deg_size_dedup1_list[[i]], oph_dedup1_list[[i]], loglike_di_list[[i]],
                                                    index_dedup2_list_cpp[[i]], n_gtc_use, n_peak_use)
        }

        # Calculate log likelihood in each marker
        loglike_cpp <- calc_loglike_cpp(lpd_gtc_all_rep, gtc_logprob_use)

        return(list(Score = - sum(loglike_cpp)))
      }

      bounds_list <- make_bounds_list(n_rep, hnc, max(tph_per_mk_use))

      opt_results <- bayesOpt(
        FUN = like_func_mix,
        bounds = bounds_list,
        initPoints = 4,     # 最初にランダムに試す点数
        iters.n = 10,       # ベイズ最適化で試行する回数
        acq = "ucb",        # 獲得関数の種類（ucb, ei, poi など）
        kappa = 2.576,      # UCBのパラメータ
        parallel = TRUE,    # 並列処理するかどうか
        verbose = 1
      )
    }



    dt_pg_list <- make_dt_pg_list(markers_use, paste0("U", 1:hnc), gtc_use, gtc_logprob_use, lpd_gtc_use)

    log10like_use <- calc_log10like_use(dt_pg_list)

    dt_log10like <- data.table(Marker = c(markers_use, "Overall"),
                               Log10_Likelihood = c(log10like_use, sum(log10like_use)))

    dt_param <- data.table("Parameter" = c(rep("Template amount", n_rep), rep("Mixture proportion", hnc), rep("Degradation", hnc)),
                           "Contributor" = c(rep(NA, n_rep), paste0("U", 1:hnc), paste0("U", 1:hnc)),
                           "Replicate" = c(1:n_rep, rep(NA, hnc), rep(NA, hnc)),
                           "Estimated value" = c(tph, mp[hnc:1], deg[hnc:1]))

    # End computation
    end_time <- proc.time()
    comp_time <- end_time[3] - start_time[3]

    return(list("dt_pg_list" = dt_pg_list,
                "dt_log10like" = dt_log10like,
                "dt_param" = dt_param,
                "logprob_di_use" = logprob_di_use,
                "comp_time" = comp_time,
                "init_list" = init_list,
                "warning_message" = warning_message,
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

    return(list("dt_pg_list" = NULL,
                "dt_log10like" = NULL,
                "dt_param" = NULL,
                "logprob_di_use" = NULL,
                "comp_time" = NULL,
                "init_list" = NULL,
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
#' @param num_init The number of initial points for Pathfinder
perform_deconvo_all_noc <- function(dt_csp, dt_af, dt_kit_auto, dt_ignore_mk, dt_anameth, list_dt_average_sr, hnc_all, theta, num_init){
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
                                      num_init))
      if(!inherits(data_hnc, "try-error")){
        data_hnc_list[[i]] <- data_hnc
      }else{
        data_hnc_list[[i]] <- list("dt_pg_list" = NULL,
                                   "dt_log10like" = NULL,
                                   "dt_param" = NULL,
                                   "logprob_di_use" = NULL,
                                   "comp_time" = NULL,
                                   "init_list" = NULL,
                                   "warning_message" = NULL,
                                   "error_message" = "Unexpected error")
      }
    }
  }else{
    peaks_list <- NULL
    heights_list <- NULL
    pop_al_list <- NULL
    pop_prob_list <- NULL
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
