perform_deconvo <- function(path_pack, dt_csp, dt_af, dt_kit, dt_at, dt_param, dt_alcor, dt_gtc_fltr, dt_gtc_prob, hnc,
                            warmup, sampling, chains, di_alpha = 1.22, di_beta = 0.0503, logprob_di = log(0.00335), di_min = 49){
  # Extract the candidate mixture proportions
  mp_cand <- extract_mp_cand(hnc)
  mp_cand_tmp <- mp_cand * 0
  for(i in 1:hnc){
    mp_cand_tmp[, i] <- mp_cand[, hnc - i + 1]
  }
  mp_cand <- asplit(mp_cand, 1)

  # Extract autosomal markers
  dt_kit_auto <- extract_kit_auto(dt_kit)
  markers_auto <- unique(dt_kit_auto[, Marker])
  n_mk <- length(markers_auto)
  n_rep <- length(unique(dt_csp[, SampleName]))

  # Extract AT
  at_all_mk <- extract_at_all_mk(markers_auto, dt_at)

  # Extract peak data from crime stain profiles
  tmp <- make_peak_data(dt_csp, dt_kit_auto, markers_auto, at_all_mk)
  peaks_list <- tmp$peaks_list
  sizes_list <- tmp$sizes_list
  heights_list <- tmp$heights_list

  # Information on consideration of each stutter
  sr_consider_all_mk <- make_sr_consider_all_mk(markers_auto, dt_param)

  # Extract parameters for AE and M2SR
  ae_alpha_all_mk <- extract_param_all_mk("AE", "Alpha", markers_auto, dt_param)
  m2sr_all_mk <- extract_param_all_mk("M2SR", "Myu", markers_auto, dt_param)

  # Extract filters to remove unrealistic genotype combinations
  gtc_fltr_type <- dt_gtc_fltr[, Type]
  gtc_fltr_value <- dt_gtc_fltr[, Value]
  bsr_fltr <- gtc_fltr_value[gtc_fltr_type == "BSR"]
  fsr_fltr <- gtc_fltr_value[gtc_fltr_type == "FSR"]
  dsr_fltr <- gtc_fltr_value[gtc_fltr_type == "DSR"]
  m2sr_fltr <- gtc_fltr_value[gtc_fltr_type == "M2SR"]
  n_dropin_fltr <- gtc_fltr_value[gtc_fltr_type == "N_DropIn"]
  h_dropin_fltr <- gtc_fltr_value[gtc_fltr_type == "H_DropIn"]
  #phr_fltr <- gtc_fltr_value[gtc_fltr_type == "PHR"]
  #mr_fltr <- gtc_fltr_value[gtc_fltr_type == "MR"]
  #st <- gtc_fltr_value[gtc_fltr_type == "ST"]

  # Extract methods to estimate genotype probabilities
  gtc_prob_param <- dt_gtc_prob[, Parameter]
  gtc_prob_value <- dt_gtc_prob[, Value]
  theta <- as.numeric(gtc_prob_value[gtc_prob_param == "Theta"])

  # Set allele probabilities
  tmp <- set_al_prob(markers_auto, dt_af, peaks_list)
  pop_al_list <- tmp$pop_al_list
  pop_prob_list <- tmp$pop_prob_list

  # Extract gamma parameter
  gamma_beta <- extract_beta(dt_param)

  setkey(dt_alcor, Marker)

  n_eph_base_adopt <- rep(0, n_mk)
  eph_base <- vector(mode = "list", length = n_mk)
  n_dedup1_adopt <- n_dedup1_mcmc <- rep(0, n_mk)
  index_dedup1_adopt <- index_dedup1_mcmc <- NULL
  size_dedup1_adopt <- size_dedup1_mcmc <- vector(mode = "list", length = n_mk)
  oph_dedup1_adopt <- oph_dedup1_mcmc <- vector(mode = "list", length = n_mk)
  n_dedup2_adopt <- n_dedup2_mcmc <- rep(0, n_mk)
  index_dedup2_adopt <- index_dedup2_mcmc <- vector(mode = "list", length = n_mk)
  n_peak_adopt <- rep(0, n_mk)
  n_gtc_adopt <- n_gtc_mcmc <- rep(0, n_mk)
  gtc_logprob_adopt <- gtc_logprob_mcmc <- vector(mode = "list", length = n_mk)
  gtc_adopt <- vector(mode = "list", length = n_mk)
  names(gtc_adopt) <- markers_auto

  sumup_index_dedup1 <- 0
  max_tph <- 0
  for(i in 1:n_mk){
    # Extract peak data of one marker
    mk <- markers_auto[i]
    peaks <- peaks_list[[mk]]
    sizes <- sizes_list[[mk]]
    heights_all_rep <- heights_list[[mk]]
    n_peak <- length(peaks)

    # Extract AT of one marker
    at <- at_all_mk[mk]

    # Extract parameters of one marker
    sr_consider <- sr_consider_all_mk[mk, ]
    dt_alcor_mk <- dt_alcor[mk]
    alcor_alleles <- dt_alcor_mk[, Allele]
    alcor_bsrs <- dt_alcor_mk[, Average_BSR]
    alcor_fsrs <- dt_alcor_mk[, Average_FSR]
    alcor_dsrs <- dt_alcor_mk[, Average_DSR]
    index_alcor_peaks <- match(peaks, alcor_alleles)
    bsrs <- alcor_bsrs[index_alcor_peaks]
    bsrs[is.na(bsrs)] <- 0
    fsrs <- alcor_fsrs[index_alcor_peaks]
    fsrs[is.na(fsrs)] <- 0
    dsrs <- alcor_dsrs[index_alcor_peaks]
    dsrs[is.na(dsrs)] <- 0
    m2sr <- m2sr_all_mk[mk]

    # Extract allele probabilities
    pop_al <- pop_al_list[[mk]]
    pop_prob <- pop_prob_list[[mk]]

    # The base of the expected peak heights
    eph_base_mk <- make_eph_base(peaks, bsrs, fsrs, dsrs, m2sr)
    n_eph_base <- length(eph_base_mk)

    # All candidate genotypes of one contributor
    cand_gt <- combinations(n = n_peak, r = 2, v = peaks, repeats.allowed = TRUE)
    cand_gt <- asplit(cand_gt, 1)

    # Indices of genotypes of each contributor in all candidate genotype combinations
    index_gtc <- permutations(n = length(cand_gt), r = hnc, repeats.allowed = TRUE)
    index_gtc_cpp <- index_gtc - 1 # C++ index
    index_gtc_cpp <- asplit(index_gtc_cpp, 1)

    # Indices of genotypes of each contributor in each adopted genotype combination
    index_gtc_adopt_cpp <- judge_gtc_mk(cand_gt, index_gtc_cpp, peaks, heights_all_rep, sr_consider, bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr, n_dropin_fltr, h_dropin_fltr)
    n_gtc <- length(index_gtc_adopt_cpp)

    if(n_gtc != 0){
      # Extract genotype combinations
      gtc_mk <- extract_gtc_adopt(cand_gt, index_gtc_adopt_cpp)

      # Calculate genotype probabilities
      n_obs <- make_n_obs(pop_al, prof_ref = numeric(0))
      gtc_prob_mk <- calc_gtc_prob_mk(gtc_mk, pop_al, pop_prob, theta, prof_ref = numeric(0), hyp_id = c(0, 0), n_obs)
      gtc_logprob_mk <- log(gtc_prob_mk)

      # Connect eph_base and genotype combinations
      index_eph_gtc_cpp <- connect_eph_gtc_mk(cand_gt, index_gtc_adopt_cpp, peaks, sr_consider)
      index_eph_gtc <- matrix(unlist(index_eph_gtc_cpp) + 1, ncol = 2 * hnc, byrow = TRUE) # R and stan index

      # Delete duplicated indices -> index for deduplication (pattern 1: Before summation of the eph_base)
      bool_duplicated_adopt <- !duplicated(index_eph_gtc)
      index_dedup1_mk_adopt <- index_eph_gtc[bool_duplicated_adopt, , drop = FALSE]

      # Sizes for de-duplicated objects (pattern 1)
      size_tmp <- rep(sizes, n_gtc)
      size_dedup1_mk_adopt <- size_tmp[bool_duplicated_adopt]

      # Index for deduplication (pattern 2: After summation of the eph_base)
      index_dedup2_mk_adopt <- make_index_dedup2_mk(index_eph_gtc, index_dedup1_mk_adopt)

      # Narrow down genotype combinations for MCMC
      index_dedup1_mk_adopt_cpp <- asplit(index_dedup1_mk_adopt - 1, 2) # C++ index
      index_dedup2_mk_adopt_cpp <- index_dedup2_mk_adopt - 1 # C++ index
      index_gtc_mk_mcmc <- integer(0)
      oph_dedup1_mk_adopt <- numeric(0)
      for(j in 1:n_rep){
        heights <- heights_all_rep[[j]]
        oph_sum <- sum(heights)
        if(oph_sum > max_tph){
          max_tph <- oph_sum
        }
        oph_tmp <- rep(heights, n_gtc)
        oph_mk <- oph_tmp[bool_duplicated_adopt]
        oph_dedup1_mk_adopt <- c(oph_dedup1_mk_adopt, oph_mk)
        loglike_gtc_all_mp <- calc_loglike_gtc_all_mp(oph_sum, mp_cand, eph_base_mk, index_dedup1_mk_adopt_cpp, oph_mk, at, gamma_beta, di_alpha, di_beta, logprob_di, di_min, n_peak, n_gtc, index_dedup2_mk_adopt_cpp, gtc_logprob_mk)
        loglike_gtc_all_mp <- matrix(unlist(loglike_gtc_all_mp), ncol = n_gtc, byrow = TRUE)
        index_tmp <- which(loglike_gtc_all_mp >= max(loglike_gtc_all_mp) - 10, arr.ind = TRUE)
        index_gtc_mk_mcmc <- c(index_gtc_mk_mcmc, unique(index_tmp[, 2]))
      }
      index_gtc_mk_mcmc <- sort(unique(index_gtc_mk_mcmc))

      # Extract indices for deduplication used in MCMC
      n_gtc_nd <- length(index_gtc_mk_mcmc)
      index_extract_mcmc <- rep(0, n_peak * n_gtc_nd)
      for(k in 1:n_gtc_nd){
        index_extract_mcmc[n_peak * (k - 1) + 1:n_peak] <- n_peak * (index_gtc_mk_mcmc[k] - 1) + 1:n_peak
      }
      index_eph_gtc_mcmc <- index_eph_gtc[index_extract_mcmc, , drop = FALSE]
      bool_duplicated_mcmc <- !duplicated(index_eph_gtc_mcmc)
      index_dedup1_mk_mcmc <- index_eph_gtc_mcmc[bool_duplicated_mcmc, , drop = FALSE]
      index_dedup2_mk_mcmc <- make_index_dedup2_mk(index_eph_gtc_mcmc, index_dedup1_mk_mcmc)

      # Sizes for de-duplicated objects in MCMC
      size_tmp <- rep(sizes, n_gtc_nd)
      size_dedup1_mk_mcmc <- size_tmp[bool_duplicated_mcmc]

      # Observed peak heights for de-duplicated objects in MCMC
      oph_dedup1_mk_mcmc <- numeric(0)
      for(j in 1:n_rep){
        heights <- heights_all_rep[[j]]
        oph_tmp <- rep(heights, n_gtc_nd)
        oph_dedup1_mk_mcmc <- c(oph_dedup1_mk_mcmc, oph_tmp[bool_duplicated_mcmc])
      }

      # Sequential numbering of index_dedup1_mk
      index_dedup1_mk_adopt <- index_dedup1_mk_adopt + sumup_index_dedup1 # index for combining all loci
      index_dedup1_adopt <- rbind(index_dedup1_adopt, index_dedup1_mk_adopt)
      index_dedup1_mk_mcmc <- index_dedup1_mk_mcmc + sumup_index_dedup1 # index for combining all loci
      index_dedup1_mcmc <- rbind(index_dedup1_mcmc, index_dedup1_mk_mcmc)
      sumup_index_dedup1 <- sumup_index_dedup1 + n_eph_base

      # Assign objects
      n_eph_base_adopt[i] <- n_eph_base
      eph_base[[i]] <- eph_base_mk

      n_dedup1_adopt[i] <- nrow(index_dedup1_mk_adopt)
      n_dedup1_mcmc[i] <- nrow(index_dedup1_mk_mcmc)

      size_dedup1_adopt[[i]] <- size_dedup1_mk_adopt
      size_dedup1_mcmc[[i]] <- size_dedup1_mk_mcmc

      oph_dedup1_adopt[[i]] <- oph_dedup1_mk_adopt
      oph_dedup1_mcmc[[i]] <- oph_dedup1_mk_mcmc

      n_dedup2_adopt[i] <- length(index_dedup2_mk_adopt)
      n_dedup2_mcmc[i] <- length(index_dedup2_mk_mcmc)

      index_dedup2_adopt[[i]] <- index_dedup2_mk_adopt
      index_dedup2_mcmc[[i]] <- index_dedup2_mk_mcmc

      n_peak_adopt[i] <- n_peak
      n_gtc_adopt[i] <- n_gtc
      n_gtc_mcmc[i] <- n_gtc_nd

      gtc_logprob_adopt[[i]] <- gtc_logprob_mk
      gtc_logprob_mcmc[[i]] <- gtc_logprob_mk[index_gtc_mk_mcmc]

      gtc_adopt[[i]] <- gtc_mk
    }
  }

  ## check n_gt_comb_all == 0

  min_size <- min(unlist(sizes_list))

  if(n_rep == 1){
    data_list_main <- list(Hnc = hnc, N_Mk = n_mk, Ltph = max_tph / 4, Htph = max_tph / 2,
                           N_Eph_Base_All = sum(n_eph_base_adopt), Eph_Base = unlist(eph_base), AE_Alpha_All_Mk = ae_alpha_all_mk,
                           N_Dedup1_All = sum(n_dedup1_mcmc), N_Dedup1_Mcmc = n_dedup1_mcmc, Index_Dedup1_Mcmc = t(index_dedup1_mcmc),
                           Size_Dedup1_Mcmc = unlist(size_dedup1_mcmc) - min_size, Oph_Dedup1_Mcmc = unlist(oph_dedup1_mcmc),
                           N_Dedup2_All = sum(n_dedup2_mcmc), N_Dedup2_Mcmc = n_dedup2_mcmc, Index_Dedup2_Mcmc = unlist(index_dedup2_mcmc),
                           N_Peak_Adopt = n_peak_adopt, N_Gtc_Mcmc = n_gtc_mcmc, Gtc_Logprob_Mcmc = unlist(gtc_logprob_mcmc),
                           At_All_Mk = at_all_mk, Beta = gamma_beta,
                           Di_Alpha = di_alpha, Di_Beta = di_beta, Logprob_Di = logprob_di, Di_Min = di_min)

#    init_fun <- function(){
#      mp_init_pre <- runif(hnc)
#      mp_init <- sort(mp_init_pre, decreasing = TRUE) / sum(mp_init_pre)
#      deg_init <- runif(hnc, -0.0025, 0)
#      return(list(tph = max_tph / 3,
#                  mp = mp_init,
#                  deg = deg_init,
#                  ae = ae_alpha_all_mk / sum(ae_alpha_all_mk)))
#    }

    # Perform MCMC
    model_main <- readRDS(paste0(path_pack, "/model/model_main.rds"))
    fit_main <- model_main$sample(
      data_list_main,
      iter_warmup = warmup,
      iter_sampling = sampling,
      chains = chains,
      parallel_chains = detectCores()
    )
  }else{

  }

  # Trace plot
  draws <- fit_main$draws()
  draws %>% mcmc_trace(pars = c("lp__")) -> traceplot

  # Extract parameters
  fit_main_summary <- fit_main$summary()
  variables <- fit_main_summary$variable
  medians <- fit_main_summary$median
  rhats <- fit_main_summary$rhat
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

###############################################################################################################

  # Calculate likelihood
  index_dedup1_adopt_cpp <- asplit(index_dedup1_adopt - 1, 2) # C++ index
  eph_dedup1_adopt <- calc_eph_dedup1_cpp(mp, deg, hnc, sum(n_dedup1_adopt), unlist(eph_base), index_dedup1_adopt_cpp, unlist(size_dedup1_adopt))
  index_dedup2_adopt_cpp <- vector(mode = "list", length = n_mk)
  for(i in 1:n_mk){
    index_dedup2_adopt_cpp[[i]] <- index_dedup2_adopt[[i]] - 1
  }

  lpd_gtc_adopt <- calc_lpd_gtc_adopt(n_mk, tph, n_dedup1_adopt, unlist(oph_dedup1_adopt), eph_dedup1_adopt, at_all_mk, gamma_beta, di_alpha, di_beta, logprob_di, di_min, n_peak_adopt, n_gtc_adopt, index_dedup2_adopt_cpp)

  dt_pg_list <- make_dt_pg_list(markers_auto, gtc_adopt, lpd_gtc_adopt, gtc_logprob_adopt, paste0("U", 1:hnc))

  log10like_adopt <- calc_log10like_all_mk(lpd_gtc_adopt, gtc_logprob_adopt)

  dt_log10like <- make_dt_log10like(markers_auto, log10like_adopt)

  dt_mpdeg <- data.table(Contributor = paste0("U", 1:hnc),
                         MP = mp[hnc:1],
                         Deg = deg[hnc:1])

  dt_ae <- data.table(Marker = markers_auto,
                      AE = ae)

  dt_rhat <- data.table(Parameter = c("lp__", "tph", paste0("lambda", hnc:1), paste0("mp", hnc:1), paste0("deg", hnc:1), paste0("pae", 1:n_mk)),
                        Rhat = c(rhat_lp, rhat_tph, rhat_lambda[hnc:1], rhat_mp[hnc:1], rhat_deg[hnc:1], rhat_pae))

  return(list("hnc" = hnc,
              "markers_auto" = markers_auto,
              "pop_al_list" = pop_al_list,
              "pop_prob_list" = pop_prob_list,
              "theta" = theta,
              "gtc_adopt" = gtc_adopt,
              "lpd_gtc_adopt" = lpd_gtc_adopt,
              "dt_pg_list" = dt_pg_list,
              "dt_log10like" = dt_log10like,
              "dt_mpdeg" = dt_mpdeg,
              "tph" = tph,
              "lambda" = lambda,
              "dt_ae" = dt_ae,
              "dt_rhat" = dt_rhat,
              "traceplot" = traceplot))
}
