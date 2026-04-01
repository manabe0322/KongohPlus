#' create_dt_ae_data
#'
#' @description The function to create a data.table of the locus-specific amplification efficiencies
#' @param peak_list A list of peaks
#' @param height_list A list of heights
create_dt_ae_data <- function(peak_list, height_list){
  sample_names <- names(peak_list)
  n_sample <- length(sample_names)
  markers <- names(peak_list[[1]])
  n_mk <- length(markers)

  n_data <- n_sample * n_mk
  sample_name_all <- rep("", n_data)
  marker_all <- rep("", n_data)
  temp_amount_all <- rep(NA, n_data)
  prop_ae_all <- rep(NA, n_data)
  ae_all <- rep(NA, n_data)

  current_index <- 1
  for(i in 1:n_sample){
    peak_list_one <- peak_list[[i]]
    height_list_one <- height_list[[i]]

    index_assign <- current_index:(current_index + n_mk - 1)
    sample_name_all[index_assign] <- sample_names[i]
    marker_all[index_assign] <- markers

    temp_amount <- sapply(height_list_one, sum)
    sum_temp_amount <- sum(temp_amount)
    temp_amount_all[index_assign] <- temp_amount
    prop_ae_all[index_assign] <- temp_amount / sum_temp_amount
    ae_all[index_assign] <- n_mk * temp_amount / sum_temp_amount

    current_index <- current_index + n_mk
  }

  dt_ae_data <- data.table(SampleName = sample_name_all,
                           Marker = marker_all,
                           TempAmount = temp_amount_all,
                           PropAE = prop_ae_all,
                           AE = ae_all)

  return(dt_ae_data)
}

#' run_mle_ae
#'
#' @description The function to perform MLE for the locus-specific amplification efficiency
#' @param prop_ae_data_use A matrix of the proportion of the amplification efficiencies (only data used for MLE)
#' @param min_alpha The minimum alpha
#' @param max_alpha The maximum alpha
#' @param max_iter Maximum number of iterations in one MLE run
#' @param num_consistent Minimum number of consistent runs
#' @param max_mle_runs Maximum number of MLE attempts
#' @param max_diff_like Reproducibility threshold (log10-likelihood difference)
#' @param waiter_ui Information on UI for waiter
run_mle_ae <- function(prop_ae_data_use, min_alpha, max_alpha, max_iter, num_consistent, max_mle_runs, max_diff_like, waiter_ui = NULL){

  # --- Define the likelihood function --- #

  like_fun_ae <- function(z){

    # alpha
    u_alpha <- plogis(z)
    alpha <- min_alpha + u_alpha * (max_alpha - min_alpha)

    n_data <- ncol(prop_ae_data_use)
    dens <- rep(Inf, n_data)

    for(i in 1:n_data){
      dens[i] <- ddirichlet(prop_ae_data_use[, i], alpha)
    }

    log10like <- - sum(log(dens)) / log(10)

    if(is.na(log10like) || is.nan(log10like) || is.infinite(log10like)){
      log10like <- Inf
    }

    return(log10like)
  }

  # --- Perform MLE --- #

  n_mk <- nrow(prop_ae_data_use)
  current_opt <- NULL
  current_score <- Inf
  count_solution <- 0
  init_list <- list()
  opt_list <- list()
  est_par_list <- list()
  index_converged <- integer(0)

  for(count_mle in 1:max_mle_runs){

    if(!is.null(waiter_ui)){
      waiter_ui$update(
        html = tagList(
          spin_3k(),
          h3("Performing MLE for LSAE..."),
          h3(paste0("Running MLE #", count_mle, " of ", max_mle_runs)),
          h3(paste0("Number of convergences so far: ", count_solution))
        )
      )
    }

    # Initial values in (0,1) from a Dirichlet distribution
    u0 <- rdirichlet(1, rep(1, n_mk))[1, ]

    # Avoid exact 0 or 1 to prevent divergence in the logit transform since qlogis(0) = -Inf and qlogis(1) = +Inf
    u0 <- pmin(pmax(u0, 1e-6), 1-1e-6)

    # Transform to the unconstrained real space for CMA-ES
    init_z <- qlogis(u0)
    init_list[[count_mle]] <- init_z

    # Perform optimization
    opt_list[[count_mle]] <- opt <- cma_es(
      par = init_z,
      fn = like_fun_ae,
      control = list(maxit = max_iter)
    )

    score <- opt$value

    # Extract unconstrained parameters
    est_z <- opt$par

    # Back-transform to the original scales for alpha
    u_alpha <- plogis(est_z)  # 0-1 scale via logistic transform
    alpha <- min_alpha + u_alpha * (max_alpha - min_alpha)  # back-transform to the original scales

    est_par_list[[count_mle]] <- c(- score, alpha)

    if(score <= current_score){

      if(current_score - score > max_diff_like){
        count_solution <- 1
        index_converged <- count_mle
      }else{
        count_solution <- count_solution + 1
        index_converged <- c(index_converged, count_mle)
      }

      current_opt <- opt
      current_score <- score
    }

    else if(score - current_score <= max_diff_like){
      count_solution <- count_solution + 1
      index_converged <- c(index_converged, count_mle)
    }

    if(length(index_converged) == num_consistent){
      break
    }
  }

  # --- Create the warning message --- #

  if(length(index_converged) == num_consistent){
    warning_message <- NULL
  }else{
    warning_message <- paste0("- The number of consistent runs is below ", num_consistent, ". This may indicate the presence of local optima.<br>")
  }

  # --- Extract the estimated parameters --- #

  dt_all_mle <- as.data.table(est_par_list)
  index_mle_best <- which.max(as.numeric(dt_all_mle[1, ]))
  mle_best <- dt_all_mle[[index_mle_best]]
  est_param <- mle_best[2:(n_mk + 1)]
  names(est_param) <- paste0("alpha[", 1:n_mk, "]")

  # --- Arrange dt_all_mle --- #

  names(dt_all_mle) <- paste0("MLE ", 1:count_mle)
  dt_all_mle[, Parameter := c("Log10 (likelihood)", paste0("alpha[", 1:n_mk, "]"))]
  setcolorder(dt_all_mle, c("Parameter", setdiff(names(dt_all_mle), "Parameter")))

  return(list("dt_all_mle" = dt_all_mle,
              "est_param" = est_param,
              "init_list" = init_list,
              "opt_list" = opt_list,
              "index_mle_best" = index_mle_best,
              "index_converged" = index_converged,
              "warning_message" = warning_message))
}

#' estimate_lsae
#'
#' @description The function to estimate parameters of the locus-specific amplification efficiencies
#' @param dt_exp The data.table for experimental data
#' @param dt_gt The data.table for genotype data
#' @param list_dt_average_sr A list of the data.table for average stutter ratios in all markers
#' @param dt_setting_est_par The data.table for the parameter settings
#' @param min_alpha The minimum alpha
#' @param max_alpha The maximum alpha
#' @param max_iter Maximum number of iterations in one MLE run
#' @param num_consistent Minimum number of consistent runs
#' @param max_mle_runs Maximum number of MLE attempts
#' @param max_diff_like Reproducibility threshold (log10-likelihood difference)
#' @param waiter_ui Information on UI for waiter
estimate_lsae <- function(dt_exp, dt_gt, list_dt_average_sr, dt_setting_est_par,
                          min_alpha, max_alpha,
                          max_iter, num_consistent, max_mle_runs, max_diff_like,
                          waiter_ui = NULL){

  # --- Extract the specified conditions --- #

  markers_auto <- dt_setting_est_par$Marker
  n_mk_auto <- length(markers_auto)
  at_auto <- dt_setting_est_par$AT
  names(at_auto) <- markers_auto

  # --- Information on considering the stutter ratios --- #

  sr_consider_auto <- matrix(TRUE, n_mk_auto, 4)
  rownames(sr_consider_auto) <- markers_auto
  colnames(sr_consider_auto) <- c("BS", "FS", "DS", "M2S")
  for(i in 1:n_mk_auto){
    mk <- markers_auto[i]
    dt_average_sr <- list_dt_average_sr[[mk]]
    sr_consider_auto[i, "BS"] <- any(dt_average_sr[["Average BSR"]] != 0)
    sr_consider_auto[i, "FS"] <- any(dt_average_sr[["Average FSR"]] != 0)
    sr_consider_auto[i, "DS"] <- any(dt_average_sr[["Average DSR"]] != 0)
    sr_consider_auto[i, "M2S"] <- any(dt_average_sr[["Average M2SR"]] != 0)
  }

  # --- Arrange experimental data --- #

  if(!is.null(waiter_ui)){
    waiter_ui$update(
      html = tagList(
        spin_3k(),
        h3("Preparing data...")
      )
    )
  }

  tmp <- arrange_exp_data(dt_exp, dt_gt, markers_auto, sr_consider_auto, at_auto)
  peak_list <- tmp$peak_list
  height_list <- tmp$height_list
  gt_list <- tmp$gt_list

  # --- Prepare data for LSAE --- #

  dt_ae_data <- create_dt_ae_data(peak_list, height_list)
  prop_ae_data <- matrix(dt_ae_data$PropAE, nrow = n_mk_auto)
  index_use <- which(apply(prop_ae_data, 2, prod) > 0)
  prop_ae_data_use <- prop_ae_data[, index_use, drop = FALSE]

  sn_all <- unique(dt_ae_data$SampleName)
  dt_prop_ae_data_use <- data.table(SampleName = sn_all[index_use])
  dt_prop_ae_data_use <- cbind(dt_prop_ae_data_use, as.data.table(t(prop_ae_data_use * n_mk_auto)))
  names(dt_prop_ae_data_use) <- c("SampleName", markers_auto)

  # --- Perform MLE for LSAE --- #

  tmp <- run_mle_ae(prop_ae_data_use, min_alpha, max_alpha, max_iter, num_consistent, max_mle_runs, max_diff_like, waiter_ui)
  dt_all_mle <- tmp$dt_all_mle
  alpha <- tmp$est_param
  init_list <- tmp$init_list
  opt_list <- tmp$opt_list
  index_mle_best <- tmp$index_mle_best
  index_converged <- tmp$index_converged
  warning_message <- tmp$warning_message

  # --- Create the data.table for LSAE --- #

  expect_ae <- n_mk_auto * (alpha + 1) / (sum(alpha) + n_mk_auto)
  dt_expect_ae <- data.table(Marker = markers_auto, Alpha = alpha, LSAE = expect_ae)

  return(list("dt_expect_ae" = dt_expect_ae,
              "dt_all_mle" = dt_all_mle,
              "init_list" = init_list,
              "opt_list" = opt_list,
              "index_mle_best" = index_mle_best,
              "index_converged" = index_converged,
              "warning_message" = warning_message,
              "dt_prop_ae_data_use" = dt_prop_ae_data_use))
}
