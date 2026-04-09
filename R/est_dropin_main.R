#' check_est_dropin
#'
#' @description The function to check files for estimating dropin parameters
#' @param dt_ntc The data.table of the non template control data
#' @param dt_kit The data.table for information on the kit
check_est_dropin <- function(dt_ntc, dt_kit){
  error_message <- ""

  # --- Check whether all required data is loaded --- #

  if(is.null(dt_ntc)){
    error_message <- paste0(error_message, "Load non template control data!<br>")
  }

  if(!is.null(dt_ntc)){

    # --- Check the marker set of the experimental data --- #

    sn_mk_problem_dropin <- check_dt_mk(dt_ntc, unique(dt_kit$Marker)) # The function is defined in est_sr_main.R
    if(length(sn_mk_problem_dropin) > 0){
      error_message <- paste0(error_message, "Inappropriate marker set of the following sample(s) in the non template control data: ", paste0(sn_mk_problem_dropin, collapse = ", "), ".<br>")
    }
  }

  return(error_message)
}

#' run_mle_dropin
#'
#' @description The function to perform MLE for the drop-in model
#' @param n_dropin The number of observed drop-in peaks
#' @param heights A vector of drop-in peak heights
#' @param min_shape_di The minimum shape
#' @param max_shape_di The maximum shape
#' @param min_scale_di The minimum scale
#' @param max_scale_di The maximum scale
#' @param max_iter Maximum number of iterations in one MLE run
#' @param num_consistent Minimum number of consistent runs
#' @param max_mle_runs Maximum number of MLE attempts
#' @param max_diff_like Reproducibility threshold (log10-likelihood difference)
#' @param waiter_ui Information on UI for waiter
run_mle_dropin <- function(n_dropin, heights, min_shape_di, max_shape_di, min_scale_di, max_scale_di, max_iter, num_consistent, max_mle_runs, max_diff_like, waiter_ui = NULL){

  like_func_dropin <- function(z){

    # shape_di
    u_shape_di <- plogis(z[1])
    shape_di <- min_shape_di + u_shape_di * (max_shape_di - min_shape_di)

    # scale_di
    u_scale_di <- plogis(z[2])
    scale_di <- min_scale_di + u_scale_di * (max_scale_di - min_scale_di)

    dens <- dgamma(heights, shape = shape_di, scale = scale_di)
    log10like <- - sum(log(dens)) / log(10)

    if(is.na(log10like) || is.nan(log10like) || is.infinite(log10like)){
      log10like <- Inf
    }

    return(log10like)
  }

  # Perform simulated annealing
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
          h3("Performing MLE for drop-in parameters..."),
          h3(paste0("Running MLE #", count_mle, " of ", max_mle_runs)),
          h3(paste0("Number of consistent runs so far: ", count_solution))
        )
      )
    }

    # Initial values in (0,1) from a Beta(1,1) distribution
    u0 <- rbeta(2, 1, 1)

    # Avoid exact 0 or 1 to prevent divergence in the logit transform since qlogis(0) = -Inf and qlogis(1) = +Inf
    u0 <- pmin(pmax(u0, 1e-6), 1-1e-6)

    # Transform to the unconstrained real space for CMA-ES
    init_z <- qlogis(u0)
    init_list[[count_mle]] <- init_z

    # Perform optimization
    opt_list[[count_mle]] <- opt <- cma_es(
      par = init_z,
      fn = like_func_dropin,
      control = list(maxit = max_iter)
    )

    score <- opt$value

    # Extract unconstrained parameters
    est_z <- opt$par

    # Back-transform to the original scales for shape_di
    u_shape_di  <- plogis(est_z[1])  # 0-1 scale via logistic transform
    shape_di <- min_shape_di + u_shape_di * (max_shape_di - min_shape_di)  # back-transform to the original scales

    # Back-transform to the original scales for scale_di
    u_scale_di <- plogis(est_z[2])  # 0-1 scale via logistic transform
    scale_di <- min_scale_di + u_scale_di * (max_scale_di - min_scale_di)  # back-transform to the original scales

    est_par_list[[count_mle]] <- c(- score, shape_di, scale_di)

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
  est_param <- mle_best[2:3]
  names(est_param) <- c("Shape parameter", "Scale parameter")

  # --- Arrange dt_all_mle --- #

  names(dt_all_mle) <- paste0("MLE ", 1:count_mle)
  dt_all_mle[, Parameter := c("Log10 (likelihood)", "Shape parameter", "Scale parameter")]
  setcolorder(dt_all_mle, c("Parameter", setdiff(names(dt_all_mle), "Parameter")))

  return(list("dt_all_mle" = dt_all_mle,
              "est_param" = est_param,
              "init_list" = init_list,
              "opt_list" = opt_list,
              "index_mle_best" = index_mle_best,
              "index_converged" = index_converged,
              "warning_message" = warning_message))
}

#' estimate_dropin_params
#'
#' @description The function to estimate drop-in parameters
#' @param dt_ntc The data.table for non-template controls
#' @param dt_kit The data.table for information on the kit
#' @param at_dropin The analytical threshold for drop-in peaks
#' @param min_shape_di The minimum shape
#' @param max_shape_di The maximum shape
#' @param min_scale_di The minimum scale
#' @param max_scale_di The maximum scale
#' @param max_iter Maximum number of iterations in one MLE run
#' @param num_consistent Minimum number of consistent runs
#' @param max_mle_runs Maximum number of MLE attempts
#' @param max_diff_like Reproducibility threshold (log10-likelihood difference)
#' @param waiter_ui Information on UI for waiter
estimate_dropin_params <- function(dt_ntc, dt_kit, at_dropin, min_shape_di, max_shape_di, min_scale_di, max_scale_di, max_iter, num_consistent, max_mle_runs, max_diff_like, waiter_ui){

  # --- Extract the specified conditions --- #

  dt_kit_auto <- extract_kit_auto(dt_kit)
  markers_auto <- unique(dt_kit_auto[, Marker])
  n_mk_di <- length(markers_auto)
  dt_ntc_auto <- dt_ntc[which(is.element(dt_ntc$Marker, markers_auto)), ]

  # --- Create the data.table for the experimental drop-in data --- #

  mat_alleles <- as.matrix(dt_ntc_auto[, grep("Allele", names(dt_ntc_auto)), with = FALSE])
  mat_heights <- as.matrix(dt_ntc_auto[, grep("Height", names(dt_ntc_auto)), with = FALSE])
  index_tmp <- which(!is.na(mat_heights), arr.ind = TRUE)
  sn_tmp <- dt_ntc_auto$SampleName[index_tmp[, 1]]
  al_tmp <- mat_alleles[index_tmp]
  ph_tmp <- mat_heights[index_tmp]

  # --- Extract data above AT --- #

  index_tmp <- which(ph_tmp >= at_dropin)
  sn_tmp <- sn_tmp[index_tmp]
  al_tmp <- al_tmp[index_tmp]
  ph_tmp <- ph_tmp[index_tmp]

  dt_dropin_data <- data.table(SampleName = sn_tmp,
                               Allele = al_tmp,
                               Height = ph_tmp)

  # --- Perform MLE --- #

  heights <- dt_dropin_data$Height
  n_dropin <- length(heights)

  tmp <- run_mle_dropin(n_dropin, heights, min_shape_di, max_shape_di, min_scale_di, max_scale_di, max_iter, num_consistent, max_mle_runs, max_diff_like, waiter_ui)
  dt_all_mle <- tmp$dt_all_mle
  est_param <- tmp$est_param
  shape_di <- est_param["Shape parameter"]
  scale_di <- est_param["Scale parameter"]
  init_list <- tmp$init_list
  opt_list <- tmp$opt_list
  index_mle_best <- tmp$index_mle_best
  index_converged <- tmp$index_converged
  warning_message <- tmp$warning_message

  # --- Create the data.table for estimated parameters --- #

  n_sample_di <- length(unique(dt_ntc$SampleName))
  freq_dropin <- n_dropin / (n_sample_di * n_mk_di)
  dt_est_dropin <- data.table(Parameter = c("Number of samples", "Number of markers", "Number of drop-in peaks", "Frequency", "Shape parameter", "Scale parameter"),
                          Value = c(n_sample_di, n_mk_di, n_dropin, freq_dropin, shape_di, scale_di))

  return(list("dt_all_mle" = dt_all_mle,
              "dt_dropin_data" = dt_dropin_data,
              "dt_est_dropin" = dt_est_dropin,
              "init_list" = init_list,
              "opt_list" = opt_list,
              "index_mle_best" = index_mle_best,
              "index_converged" = index_converged,
              "warning_message" = warning_message))
}
