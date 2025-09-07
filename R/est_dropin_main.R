#' check_est_dropin
#'
#' @description The function to check files for estimating dropin parameters
#' @param dt_ntc The data.table of the non template control data
#' @param dt_kit The data.table for information on the kit
check_est_dropin <- function(dt_ntc, dt_kit){
  error_message <- ""

  # Whether all required data is loaded or not
  if(is.null(dt_ntc)){
    error_message <- paste0(error_message, "Load non template control data!<br>")
  }

  if(!is.null(dt_ntc)){
    # Check the marker set of the experimental data
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
#' @param min_rate_di The minimum rate
#' @param max_rate_di The maximum rate
#' @param stop_iter Number of iterations without improvement before stopping one MLE run
#' @param max_iter Maximum number of iterations in one MLE run
#' @param num_converged Number of converged MLE runs before stopping
#' @param max_mle_runs Maximum number of MLE attempts
run_mle_dropin <- function(n_dropin, heights, min_shape_di, max_shape_di, min_rate_di, max_rate_di, stop_iter, max_iter, num_converged, max_mle_runs){

  like_func_dropin <- function(par){
    shape_di <- par[1] * (max_shape_di - min_shape_di) + min_shape_di
    rate_di <- par[2] * (max_rate_di - min_rate_di) + min_rate_di

    dens <- dgamma(heights, shape = shape_di, rate = rate_di)
    loglike <- - sum(log(dens))

    if(is.na(loglike) || is.nan(loglike) || is.infinite(loglike)){
      loglike <- Inf
    }

    return(loglike)
  }

  lower_bounds <- c(0, 0)
  upper_bounds <- c(1, 1)

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

    init_list[[count_mle]] <- init <- c(runif(1, min = min_shape_di, max = max_shape_di), runif(1, min = min_rate_di, max = max_rate_di))

    opt_list[[count_mle]] <- opt <- GenSA(
      par = init,
      fn = like_func_dropin,
      lower = lower_bounds,
      upper = upper_bounds,
      control = list(verbose = TRUE, maxit = max_iter)
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
      current_opt <- NULL
      break
    }
  }

  # Extract parameters
  est_param <- rep(NA, 2)
  names(est_param) <- c("shape_di", "rate_di")
  if(!is.null(current_opt)){
    tmp <- current_opt$par
    est_param[1] <- tmp[1] * (max_shape_di - min_shape_di) + min_shape_di
    est_param[2] <- tmp[2] * (max_rate_di - min_rate_di) + min_rate_di
  }

  return(list(est_param = est_param, init_list = init_list, opt_list = opt_list))
}

#' estimate_dropin_params
#'
#' @description The function to estimate drop-in parameters
#' @param dt_ntc The data.table for non-template controls
#' @param dt_kit The data.table for information on the kit
#' @param at_dropin The analytical threshold for drop-in peaks
#' @param min_shape_di The minimum shape
#' @param max_shape_di The maximum shape
#' @param min_rate_di The minimum rate
#' @param max_rate_di The maximum rate
#' @param stop_iter Number of iterations without improvement before stopping one MLE run
#' @param max_iter Maximum number of iterations in one MLE run
#' @param num_converged Number of converged MLE runs before stopping
#' @param max_mle_runs Maximum number of MLE attempts
estimate_dropin_params <- function(dt_ntc, dt_kit, at_dropin, min_shape_di, max_shape_di, min_rate_di, max_rate_di, stop_iter, max_iter, num_converged, max_mle_runs){

    # --- Extract data for autosomal markers --- #

  dt_kit_auto <- extract_kit_auto(dt_kit)
  markers_auto <- unique(dt_kit_auto[, Marker])
  n_mk_di <- length(markers_auto)
  dt_ntc_auto <- dt_ntc[which(is.element(dt_ntc$Marker, markers_auto)), ]

  # --- Create the data.table for experimental drop-in data --- #

  mat_alleles <- as.matrix(dt_ntc_auto[, grep("Allele", names(dt_ntc_auto)), with = FALSE])
  mat_heights <- as.matrix(dt_ntc_auto[, grep("Height", names(dt_ntc_auto)), with = FALSE])
  index_tmp <- which(!is.na(mat_heights), arr.ind = TRUE)
  sn_tmp <- dt_ntc_auto$SampleName[index_tmp[, 1]]
  al_tmp <- mat_alleles[index_tmp]
  ph_tmp <- mat_heights[index_tmp]

  # Extract data above AT
  index_tmp <- which(ph_tmp >= at_dropin)
  sn_tmp <- sn_tmp[index_tmp]
  al_tmp <- al_tmp[index_tmp]
  ph_tmp <- ph_tmp[index_tmp]

  dt_dropin_data <- data.table(SampleName = sn_tmp,
                               Allele = al_tmp,
                               Height = ph_tmp)

  # --- Perform MLE --- #

  heights <- dt_dropin_data$Heights
  n_dropin <- length(heights)

  tmp <- run_mle_dropin(n_dropin, heights, min_shape_di, max_shape_di, min_rate_di, max_rate_di, stop_iter, max_iter, num_converged, max_mle_runs)
  est_param <- tmp$est_param
  shape_di <- est_param["shape_di"]
  rate_di <- est_param["rate_di"]
  init_list_di <- tmp$init_list
  opt_list_di <- tmp$opt_list

  # --- Create the data.table for estimated parameters --- #

  n_sample_di <- length(unique(dt_ntc$SampleName))
  freq_dropin <- n_dropin / (n_sample * n_mk)
  dt_est_dropin <- data.table(Parameter = c("Number of samples", "Number of markers", "Number of drop-in peaks", "Frequency", "Shape parameter", "Rate parameter"),
                          Value = c(n_sample_di, n_mk_di, n_dropin, freq_dropin, shape_di, rate_di))

  return(list("dt_dropin_data" = dt_dropin_data,
              "dt_est_dropin" = dt_dropin,
              "init_list_di" = init_list_di,
              "opt_list_di" = opt_list_di))
}
