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

#' generate_init_sobol_dropin
#'
#' @description The function to generate initial values based on Sobol sequence for the drop-in model
#' @param x Sobol point
#' @param min_shape_di The minimum shape
#' @param max_shape_di The maximum shape
#' @param min_rate_di The minimum rate
#' @param max_rate_di The maximum rate
generate_init_sobol_dropin <- function(x, min_shape_di, max_shape_di, min_rate_di, max_rate_di){
  # Get only one Sobol point (x-th)
  sobol_vec <- sobol(n = x, d = 2, randomize = "digital.shift")[x, ]

  # Index partitioning and scaling
  shape_di <- sobol_vec[1] * (max_shape_di - min_shape_di) + min_shape_di
  rate_di <- sobol_vec[2] * (max_rate_di - min_rate_di) + min_rate_di

  return(list(list(shape_di = shape_di, rate_di = rate_di)))
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
#' @param max_iter The number of iterations for GenSA
#' @param num_converged_solutions The number of converged solutions required to terminate the MLE
#' @param max_mle_runs The maximum number of MLE attempts with different initial values
run_mle_dropin <- function(n_dropin, heights, min_shape_di, max_shape_di, min_rate_di, max_rate_di, max_iter, num_converged_solutions, max_mle_runs){

  like_func_dropin <- function(par){
    dens <- dgamma(heights, shape = par[1], rate = par[2])
    loglike <- - sum(log(dens))

    if(is.na(loglike) || is.nan(loglike) || is.infinite(loglike)){
      loglike <- Inf
    }

    return(loglike)
  }

  lower_bounds <- c(min_shape_di, min_rate_di)
  upper_bounds <- c(max_shape_di, max_rate_di)

  # Perform simulated annealing
  current_opt <- NULL
  current_score <- Inf
  current_score_round <- Inf
  count_mle <- 0
  count_solution <- 0
  init_list <- list()
  record_accept <- rep(NA, max_mle_runs)
  repeat{
    count_mle <- count_mle + 1

    init_list[[count_mle]] <- init <- c(runif(1, min = min_shape_di, max = max_shape_di), runif(1, min = min_rate_di, max = max_rate_di))

    opt <- GenSA(
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
    record_accept[count_mle] <- count_solution

    if(count_solution == num_converged_solutions){
      break
    }else if(count_mle == max_mle_runs){
      current_opt <- NULL
      break
    }
  }

  # Extract parameters
  if(!is.null(current_opt)){
    est_param <- current_opt$par
  }else{
    est_param <- rep(NA, 2)
  }
  names(est_param) <- c("shape_di", "rate_di")

  return(est_param)
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
#' @param max_iter The number of iterations for GenSA
#' @param num_converged_solutions The number of converged solutions required to terminate the MLE
#' @param max_mle_runs The maximum number of MLE attempts with different initial values
estimate_dropin_params <- function(dt_ntc, dt_kit, at_dropin, min_shape_di, max_shape_di, min_rate_di, max_rate_di, max_iter, num_converged_solutions, max_mle_runs){
  # Extract autosomal markers
  dt_kit_auto <- extract_kit_auto(dt_kit)
  markers_auto <- unique(dt_kit_auto[, Marker])
  n_mk <- length(markers_auto)

  # Extract drop-in data
  dt_ntc_auto <- dt_ntc[which(is.element(dt_ntc$Marker, markers_auto)), ]
  col_ntc <- names(dt_ntc_auto)
  heights <- as.numeric(unlist(dt_ntc_auto[, grep("Height", col_ntc), with = FALSE]))
  heights <- heights[!is.na(heights)]
  heights <- heights[heights >= at_dropin]
  n_dropin <- length(heights)
  n_sample <- length(unique(dt_ntc$SampleName))

  # Perform MLE
  tmp <- run_mle_dropin(n_dropin, heights, min_shape_di, max_shape_di, min_rate_di, max_rate_di, max_iter, num_converged_solutions, max_mle_runs)
  shape_di <- tmp["shape_di"]
  rate_di <- tmp["rate_di"]

  return(list("heights" = heights,
              "n_sample" = n_sample,
              "n_mk" = n_mk,
              "n_dropin" = n_dropin,
              "freq_dropin" = n_dropin / (n_sample * n_mk),
              "shape_di" = shape_di,
              "rate_di" = rate_di))
}
