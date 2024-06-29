#' KongohPlus: Probabilistic genotyping software based on a quantitative continuous model
#'
#' @importFrom bayesplot mcmc_trace
#' @importFrom cmdstanr cmdstan_model
#' @importFrom data.table as.data.table data.table copy fread setDT setkey
#' @importFrom dplyr full_join
#' @importFrom DT datatable dataTableOutput formatRound formatSignif formatStyle JS renderDataTable
#' @import gluedown
#' @importFrom gtools combinations permutations
#' @importFrom instantiate stan_package_model
#' @importFrom pandoc pandoc_activate pandoc_install
#' @importFrom magrittr %>%
#' @importFrom parallel detectCores
#' @importFrom posterior as_draws_matrix rhat
#' @importFrom Rcpp sourceCpp
#' @import rmarkdown
#' @importFrom shinyFeedback showFeedbackDanger hideFeedback useShinyFeedback
#' @importFrom shinyjs disable disabled enable useShinyjs
#' @importFrom shinythemes shinytheme
"_PACKAGE"

#' @useDynLib KongohPlus, .registration = TRUE
NULL
