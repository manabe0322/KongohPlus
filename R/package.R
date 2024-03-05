#' KongohPlus: Probabilistic genotyping software based on a quantitative continuous model
#'
#' @importFrom bayesplot mcmc_trace
#' @importFrom cmdstanr cmdstan_model
#' @importFrom data.table as.data.table data.table copy fread setDT setkey
#' @importFrom dplyr full_join
#' @importFrom DT datatable dataTableOutput formatRound formatSignif JS renderDataTable
#' @importFrom gtools combinations permutations
#' @importFrom magrittr %>%
#' @importFrom MCMCpack rdirichlet
#' @importFrom parallel detectCores
#' @importFrom posterior as_draws_matrix rhat
#' @importFrom Rcpp sourceCpp
#' @importFrom shinyFeedback showFeedbackDanger hideFeedback useShinyFeedback
#' @importFrom shinyjs disable disabled enable useShinyjs
#' @importFrom shinythemes shinytheme
"_PACKAGE"

#' @useDynLib KongohPlus, .registration = TRUE
NULL
