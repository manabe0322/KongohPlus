#' KongohPlus: Probabilistic genotyping software based on a quantitative continuous model
#'
#' @importFrom bayesplot mcmc_trace
#' @importFrom cmdstanr cmdstan_model
#' @importFrom data.table as.data.table data.table copy fread setDT setindex
#' @importFrom dplyr arrange desc full_join group_by_at summarise
#' @importFrom DT datatable dataTableOutput formatRound formatSignif formatStyle JS renderDataTable
#' @import gluedown
#' @importFrom gtools combinations permutations
#' @importFrom instantiate stan_package_model
#' @importFrom pandoc pandoc_activate pandoc_install
#' @importFrom magrittr %>%
#' @importFrom parallel detectCores
#' @importFrom Rcpp sourceCpp
#' @importFrom rhandsontable hot_cell hot_col hot_table hot_to_r renderRHandsontable rhandsontable rHandsontableOutput
#' @import rmarkdown
#' @importFrom shinyFeedback showFeedbackDanger hideFeedback useShinyFeedback
#' @importFrom shinyFiles parseFilePaths shinyFilesButton shinyFileChoose
#' @importFrom shinyjs disable disabled enable refresh useShinyjs
#' @importFrom shinythemes shinytheme
#' @importFrom waiter useWaiter waiter_show waiter_hide spin_3k
"_PACKAGE"

#' @useDynLib KongohPlus, .registration = TRUE
NULL
