#' KongohPlus: Probabilistic genotyping software based on a quantitative continuous model
#'
#' @importFrom combinat permn
#' @importFrom data.table as.data.table data.table copy fread setDT setindex
#' @importFrom dplyr full_join group_by_at summarise
#' @importFrom DT datatable dataTableOutput formatRound formatSignif formatStyle JS renderDataTable
#' @importFrom GenSA GenSA
#' @importFrom ggplot2 aes element_text geom_boxplot geom_col geom_line geom_point ggplot labs margin scale_color_manual scale_x_continuous scale_x_log10 scale_y_continuous scale_y_log10 sec_axis stat_function theme theme_minimal xlim ylim
#' @importFrom gridExtra grid.arrange
#' @importFrom gtools combinations ddirichlet permutations rdirichlet
#' @importFrom magrittr %>%
#' @importFrom parallel detectCores
#' @importFrom Rcpp sourceCpp
#' @importFrom rhandsontable hot_cell hot_col hot_table hot_to_r renderRHandsontable rhandsontable rHandsontableOutput
#' @importFrom shinyFeedback showFeedbackDanger hideFeedback useShinyFeedback
#' @importFrom shinyFiles parseFilePaths shinyFilesButton shinyFileChoose
#' @importFrom shinyjs disable disabled enable refresh useShinyjs
#' @importFrom shinythemes shinytheme
#' @importFrom stringr str_replace
#' @importFrom waiter useWaiter waiter_show waiter_hide spin_3k
"_PACKAGE"

#' @useDynLib KongohPlus, .registration = TRUE
NULL
