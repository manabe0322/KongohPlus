#' KongohPlus: Probabilistic genotyping software based on a quantitative continuous model
#'
#' @importFrom cmaes cma_es
#' @importFrom combinat permn
#' @importFrom data.table as.data.table data.table copy fread setcolorder setDT setindex
#' @importFrom dplyr arrange desc full_join group_by_at summarise
#' @importFrom DT datatable dataTableOutput formatRound formatSignif formatStyle JS renderDataTable
#' @importFrom ggplot2 aes coord_cartesian element_text geom_abline geom_boxplot geom_col geom_hline geom_line geom_point geom_ribbon ggplot labs margin scale_color_manual scale_fill_manual scale_x_continuous scale_x_log10 scale_y_continuous scale_y_log10 sec_axis stat_function theme theme_minimal xlab xlim ylab ylim
#' @importFrom gridExtra grid.arrange
#' @importFrom gtools combinations ddirichlet permutations rdirichlet
#' @importFrom htmltools tags withTags
#' @importFrom magrittr %>%
#' @importFrom parallel detectCores
#' @importFrom Rcpp sourceCpp
#' @importFrom rhandsontable hot_cell hot_col hot_table hot_to_r renderRHandsontable rhandsontable rHandsontableOutput
#' @importFrom rmarkdown find_pandoc html_document render
#' @importFrom shinyFeedback showFeedbackDanger hideFeedback useShinyFeedback
#' @importFrom shinyFiles parseFilePaths shinyFilesButton shinyFileChoose
#' @importFrom shinyjs disable disabled enable refresh useShinyjs
#' @importFrom shinythemes shinytheme
#' @importFrom stringr str_replace
#' @importFrom waiter useWaiter spin_3k transparent Waiter
"_PACKAGE"

#' @useDynLib KongohPlus, .registration = TRUE
NULL
