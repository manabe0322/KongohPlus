#' KongohPlus: Probabilistic genotyping software based on a quantitative continuous model
#'
#' @importFrom combinat permn
#' @importFrom data.table as.data.table data.table copy fread setDT setindex
#' @importFrom dplyr arrange desc full_join group_by_at summarise
#' @importFrom DT datatable dataTableOutput formatRound formatSignif formatStyle JS renderDataTable
#' @importFrom GenSA GenSA
#' @importFrom ggplot2 aes element_text geom_col ggplot labs margin scale_x_continuous theme theme_minimal
#' @importFrom gridExtra grid.arrange
#' @importFrom gt cell_text cells_body cols_label gt gt_output opt_all_caps px render_gt sub_missing tab_options tab_spanner tab_style
#' @importFrom gtools combinations ddirichlet permutations rdirichlet
#' @importFrom magrittr %>%
#' @importFrom parallel detectCores
#' @importFrom Rcpp sourceCpp
#' @importFrom rhandsontable hot_cell hot_col hot_table hot_to_r renderRHandsontable rhandsontable rHandsontableOutput
#' @importFrom rmarkdown render
#' @importFrom shinyFeedback showFeedbackDanger hideFeedback useShinyFeedback
#' @importFrom shinyFiles parseFilePaths shinyFilesButton shinyFileChoose
#' @importFrom shinyjs disable disabled enable refresh useShinyjs
#' @importFrom shinythemes shinytheme
#' @importFrom stringr str_replace
#' @importFrom waiter useWaiter waiter_show waiter_hide spin_3k
"_PACKAGE"

#' @useDynLib KongohPlus, .registration = TRUE
NULL
