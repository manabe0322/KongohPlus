#' est_lsae_ui
#'
#' @description The function to create the ui module for estimating locus-specific amplification efficiencies
#' @param id The ui ID
est_lsae_ui <- function(id){
  ns <- NS(id)

  tabPanel("4. Estimate the locus-specific amplification efficiencies",
           titlePanel("4. Estimate the locus-specific amplification efficiencies"),
           br(),
           tabsetPanel(id = ns("tabset_est_lsae"),
                       type = "tabs",
                       est_lsae_setting_ui(ns("est_lsae_setting")),
                       est_lsae_result_ui(ns("est_lsae_result_lsae"))
           )
  )
}

#' est_lsae_server
#'
#' @description The function to create the server module for estimating locus-specific amplification efficiencies
#' @param id The server ID
#' @param path_pack The path of the package KongohPlus
est_lsae_server <- function(id, path_pack){
  moduleServer(
    id,
    function(input, output, session){
      rv_est_lsae <- est_lsae_setting_server("est_lsae_setting", path_pack, session)

      observe({
        req(rv_est_lsae)
        est_lsae_result_server("est_lsae_result_lsae", path_pack, rv_est_lsae)
      })
    }
  )
}
