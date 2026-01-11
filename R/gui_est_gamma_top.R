#' est_gamma_ui
#'
#' @description The function to create the ui module for estimating the gamma scale parameter
#' @param id The ui ID
est_gamma_ui <- function(id){
  ns <- NS(id)

  tabPanel("5. Estimate the gamma scale parameter",
           titlePanel("5. Estimate the gamma scale parameter"),
           br(),
           tabsetPanel(id = ns("tabset_est_gamma"),
                       type = "tabs",
                       est_gamma_setting_ui(ns("est_gamma_setting")),
                       est_gamma_result_ui(ns("est_gamma_result"))
           )
  )
}

#' est_gamma_server
#'
#' @description The function to create the server module for estimating the gamma scale parameter
#' @param id The server ID
#' @param path_pack The path of the package KongohPlus
est_gamma_server <- function(id, path_pack){
  moduleServer(
    id,
    function(input, output, session){
      rv_est_gamma <- est_gamma_setting_server("est_gamma_setting", path_pack, session)

      observe({
        req(rv_est_gamma)
        est_gamma_result_server("est_gamma_result", path_pack, rv_est_gamma)
      })
    }
  )
}
