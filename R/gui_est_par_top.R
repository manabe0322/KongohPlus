#' est_par_ui
#'
#' @description The function to create the ui module for estimating parameters
est_par_ui <- function(id){
  ns <- NS(id)

  tabPanel("4. Estimate LSAE & peak height parameters",
           titlePanel("4. Estimate LSAE & peak height parameters"),
           br(),
           tabsetPanel(id = ns("tabset_est_par"),
                       type = "tabs",
                       est_par_setting_ui(ns("est_par_setting")),
                       est_par_result_ui(ns("est_par_result"))
           )
  )
}

#' est_par_server
#'
#' @description The function to create the server module for estimating parameters
est_par_server <- function(id, path_pack){
  moduleServer(
    id,
    function(input, output, session){
      rv_est_par <- est_par_setting_server("est_par_setting", path_pack, session)

      observe({
        req(rv_est_par)
        est_par_result_server("est_par_result", path_pack, rv_est_par)
      })
    }
  )
}
