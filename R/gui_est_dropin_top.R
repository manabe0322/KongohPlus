#' est_dropin_ui
#'
#' @description The function to create the ui module for estimating drop-in parameters
est_dropin_ui <- function(id){
  ns <- NS(id)

  tabPanel("6. Estimate drop-in parameters",
           titlePanel("6. Estimate drop-in parameters"),
           br(),
           tabsetPanel(id = ns("tabset_est_dropin"),
                       type = "tabs",
                       est_dropin_setting_ui(ns("est_dropin_setting")),
                       est_dropin_result_ui(ns("est_dropin_result"))
           )
  )
}

#' est_dropin_server
#'
#' @description The function to create the server module for estimating drop-in parameters
est_dropin_server <- function(id, path_pack){
  moduleServer(
    id,
    function(input, output, session){
      rv_est_dropin <- est_dropin_setting_server("est_dropin_setting", path_pack, session)

      observe({
        req(rv_est_dropin)
        est_dropin_result_server("est_dropin_result", path_pack, rv_est_dropin)
      })

      return(rv_est_dropin)
    }
  )
}
