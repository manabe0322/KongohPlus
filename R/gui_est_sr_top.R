#' est_sr_ui
#'
#' @description The function to create the ui module for estimating the average stutter ratios
#' @param id The ui ID
est_sr_ui <- function(id){
  ns <- NS(id)

  tabPanel("3. Estimate the stutter ratios",
           titlePanel("3. Estimate the stutter ratios"),
           br(),
           tabsetPanel(id = ns("tabset_est_sr"),
                       type = "tabs",
                       est_sr_setting_ui(ns("est_sr_setting")),
                       est_sr_result_ui(ns("est_sr_result"))
           )
  )
}

#' est_sr_server
#'
#' @description The function to create the server module for estimating the average stutter ratios
#' @param id The server ID
#' @param path_pack The path of the package KongohPlus
est_sr_server <- function(id, path_pack){
  moduleServer(
    id,
    function(input, output, session){
      rv_est_sr <- est_sr_setting_server("est_sr_setting", path_pack, session)

      observe({
        req(rv_est_sr)
        est_sr_result_server("est_sr_result", path_pack, rv_est_sr)
      })

      return(rv_est_sr)
    }
  )
}
