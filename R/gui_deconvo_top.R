#' deconvo_ui
#'
#' @description The function to create the ui module for deconvolution
#' @param id The ui ID
deconvo_ui <- function(id){
  ns <- NS(id)

  tabPanel("Deconvolution",
           tabsetPanel(id = ns("tabset_deconvo"),
                       type = "tabs",
                       deconvo_setting_ui(ns("deconvo_setting")),
                       deconvo_result_ui(ns("deconvo_result")),
                       deconvo_pg_ui(ns("deconvo_pg"))
           )
  )
}

#' deconvo_server
#'
#' @description The function to create the server module for deconvolution
#' @param id The server ID
#' @param path_pack The path of the package KongohPlus
#' @param session_top The top session
#' @param rv_input The reactive values derived from the input file server
deconvo_server <- function(id, path_pack, session_top, rv_input){
  moduleServer(
    id,
    function(input, output, session){
      rv_deconvo <- deconvo_setting_server("deconvo_setting", path_pack, rv_input)

      observe({
        req(rv_deconvo)
        deconvo_result_server("deconvo_result", session_top, rv_deconvo)
        deconvo_pg_server("deconvo_pg", rv_deconvo)
        if(!is.null(rv_deconvo$data_deconvo)){
          updateTabsetPanel(session, "tabset_deconvo", selected = "Result")
        }
      })

      return(rv_deconvo)
    }
  )
}
