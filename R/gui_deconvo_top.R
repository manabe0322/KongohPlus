#' deconvo_ui
#'
#' @description The function to create the ui module for deconvolution
deconvo_ui <- function(id){
  ns <- NS(id)

  tabPanel("Deconvolution",
           titlePanel("Deconvolution"),
           br(),
           tabsetPanel(id = ns("tabset_deconvo"),
                       type = "tabs",
                       deconvo_setting_ui(ns("deconvo_setting")),
                       deconvo_result_ui(ns("deconvo_result"))
           )
  )
}

#' deconvo_server
#'
#' @description The function to create the server module for deconvolution
deconvo_server <- function(id, path_pack){
  moduleServer(
    id,
    function(input, output, session){
      deconvo_setting_server("deconvo_setting", path_pack)
    }
  )
}
