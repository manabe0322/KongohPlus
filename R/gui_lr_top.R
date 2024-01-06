#' lr_ui
#'
#' @description The function to create the ui module for the likelihood ratio
lr_ui <- function(id){
  ns <- NS(id)

  tabPanel("Likelihood ratio",
           titlePanel("Likelihood ratio"),
           br(),
           tabsetPanel(id = ns("tabset_lr"),
                       type = "tabs",
                       lr_setting_ui(ns("lr_setting")),
                       lr_result_ui(ns("lr_result"))
           )
  )
}

#' lr_server
#'
#' @description The function to create the server module for the likelihood ratio
lr_server <- function(id, path_pack, rv_file){
  moduleServer(
    id,
    function(input, output, session){
      lr_setting_server("lr_setting", path_pack, rv_file)
    }
  )
}
