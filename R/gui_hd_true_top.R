#' hd_true_ui
#'
#' @description The function to create the ui module for the Hd true test
hd_true_ui <- function(id){
  ns <- NS(id)

  tabPanel("Hd true test",
           titlePanel("Hd true test"),
           tabsetPanel(id = ns("tabset_hd_true"),
                       type = "tabs",
                       hd_true_setting_ui(ns("hd_true_setting")),
                       hd_true_result_ui(ns("hd_true_result"))
           )
  )
}

#' hd_true_server
#'
#' @description The function to create the server module for the Hd true test
hd_true_server <- function(id, path_pack, rv_file, rv_deconvo){
  moduleServer(
    id,
    function(input, output, session){
      rv_hd_true <- hd_true_setting_server("hd_true_setting", path_pack, rv_file, rv_deconvo)

      observe({
        req(rv_hd_true)
        hd_true_result_server("hd_true_result", rv_hd_true)
        if(!is.null(rv_hd_true$data_hd_true)){
          updateTabsetPanel(session, "tabset_hd_true", selected = "Result")
        }
      })

      return(rv_hd_true)
    }
  )
}
