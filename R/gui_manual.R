#' manual_ui
#'
#' @description The function to create the ui module for the manual
manual_ui <- function(id){
  ns <- NS(id)

  tabPanel(
    "User manual",
    titlePanel("User manual"),
    br(),
    downloadButton(ns("download_manual_jp"), "Manual (Japanese)", class = "btn btn-primary btn-lg")
  )
}

#' manual_server
#'
#' @description The function to create the server module for the manual
manual_server <- function(id, ver_soft, path_pack){
  moduleServer(
    id,
    function(input, output, session){

      output$download_manual_jp <- downloadHandler(
        filename = paste0("KongohPlus_", ver_soft, "_Manual.pdf"),
        content = function(file){
          file.copy(paste0(path_pack, "/extdata/manual/", "KongohPlus_", ver_soft, "_Manual_JP.pdf"), file)
        }
      )
    }
  )
}
