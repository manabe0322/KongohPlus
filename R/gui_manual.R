#' manual_ui
#'
#' @description The function to create the ui module for the manual
manual_ui <- function(id){
  ns <- NS(id)

  tabPanel(
    "User manual",
    titlePanel("User manual"),
    br(),
    downloadButton(ns("download_manual_en"), "Manual (English)", class = "btn btn-primary btn-lg"),
    br(),
    br(),
    downloadButton(ns("download_manual_jp"), "Manual (Japanese)", class = "btn btn-primary btn-lg")
  )
}

#' manual_server
#'
#' @description The function to create the server module for the manual
#' @param id The server ID
#' @param ver_soft The package version
#' @param path_pack The path of the package KongohPlus
manual_server <- function(id, ver_soft, path_pack){
  moduleServer(
    id,
    function(input, output, session){

      output$download_manual_en <- downloadHandler(
        filename = paste0("UserManual_KongohPlus_EN_v", ver_soft, ".pdf"),
        content = function(file){
          file.copy(paste0(path_pack, "/extdata/manual/", "UserManual_KongohPlus_EN_v", ver_soft, ".pdf"), file)
        }
      )

      output$download_manual_jp <- downloadHandler(
        filename = paste0("UserManual_KongohPlus_JP_v", ver_soft, ".pdf"),
        content = function(file){
          file.copy(paste0(path_pack, "/extdata/manual/", "UserManual_KongohPlus_JP_v", ver_soft, ".pdf"), file)
        }
      )
    }
  )
}
