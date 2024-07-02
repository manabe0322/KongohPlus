#' kit_ui
#'
#' @description The function to create the ui module for the kit
kit_ui <- function(id){
  ns <- NS(id)

  tabPanel("Kit",
           titlePanel("Kit"),
           br()
  )
}

#' kit_server
#'
#' @description The function to create the server module for the kit
kit_server <- function(id, path_pack){
  moduleServer(
    id,
    function(input, output, session){

    }
  )
}
