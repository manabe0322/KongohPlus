#' view_dropin_ui
#'
#' @description The function to create the ui module for viewing drop-in parameters
view_dropin_ui <- function(id){
  ns <- NS(id)

  tabPanel(
    "View drop-in parameters",
    titlePanel("View drop-in parameters"),
    br()
  )
}

#' view_dropin_server
#'
#' @description The function to create the server module for viewing drop-in parameters
view_dropin_server <- function(id, path_pack){
  moduleServer(
    id,
    function(input, output, session){

    }
  )
}
