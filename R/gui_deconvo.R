#' deconvo_ui
#'
#' @description The function to create the ui module for deconvolution
deconvo_ui <- function(id){
  ns <- NS(id)

  tabPanel("Deconvolution",
           sidebarLayout(
             sidebarPanel(width = 3,
                          numericInput(ns("noc_deconvo"), label = h4("Number of contributors"), value = NULL)
             ),
             mainPanel(width = 9,
                       selectInput(ns("analysis_method_deconvo"), label = h4("Analysis method"), choices = LETTERS[1:3])
             )
           )
  )
}

#' deconvo_server
#'
#' @description The function to create the server module for deconvolution
deconvo_server <- function(id){
  moduleServer(
    id,
    function(input, output, session){

    }
  )
}
