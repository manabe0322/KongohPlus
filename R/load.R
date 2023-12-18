#' load_ui
#'
#' @description The function to create the ui module for loading files
load_ui <- function(id){
  ns <- NS(id)

  tabPanel("Load",
           fluidRow(
             column(4,
                    h2("Crime stain profile"),
                    fileInput(ns("file_csp"), label = NULL, accept = ".csv"),
                    dataTableOutput(ns("info_csp"))
             ),
             column(4,
                    h2("Reference profile"),
                    fileInput(ns("file_ref"), label = NULL, accept = ".csv"),
                    dataTableOutput(ns("info_ref"))
             ),
             column(4,
                    h2("Allele probability"),
                    fileInput(ns("file_af"), label = NULL, accept = ".csv")
             )
           ),
           br(),
           actionButton(ns("act_deconvo"), label = "Deconvolution", class = "btn btn-primary btn-lg"),
           br(),
           br(),
           actionButton(ns("act_lr"), label = "Likelihood ratio", class = "btn btn-primary btn-lg")
  )
}

#' load_server
#'
#' @description The function to create the server module for loading files
load_server <- function(id){
  moduleServer(
    id,
    function(input, output, session){

    }
  )
}
