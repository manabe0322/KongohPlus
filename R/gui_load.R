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
load_server <- function(id, session_top){
  moduleServer(
    id,
    function(input, output, session){
      rv_file <- reactiveValues()
      rv_file$dt_csp <- NULL
      rv_file$dt_ref <- NULL
      rv_file$dt_af <- NULL

      # Crime stain profile
      observeEvent(input$file_csp, {
        file_input <- input$file_csp
        if(!is.null(file_input)){
          dt_csp <- fread(file_input$datapath)
          dt_csp <- change_data_type_csp(dt_csp)
          rv_file$dt_csp <- dt_csp
        }
      })

      # Reference profile
      observeEvent(input$file_ref, {
        file_input <- input$file_ref
        if(!is.null(file_input)){
          dt_ref <- fread(file_input$datapath)
          dt_ref <- change_data_type_ref(dt_ref)
          rv_file$dt_ref <- dt_ref
        }
      })

      # Allele frequency
      observeEvent(input$file_af, {
        file_input <- input$file_af
        if(!is.null(file_input)){
          dt_af <- fread(file_input$datapath)
          dt_af <- change_data_type_af(dt_af)
          rv_file$dt_af <- dt_af
        }
      })

      observeEvent(input$act_deconvo, {
        updateNavbarPage(session = session_top, inputId = "navbar", selected = "Deconvolution")
      })

      observeEvent(input$act_lr, {
        updateNavbarPage(session = session_top, inputId = "navbar", selected = "Likelihood ratio")
      })

      return(rv_file)
    }
  )
}
