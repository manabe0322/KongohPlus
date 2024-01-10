#' load_ui
#'
#' @description The function to create the ui module for loading files
load_ui <- function(id){
  ns <- NS(id)

  tabPanel("Load",
           titlePanel("Load data"),
           br(),
           sidebarLayout(
             sidebarPanel(width = 3,
                          fileInput(ns("file_csp"), label = h4("Crime stain profile"), accept = ".csv", buttonLabel = "Load"),
                          fileInput(ns("file_ref"), label = h4("Reference profile"), accept = ".csv", buttonLabel = "Load"),
                          fileInput(ns("file_ladder"), label = h4("Allelic ladder"), accept = ".csv", buttonLabel = "Load"),
                          fileInput(ns("file_af"), label = h4("Allele frequency"), accept = ".csv", buttonLabel = "Load"),
                          br(),
                          actionButton(ns("act_deconvo"), label = "Deconvolution", class = "btn btn-primary btn-lg"),
                          br(),
                          br(),
                          br(),
                          actionButton(ns("act_lr"), label = "Likelihood ratio", class = "btn btn-primary btn-lg")
             ),
             mainPanel(width = 9,
                       tabsetPanel(id = ns("tabset_load_file"),
                                   type = "tabs",
                                   tabPanel("Crime stain profile",
                                            dataTableOutput(ns("dt_csp"))
                                   ),
                                   tabPanel("Reference profile",
                                            dataTableOutput(ns("dt_ref"))
                                   ),
                                   tabPanel("Allelic ladder",
                                            dataTableOutput(ns("dt_ladder"))
                                   ),
                                   tabPanel("Allele frequency",
                                            dataTableOutput(ns("dt_af"))
                                   )
                       )
             )
           )
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
      rv_file$dt_ladder <- NULL
      rv_file$dt_af <- NULL

      # Load the crime stain profile
      observeEvent(input$file_csp, {
        file_input <- input$file_csp
        if(!is.null(file_input)){
          dt_csp <- fread(file_input$datapath)
          dt_csp <- change_data_type_csp(dt_csp)
          rv_file$dt_csp <- dt_csp
          updateTabsetPanel(inputId = "tabset_load_file", selected = "Crime stain profile")
        }
      })

      # Load the reference profile
      observeEvent(input$file_ref, {
        file_input <- input$file_ref
        if(!is.null(file_input)){
          dt_ref <- fread(file_input$datapath)
          dt_ref <- change_data_type_ref(dt_ref)
          rv_file$dt_ref <- dt_ref
          updateTabsetPanel(inputId = "tabset_load_file", selected = "Reference profile")
        }
      })

      # Load the allelic ladder
      observeEvent(input$file_ladder, {
        file_input <- input$file_ladder
        if(!is.null(file_input)){
          dt_ladder <- fread(file_input$datapath)
          dt_ladder <- change_data_type_csp(dt_ladder)
          rv_file$dt_ladder <- dt_ladder
          updateTabsetPanel(inputId = "tabset_load_file", selected = "Allelic ladder")
        }
      })

      # Load the allele frequency
      observeEvent(input$file_af, {
        file_input <- input$file_af
        if(!is.null(file_input)){
          dt_af <- fread(file_input$datapath)
          dt_af <- change_data_type_af(dt_af)
          rv_file$dt_af <- dt_af
          updateTabsetPanel(inputId = "tabset_load_file", selected = "Allele frequency")
        }
      })

      # Display the crime stain profile
      output$dt_csp <- renderDataTable(server = FALSE, {
        datatable(
          rv_file$dt_csp,
          selection = "none",
          options = list(iDisplayLength = 50, ordering = FALSE, autowidth = TRUE),
          rownames = FALSE
        )
      })

      # Display the reference profile
      output$dt_ref <- renderDataTable(server = FALSE, {
        datatable(
          rv_file$dt_ref,
          selection = "none",
          options = list(iDisplayLength = 50, ordering = FALSE, autowidth = TRUE),
          rownames = FALSE
        )
      })

      # Display the allelic ladder
      output$dt_ladder <- renderDataTable(server = FALSE, {
        datatable(
          rv_file$dt_ladder,
          selection = "none",
          options = list(iDisplayLength = 50, ordering = FALSE, autowidth = TRUE),
          rownames = FALSE
        )
      })

      # Display the allele frequency
      output$dt_af <- renderDataTable(server = FALSE, {
        datatable(
          rv_file$dt_af,
          selection = "none",
          options = list(iDisplayLength = 50, ordering = FALSE, autowidth = TRUE),
          rownames = FALSE
        )
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
