#' file_ui
#'
#' @description The function to create the ui module for loading files
file_ui <- function(id){
  ns <- NS(id)

  tabPanel("File",
           sidebarLayout(
             sidebarPanel(width = 3,
                          fileInput(ns("file_csp"), label = h4("Crime stain profile"), accept = ".csv", buttonLabel = "Load"),
                          fileInput(ns("file_ref"), label = h4("Reference profile"), accept = ".csv", buttonLabel = "Load"),
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
                                   tabPanel("Allele frequency",
                                            dataTableOutput(ns("dt_af"))
                                   )
                       )
             )
           )
  )
}

#' file_server
#'
#' @description The function to create the server module for loading files
file_server <- function(id, session_top){
  moduleServer(
    id,
    function(input, output, session){
      rv_file <- reactiveValues()
      rv_file$dt_csp <- NULL
      rv_file$dt_ref <- NULL
      rv_file$dt_af <- NULL
      rv_file$fn_csp <- NULL
      rv_file$fn_ref <- NULL
      rv_file$fn_af <- NULL

      # Load the crime stain profile
      observeEvent(input$file_csp, {
        file_input <- input$file_csp
        if(!is.null(file_input)){
          dt_csp <- fread(file_input$datapath)
          dt_csp <- change_data_type_csp(dt_csp)
          rv_file$dt_csp <- dt_csp
          rv_file$fn_csp <- file_input$name
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
          rv_file$fn_ref <- file_input$name
          updateTabsetPanel(inputId = "tabset_load_file", selected = "Reference profile")
        }
      })

      # Load the allele frequency
      observeEvent(input$file_af, {
        file_input <- input$file_af
        if(!is.null(file_input)){
          dt_af <- fread(file_input$datapath)
          dt_af <- change_data_type_af(dt_af)
          rv_file$dt_af <- dt_af
          rv_file$fn_af <- file_input$name
          updateTabsetPanel(inputId = "tabset_load_file", selected = "Allele frequency")
        }
      })

      # Display the crime stain profile
      output$dt_csp <- renderDataTable(server = FALSE, {
        dt_csp <- rv_file$dt_csp
        if(!is.null(dt_csp)){
          dt_csp_display <- copy(dt_csp)
          dt_csp_display$SampleName <- as.factor(dt_csp_display$SampleName)
          dt_csp_display$Marker <- as.factor(dt_csp_display$Marker)
          name_dt_csp <- names(dt_csp_display)
          index_sn <- which(name_dt_csp == "SampleName") - 1
          index_mk <- which(name_dt_csp == "Marker") - 1
          index_search_disable <- setdiff(0:(ncol(dt_csp_display) - 1), c(index_sn, index_mk))

          datatable(
            dt_csp_display,
            filter = "top",
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_csp_display), ordering = FALSE, autoWidth = TRUE, scrollX = TRUE, scrollY = "calc(100vh - 200px)",
                           columnDefs = list(list(targets = index_search_disable, searchable = FALSE, width = "40px"),
                                             list(targets = c(index_sn, index_mk), width = "100px"))),
            rownames = FALSE
          )
        }
      })

      # Display the reference profile
      output$dt_ref <- renderDataTable(server = FALSE, {
        dt_ref <- rv_file$dt_ref
        if(!is.null(dt_ref)){
          dt_ref_display <- copy(dt_ref)
          dt_ref_display$SampleName <- as.factor(dt_ref_display$SampleName)
          dt_ref_display$Marker <- as.factor(dt_ref_display$Marker)
          name_dt_ref <- names(dt_ref_display)
          index_sn <- which(name_dt_ref == "SampleName") - 1
          index_mk <- which(name_dt_ref == "Marker") - 1
          index_search_disable <- setdiff(0:(ncol(dt_ref_display) - 1), c(index_sn, index_mk))

          datatable(
            dt_ref_display,
            filter = "top",
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_ref_display), ordering = FALSE, autoWidth = TRUE, scrollX = TRUE, scrollY = "calc(100vh - 200px)",
                           columnDefs = list(list(targets = index_search_disable, searchable = FALSE, width = "40px"),
                                             list(targets = c(index_sn, index_mk), width = "100px"))),
            rownames = FALSE
          )
        }
      })

      # Display the allele frequency
      output$dt_af <- renderDataTable(server = FALSE, {
        dt_af <- rv_file$dt_af
        datatable(
          dt_af,
          selection = "none",
          extensions = c("FixedColumns", "FixedHeader"),
          options = list(dom = "t", iDisplayLength = nrow(dt_af), ordering = FALSE, autoWidth = TRUE, scrollX = TRUE, scrollY = "calc(100vh - 200px)",
                         fixedHeader = TRUE, fixedColumns = list(leftColumns = 1, rightColumns = 0)),
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
