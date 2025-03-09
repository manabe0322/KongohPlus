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
                          actionButton(ns("act_deconvo"), label = "Deconvolution", class = "btn btn-primary btn-lg")
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
      rv_input <- reactiveValues()
      rv_input$dt_csp <- NULL
      rv_input$dt_ref <- NULL
      rv_input$dt_af <- NULL
      rv_input$fn_csp <- NULL
      rv_input$fn_ref <- NULL
      rv_input$fn_af <- NULL
      rv_input$hnc_from <- 1
      rv_input$hnc_to <- 1
      rv_input$fn_anameth <- ""

      # Load the crime stain profile
      observeEvent(input$file_csp, {
        file_input <- input$file_csp
        if(!is.null(file_input)){
          dt_csp <- fread(file_input$datapath)
          dt_csp <- change_data_type_csp(dt_csp)
          rv_input$dt_csp <- dt_csp
          rv_input$fn_csp <- file_input$name
          updateTabsetPanel(inputId = "tabset_load_file", selected = "Crime stain profile")
        }
      })

      # Load the reference profile
      observeEvent(input$file_ref, {
        file_input <- input$file_ref
        if(!is.null(file_input)){
          dt_ref <- fread(file_input$datapath)
          dt_ref <- change_data_type_ref(dt_ref)
          rv_input$dt_ref <- dt_ref
          rv_input$fn_ref <- file_input$name
          updateTabsetPanel(inputId = "tabset_load_file", selected = "Reference profile")
        }
      })

      # Load the allele frequency
      observeEvent(input$file_af, {
        file_input <- input$file_af
        if(!is.null(file_input)){
          dt_af <- fread(file_input$datapath)
          dt_af <- change_data_type_af(dt_af)
          rv_input$dt_af <- dt_af
          rv_input$fn_af <- file_input$name
          updateTabsetPanel(inputId = "tabset_load_file", selected = "Allele frequency")
        }
      })

      # Display the crime stain profile
      output$dt_csp <- renderDataTable(server = FALSE, {
        dt_csp <- rv_input$dt_csp

        if(!is.null(dt_csp)){
          name_dt_csp <- names(dt_csp)
          index_sn <- which(name_dt_csp == "SampleName") - 1
          index_mk <- which(name_dt_csp == "Marker") - 1
          index_num <- setdiff(0:(ncol(dt_csp) - 1), c(index_sn, index_mk))

          datatable(
            dt_csp,
            selection = "none",
            extensions = c("FixedColumns", "FixedHeader"),
            options = list(dom = "t", iDisplayLength = nrow(dt_csp), ordering = FALSE, autoWidth = TRUE, scrollX = TRUE, scrollY = "calc(100vh - 200px)",
                           fixedHeader = TRUE, fixedColumns = list(leftColumns = 1, rightColumns = 0),
                           columnDefs = list(list(targets = index_num, width = "40px"),
                                             list(targets = c(index_sn, index_mk), width = "100px"))),
            rownames = FALSE
          ) %>%
            formatStyle(names(dt_csp), lineHeight = "80%")
        }
      })

      # Display the reference profile
      output$dt_ref <- renderDataTable(server = FALSE, {
        dt_ref <- rv_input$dt_ref

        if(!is.null(dt_ref)){
          name_dt_ref <- names(dt_ref)
          index_sn <- which(name_dt_ref == "SampleName") - 1
          index_mk <- which(name_dt_ref == "Marker") - 1
          index_num <- setdiff(0:(ncol(dt_ref) - 1), c(index_sn, index_mk))

          datatable(
            dt_ref,
            selection = "none",
            extensions = c("FixedColumns", "FixedHeader"),
            options = list(dom = "t", iDisplayLength = nrow(dt_ref), ordering = FALSE, autoWidth = TRUE, scrollX = TRUE, scrollY = "calc(100vh - 200px)",
                           fixedHeader = TRUE, fixedColumns = list(leftColumns = 1, rightColumns = 0),
                           columnDefs = list(list(targets = index_num, width = "40px"),
                                             list(targets = c(index_sn, index_mk), width = "100px"))),
            rownames = FALSE
          ) %>%
            formatStyle(names(dt_ref), lineHeight = "80%")
        }
      })

      # Display the allele frequency
      output$dt_af <- renderDataTable(server = FALSE, {
        dt_af <- rv_input$dt_af
        if(!is.null(dt_af)){
          datatable(
            dt_af,
            selection = "none",
            extensions = c("FixedColumns", "FixedHeader"),
            options = list(dom = "t", iDisplayLength = nrow(dt_af), ordering = FALSE, autoWidth = TRUE, scrollX = TRUE, scrollY = "calc(100vh - 150px)",
                           fixedHeader = TRUE, fixedColumns = list(leftColumns = 1, rightColumns = 0)),
            rownames = FALSE
          ) %>%
            formatStyle(names(dt_af), lineHeight = "80%")
        }
      })

      observeEvent(input$act_deconvo, {
        updateNavbarPage(session = session_top, inputId = "navbar", selected = "Deconvolution")
      })

      return(rv_input)
    }
  )
}
