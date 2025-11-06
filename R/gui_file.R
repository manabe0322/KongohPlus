#' file_ui
#'
#' @description The function to create the ui module for loading files
file_ui <- function(id){
  ns <- NS(id)

  tabPanel(
    "File",
    sidebarLayout(
      sidebarPanel(
        width = 3,
        fileInput(ns("file_csp"), label = h4("Crime stain profile"), accept = ".csv", buttonLabel = "Load"),
        fileInput(ns("file_ref"), label = h4("Reference profile"), accept = ".csv", buttonLabel = "Load"),
        fileInput(ns("file_af"), label = h4("Allele frequency"), accept = ".csv", buttonLabel = "Load"),
        br(),
        actionButton(ns("act_deconvo"), label = "Deconvolution", class = "btn btn-primary btn-lg")
      ),
      mainPanel(
        width = 9,
        tabsetPanel(
          id = ns("tabset_load_file"),
          type = "tabs",
          tabPanel(
            "Crime stain profile",
            uiOutput(ns("info_csp"))
          ),
          tabPanel(
            "Reference profile",
            uiOutput(ns("info_ref"))
          ),
          tabPanel(
            "Allele frequency",
            uiOutput(ns("info_af"))
          )
        )
      )
    )
  )
}

#' file_server
#'
#' @description The function to create the server module for loading files
file_server <- function(id, session_top, bool_load_proj){
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
      rv_input$bool_csp <- TRUE
      rv_input$bool_ref <- TRUE
      rv_input$bool_af <- TRUE

      # --- Load crime stain profiles --- #

      observeEvent(input$file_csp, {
        file_input <- input$file_csp
        if(!is.null(file_input)){
          dt_csp <- fread(file_input$datapath)
          dt_csp <- try(change_data_type_csp(dt_csp))
          if(!inherits(dt_csp, "try-error")){
            rv_input$dt_csp <- dt_csp
            rv_input$fn_csp <- file_input$name
            rv_input$bool_csp <- TRUE
          }else{
            showModal(modalDialog(title = "Error", "The selected file format is invalid. Please upload a file in the correct format.", easyClose = TRUE, footer = NULL))
            rv_input$dt_csp <- NULL
            rv_input$fn_csp <- NULL
            rv_input$bool_csp <- FALSE
          }
          updateTabsetPanel(inputId = "tabset_load_file", selected = "Crime stain profile")
        }
      })

      # --- Load reference profiles --- #

      observeEvent(input$file_ref, {
        file_input <- input$file_ref
        if(!is.null(file_input)){
          dt_ref <- fread(file_input$datapath)
          dt_ref <- try(change_data_type_ref(dt_ref))
          if(!inherits(dt_ref, "try-error")){
            rv_input$dt_ref <- dt_ref
            rv_input$fn_ref <- file_input$name
            rv_input$bool_ref <- TRUE
          }else{
            showModal(modalDialog(title = "Error", "The selected file format is invalid. Please upload a file in the correct format.", easyClose = TRUE, footer = NULL))
            rv_input$dt_ref <- NULL
            rv_input$fn_ref <- NULL
            rv_input$bool_ref <- FALSE
          }
          updateTabsetPanel(inputId = "tabset_load_file", selected = "Reference profile")
        }
      })

      # --- Load allele frequencies --- #

      observeEvent(input$file_af, {
        file_input <- input$file_af
        if(!is.null(file_input)){
          dt_af <- fread(file_input$datapath)
          dt_af <- try(change_data_type_af(dt_af))
          if(!inherits(dt_af, "try-error")){
            rv_input$dt_af <- dt_af
            rv_input$fn_af <- file_input$name
            rv_input$bool_af <- TRUE
          }else{
            showModal(modalDialog(title = "Error", "The selected file format is invalid. Please upload a file in the correct format.", easyClose = TRUE, footer = NULL))
            rv_input$dt_af <- NULL
            rv_input$fn_af <- NULL
            rv_input$bool_af <- FALSE
          }
          updateTabsetPanel(inputId = "tabset_load_file", selected = "Allele frequency")
        }
      })

      # --- Display crime stain profiles --- #

      output$info_csp <- renderUI({
        if(rv_input$bool_csp){
          dataTableOutput(session$ns("dt_csp"))
        }

        else{
          tagList(
            br(),
            tags$p("The selected file format is invalid. Please upload a file in the correct format.", style = "color: red; font-weight: bold;")
          )
        }
      })

      output$dt_csp <- renderDT({
        dt_csp <- rv_input$dt_csp

        if(!is.null(dt_csp)){
          dt_display_csp <- create_dt_display_csp(dt_csp)
          sn_csp_all <- unique(dt_csp[, SampleName])
          n_rep <- length(sn_csp_all)
          group_names <- paste0("Replicate", 1:n_rep, ": ", sn_csp_all)
          groups <- define_rep_groups(group_names)
          sketch <- create_sketch_csp(dt_display_csp, groups)

          datatable(
            dt_display_csp,
            container = sketch,
            selection = "none",
            rownames = FALSE,
            extensions = c("FixedColumns", "FixedHeader", "Scroller"),
            options = list(
              ordering = FALSE,
              autoWidth = FALSE,
              paging = FALSE,
              scrollX = TRUE,
              scrollY = "650px",
              fixedHeader = FALSE,
              fixedColumns = list(leftColumns = 2),
              columnDefs = list(
                list(className = "dt-right", targets = 1),
                list(width = '100px', targets = 0),
                list(width = '100px', targets = 1)
              )
            )
          ) %>%
            formatStyle(names(dt_display_csp), lineHeight = "80%")
        }
      })

      # --- Display reference profiles --- #

      output$info_ref <- renderUI({
        if(rv_input$bool_ref){
          dataTableOutput(session$ns("dt_ref"))
        }

        else{
          tagList(
            br(),
            tags$p("The selected file format is invalid. Please upload a file in the correct format.", style = "color: red; font-weight: bold;")
          )
        }
      })

      output$dt_ref <- renderDataTable({
        dt_ref <- rv_input$dt_ref

        if(!is.null(dt_ref)){
          dt_ref_display <- NULL
          samples <- unique(dt_ref$SampleName)
          n_sample <- length(samples)
          for(i in 1:n_sample){
            dt_tmp <- dt_ref[dt_ref$SampleName == samples[i]]

            gts <- rep(NA, nrow(dt_tmp))
            for(j in 1:nrow(dt_tmp)){
              gt1 <- c(dt_tmp$Allele1[j], dt_tmp$Allele2[j])
              gt1 <- gt1[!is.na(gt1)]
              gts[j] <- paste0(gt1, collapse = ", ")
            }

            dt_tmp2 <- data.table(Marker = dt_tmp$Marker,
                                  Gt = gts)
            names(dt_tmp2) <- c("Marker", samples[i])

            if(i == 1){
              dt_ref_display <- dt_tmp2
            }else{
              dt_ref_display <- full_join(dt_ref_display, dt_tmp2, by = c("Marker"))
            }
          }

          datatable(
            dt_ref_display,
            selection = "none",
            extensions = c("FixedColumns", "FixedHeader"),
            options = list(dom = "t", iDisplayLength = nrow(dt_ref_display), ordering = FALSE, autoWidth = FALSE, scrollX = TRUE, scrollY = "calc(100vh - 200px)",
                           fixedHeader = FALSE, fixedColumns = list(leftColumns = 1, rightColumns = 0),
                           columnDefs = list(list(targets = 1:(ncol(dt_ref_display) - 1), searchable = FALSE))),
            rownames = FALSE
          ) %>%
            formatStyle(names(dt_ref_display), lineHeight = "80%")
        }
      })

      # --- Display allele frequencies --- #

      output$info_af <- renderUI({
        if(rv_input$bool_af){
          dataTableOutput(session$ns("dt_af"))
        }

        else{
          tagList(
            br(),
            tags$p("The selected file format is invalid. Please upload a file in the correct format.", style = "color: red; font-weight: bold;")
          )
        }
      })

      output$dt_af <- renderDataTable({
        dt_af <- rv_input$dt_af
        if(!is.null(dt_af)){

          datatable(
            dt_af,
            selection = "none",
            extensions = c("FixedColumns", "FixedHeader"),
            options = list(dom = "t", iDisplayLength = nrow(dt_af), ordering = FALSE, autoWidth = FALSE, scrollX = TRUE, scrollY = "calc(100vh - 150px)",
                           fixedHeader = FALSE, fixedColumns = list(leftColumns = 1, rightColumns = 0)),
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
