#' file_ui
#'
#' @description The function to create the ui module for loading files
file_ui <- function(id){
  ns <- NS(id)

  tabPanel("File",
           tags$style(HTML("
             .gt_table {
               margin-left: 0 !important;
               margin-right: auto !important;
             }
           ")),
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
                                            br(),
                                            fluidRow(
                                              column(3, uiOutput(ns("select_csp_sn"))),
                                              column(3, uiOutput(ns("select_csp_mk"))),
                                              column(3),
                                              column(3)
                                            ),
                                            div(style = "overflow-x: auto; overflow-y: auto; max-height: calc(100vh - 150px); max-width: 100%;",
                                              gt_output(ns("dt_csp"))
                                            )
                                   ),
                                   tabPanel("Reference profile",
                                            div(style = "overflow-x: auto; overflow-y: auto; max-height: calc(100vh - 150px); max-width: 100%;",
                                              gt_output(ns("dt_ref"))
                                            )
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

      # --- Load crime stain profiles --- #

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

      # --- Load reference profiles --- #

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

      # --- Load allele frequencies --- #

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

      # --- Display crime stain profiles --- #

      output$select_csp_sn <- renderUI({
        dt_csp <- rv_input$dt_csp
        if(!is.null(dt_csp)){
          selectInput(session$ns("select_csp_sn"), label = "Sample name", choices = c("All", unique(dt_csp$SampleName)), selected = "All", selectize = FALSE)
        }
      })

      output$select_csp_mk <- renderUI({
        dt_csp <- rv_input$dt_csp
        if(!is.null(dt_csp)){
          selectInput(session$ns("select_csp_mk"), label = "Marker", choices = c("All", unique(dt_csp$Marker)), selected = "All", selectize = FALSE)
        }
      })

      output$dt_csp <- render_gt({
        dt_csp <- rv_input$dt_csp
        select_csp_sn <- input$select_csp_sn
        select_csp_mk <- input$select_csp_mk

        if(all(!is.null(dt_csp), !is.null(select_csp_sn), !is.null(select_csp_mk))){
          if(select_csp_sn != "All"){
            dt_csp <- dt_csp[dt_csp$SampleName == select_csp_sn]
          }
          if(select_csp_mk != "All"){
            dt_csp <- dt_csp[dt_csp$Marker == select_csp_mk]
          }

          allele_cols <- names(dt_csp)[grep("Allele", names(dt_csp))]
          height_cols <- names(dt_csp)[grep("Height", names(dt_csp))]

          gt(dt_csp) %>%
            sub_missing(
              columns = 1:ncol(dt_csp),
              missing_text = ""
            ) %>%
            tab_spanner(
              label = "Allele",
              columns = allele_cols
            ) %>%
            cols_label(
              !!!setNames(rep("", length(allele_cols)), allele_cols)
            ) %>%
            tab_spanner(
              label = "Height",
              columns = height_cols
            ) %>%
            cols_label(
              !!!setNames(rep("", length(height_cols)), height_cols)
            ) %>%
            tab_options(
              table.font.size = px(14),
              data_row.padding = px(5),
              data_row.padding.horizontal = px(7.5),
              column_labels.padding.horizontal = px(7.5)
            ) %>%
            tab_style(
              style = cell_text(whitespace = "nowrap"), # No string wrapping
              locations = cells_body()
            )
        }
      })

      # --- Display reference profiles --- #

      output$dt_ref <- render_gt({
        dt_ref <- rv_input$dt_ref

        if(!is.null(dt_ref)){
          dt_ref_new <- NULL
          samples <- unique(dt_ref$SampleName)
          n_sample <- length(samples)
          for(i in 1:n_sample){
            dt_tmp <- dt_ref[dt_ref$SampleName == samples[i]]

            if(i == 1){
              dt_ref_new <- dt_tmp[, .(Marker, Allele1, Allele2)]
            }else{
              dt_tmp2 <- dt_tmp[, .(Allele1, Allele2)]
              names(dt_tmp2) <- paste0("Allele", c(2 * i - 1, 2 * i))
              dt_ref_new <- cbind(dt_ref_new, dt_tmp2)
            }
          }

          gt_table <- gt(dt_ref_new)
          for(i in 1:n_sample){
            col_target <- names(dt_ref_new)[c(2 * i, 2 * i + 1)]

            gt_table <- gt_table %>%
              tab_spanner(
                label = samples[i],
                columns = col_target
              ) %>%
              cols_label(
                !!!setNames(c("Allele1", "Allele2"), col_target)
              )
          }
          gt_table %>%
            sub_missing(
              columns = 1:ncol(dt_ref_new),
              missing_text = ""
            ) %>%
            tab_options(
              table.font.size = px(14),
              data_row.padding = px(5),
              data_row.padding.horizontal = px(7.5),
              column_labels.padding.horizontal = px(7.5)
            ) %>%
            tab_style(
              style = cell_text(whitespace = "nowrap"), # No string wrapping
              locations = cells_body()
            )
        }
      })

      # --- Display allele frequencies --- #

      output$dt_af <- renderDataTable({
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
