#' edit_sr_ui
#'
#' @description The function to create the ui module for the stutter ratio
edit_sr_ui <- function(id){
  ns <- NS(id)

  tabPanel("Edit average stutter ratio",
           titlePanel("Edit average stutter ratio"),
           br(),
           fluidRow(
             column(2,
                    wellPanel(
                      h4("Stutter file"),
                      shinyFilesButton(ns("load_stutter_file"), "Load", "Select a file", multiple = FALSE),
                      br(),
                      br(),
                      textOutput(ns("fn_stutter_file")),
                      br(),
                      uiOutput(ns("select_mk"))
                    ),
                    wellPanel(
                      h4("Calculate stutter ratio"),
                      br(),
                      uiOutput(ns("select_sr_type"))
                    )
             ),
             column(10,
                    tags$style(HTML("
                         .handsontable .htCore td {
                           vertical-align: middle !important;
                           text-align: center !important;
                         }"
                    )),
                    rHandsontableOutput(ns("dt_edit_sr"))
             )
           )
  )
}

#' edit_sr_server
#'
#' @description The function to create the server module for the stutter ratio
edit_sr_server <- function(id, path_pack){
  moduleServer(
    id,
    function(input, output, session){
      rv_edit_sr <- reactiveValues()
      rv_edit_sr$dt_setting_est_sr <- NULL
      rv_edit_sr$list_dt_average_sr <- NULL
      rv_edit_sr$markers_auto <- character(0)
      rv_edit_sr$current_mk <- NULL
      rv_edit_sr$current_sr_types <- character(0)

      shinyFileChoose(input, "load_stutter_file", roots = c(documents = paste0(path_pack, "/extdata/stutter_file")), session = session, filetypes = c("RData"))

      observeEvent(input$load_stutter_file, {
        req(input$load_stutter_file)
        path_stutter_file <- parseFilePaths(roots = c(documents = paste0(path_pack, "/extdata/stutter_file")), input$load_stutter_file)$datapath

        output$fn_stutter_file <- renderText({paste0("File name: ", basename(path_stutter_file))})

        if(length(path_stutter_file) > 0){
          load(path_stutter_file)
          rv_edit_sr$dt_setting_est_sr <- dt_setting_est_sr
          rv_edit_sr$list_dt_average_sr <- list_dt_average_sr
          rv_edit_sr$markers_auto <- names(list_dt_average_sr)
        }
      })

      output$select_mk <- renderUI({
        markers_auto <- rv_edit_sr$markers_auto
        if(length(markers_auto) != 0){
          selectInput(session$ns("select_mk"), label = h4("Marker"), choices = markers_auto, selected = markers_auto[1], selectize = FALSE)
        }else{
          selectInput(session$ns("select_mk"), label = h4("Marker"), choices = markers_auto, selected = NULL, selectize = FALSE)
        }
      })

      observeEvent(input$select_mk, {
        select_mk <- input$select_mk
        if(!is.null(select_mk)){
          # Update edited information
          if(!is.null(rv_edit_sr$current_mk)){
            dt_edit_sr <- hot_to_r(input$dt_edit_sr)
            rv_edit_sr$list_dt_average_sr[[rv_edit_sr$current_mk]] <- dt_edit_sr
          }

          dt_setting_est_sr <- rv_edit_sr$dt_setting_est_sr
          index_mk <- which(dt_setting_est_sr$Marker == select_mk)
          meth_bsr <- dt_setting_est_sr[["BSR model"]][index_mk]
          meth_fsr <- dt_setting_est_sr[["FSR model"]][index_mk]
          meth_dsr <- dt_setting_est_sr[["DSR model"]][index_mk]
          meth_m2sr <- dt_setting_est_sr[["M2SR model"]][index_mk]
          sr_types <- c("BSR", "FSR", "DSR", "M2SR")
          rv_edit_sr$current_sr_types <- sr_types[!is.element(c(meth_bsr, meth_fsr, meth_dsr, meth_m2sr), "Not consider")]

          dt_ave_sr <- rv_edit_sr$list_dt_average_sr[[select_mk]]

          output$dt_edit_sr <- renderRHandsontable({
            rhandsontable(dt_ave_sr, rowHeaders = FALSE, height = 750) %>%
              hot_table(fillHandle = FALSE, rowHeights = 30) %>%
              hot_col(c("Allele"), type = "text", readOnly = TRUE, width = 50) %>%
              hot_col(c("Corrected allele (BSR)", "Corrected allele (FSR)", "Corrected allele (DSR)", "Corrected allele (M2SR)"), type = "numeric", format = "0.000", readOnly = TRUE, width = 200) %>%
              hot_col(c("Average BSR", "Average FSR", "Average DSR", "Average M2SR"), type = "numeric", format = "0.0000", width = 120)
          })
        }

        rv_edit_sr$current_mk <- select_mk
      })

      output$select_sr_type <- renderUI({
        sr_types <- rv_edit_sr$current_sr_types
        if(length(sr_types) != 0){
          selectInput(session$ns("select_sr_type"), label = "SR type", choices = sr_types, selected = sr_types[1], selectize = FALSE)
        }else{
          selectInput(session$ns("select_sr_type"), label = "SR type", choices = sr_types, selected = NULL, selectize = FALSE)
        }
      })
    }
  )
}
