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
                      actionButton(ns("act_add_al"), label = "Add"),
                      br(),
                      br(),
                      actionButton(ns("act_del_al"), label = "Delete")
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

      ## Add an allele ##

      observeEvent(input$act_add_al, {
        showModal(modalDialog(
          title = "Add an allele",
          numericInput(session$ns("add_al"), label = "Allele", value = NULL),
          br(),
          checkboxInput(session$ns("ck_auto_calc_sr"), label = "Calculate stutter ratios", value = FALSE, width = NULL),
          disabled(numericInput(session$ns("alcor_bsr"), label = "Corrected allele for BSR", value = NULL)),
          disabled(numericInput(session$ns("alcor_fsr"), label = "Corrected allele for FSR", value = NULL)),
          disabled(numericInput(session$ns("alcor_dsr"), label = "Corrected allele for DSR", value = NULL)),
          disabled(numericInput(session$ns("alcor_m2sr"), label = "Corrected allele for M2SR", value = NULL)),
          footer = tagList(
            actionButton(session$ns("act_add_al_done"), "Done"),
            modalButton("Cancel")
          ),
          size = "l"
        ))
      })

      observeEvent(input$add_al, {
        add_al <- input$add_al
        if(!is.numeric(add_al)){
          showFeedbackDanger(inputId = "add_al", text = "A positive number is required!")
        }else if(add_al <= 0){
          showFeedbackDanger(inputId = "add_al", text = "A positive number is required!")
        }else{
          hideFeedback("add_al")
        }
      })

      observeEvent(input$ck_auto_calc_sr, {
        ck_auto_calc_sr <- input$ck_auto_calc_sr
        if(ck_auto_calc_sr){
          enable("alcor_bsr")
          enable("alcor_fsr")
          enable("alcor_dsr")
          enable("alcor_m2sr")
        }else{
          #updateNumericInput(session, "alcor_bsr", value = NULL)
          #updateNumericInput(session, "alcor_fsr", value = NULL)
          #updateNumericInput(session, "alcor_dsr", value = NULL)
          #updateNumericInput(session, "alcor_m2sr", value = NULL)
          disable("alcor_bsr")
          disable("alcor_fsr")
          disable("alcor_dsr")
          disable("alcor_m2sr")
        }
      })

      observeEvent(input$act_add_al_done, {
        current_mk <- rv_edit_sr$current_mk
        dt_ave_sr <- rv_edit_sr$list_dt_average_sr[[current_mk]]
        dt_ave_sr <- rbind(dt_ave_sr, list(as.character(input$add_al), NA, NA, NA, NA, NA, NA, NA, NA))
        tmp_al <- as.numeric(dt_ave_sr$Allele)
        tmp_al_sort <- sort(tmp_al)
        dt_ave_sr <- dt_ave_sr[match(tmp_al_sort, tmp_al)]
        rv_edit_sr$list_dt_average_sr[[current_mk]] <- dt_ave_sr

        output$dt_edit_sr <- renderRHandsontable({
          rhandsontable(dt_ave_sr, rowHeaders = FALSE, height = 750) %>%
            hot_table(fillHandle = FALSE, rowHeights = 30) %>%
            hot_col(c("Allele"), type = "text", readOnly = TRUE, width = 50) %>%
            hot_col(c("Corrected allele (BSR)", "Corrected allele (FSR)", "Corrected allele (DSR)", "Corrected allele (M2SR)"), type = "numeric", format = "0.000", readOnly = TRUE, width = 200) %>%
            hot_col(c("Average BSR", "Average FSR", "Average DSR", "Average M2SR"), type = "numeric", format = "0.0000", width = 120)
        })

        removeModal()
      })

      ## Delete an allele ##

      observeEvent(input$act_del_al, {
        current_mk <- rv_edit_sr$current_mk
        dt_ave_sr <- rv_edit_sr$list_dt_average_sr[[current_mk]]
        current_alleles <- dt_ave_sr$Allele

        showModal(modalDialog(
          title = "Delete an allele",
          selectInput(session$ns("del_al"), label = h4("Allele"), choices = current_alleles, selected = current_alleles[1], selectize = FALSE),
          footer = tagList(
            actionButton(session$ns("act_del_al_done"), "Done"),
            modalButton("Cancel")
          ),
          size = "l"
        ))
      })

      observeEvent(input$act_del_al_done, {
        del_al <- input$del_al
        current_mk <- rv_edit_sr$current_mk
        dt_ave_sr <- rv_edit_sr$list_dt_average_sr[[current_mk]]
        current_alleles <- dt_ave_sr$Allele
        dt_ave_sr <- dt_ave_sr[- which(current_alleles == del_al)]
        rv_edit_sr$list_dt_average_sr[[current_mk]] <- dt_ave_sr

        output$dt_edit_sr <- renderRHandsontable({
          rhandsontable(dt_ave_sr, rowHeaders = FALSE, height = 750) %>%
            hot_table(fillHandle = FALSE, rowHeights = 30) %>%
            hot_col(c("Allele"), type = "text", readOnly = TRUE, width = 50) %>%
            hot_col(c("Corrected allele (BSR)", "Corrected allele (FSR)", "Corrected allele (DSR)", "Corrected allele (M2SR)"), type = "numeric", format = "0.000", readOnly = TRUE, width = 200) %>%
            hot_col(c("Average BSR", "Average FSR", "Average DSR", "Average M2SR"), type = "numeric", format = "0.0000", width = 120)
        })

        removeModal()
      })
    }
  )
}
