#' edit_sr_ui
#'
#' @description The function to create the ui module for the stutter ratio
edit_sr_ui <- function(id){
  ns <- NS(id)

  tabPanel("4. Edit stutter ratios",
           titlePanel("4. Edit stutter ratios"),
           br(),
           fluidRow(
             column(3,
                    wellPanel(
                      h4("Stutter data"),
                      fluidRow(
                        column(3, br(), shinyFilesButton(ns("load_stutter_data"), label = "Load", title = "Select a file", multiple = FALSE)),
                        column(9, uiOutput(ns("fn_stutter_data")))
                      ),
                      h4("Marker"),
                      uiOutput(ns("select_mk")),
                      h4("Add an allele"),
                      fluidRow(
                        column(9, disabled(numericInput(ns("add_al"), label = NULL, value = NULL))),
                        column(3, disabled(actionButton(ns("act_add_al"), label = "Add")))
                      ),
                      h4("Delete an allele"),
                      fluidRow(
                        column(9, uiOutput(ns("select_del_al"))),
                        column(3, disabled(actionButton(ns("act_del_al"), label = "Delete")))
                      ),
                      br(),
                      br(),
                      disabled(actionButton(ns("act_save_as_edit_sr"), label = "Save as", class = "btn btn-primary btn-lg")),
                      br(),
                      br(),
                      disabled(actionButton(ns("act_save_edit_sr"), label = "Save", class = "btn btn-primary btn-lg"))
                    )
             ),
             column(9,
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

replace_ave_sr <- function(list_dt_average_sr){
  for(i in 1:length(list_dt_average_sr)){
    dt_ave_sr <- list_dt_average_sr[[i]]

    # Replace NA to 0
    dt_ave_sr[["Average BSR"]][is.na(dt_ave_sr[["Average BSR"]])] <- 0
    dt_ave_sr[["Average FSR"]][is.na(dt_ave_sr[["Average FSR"]])] <- 0
    dt_ave_sr[["Average DSR"]][is.na(dt_ave_sr[["Average DSR"]])] <- 0
    dt_ave_sr[["Average M2SR"]][is.na(dt_ave_sr[["Average M2SR"]])] <- 0

    # Replace <0 to 0
    dt_ave_sr[["Average BSR"]][dt_ave_sr[["Average BSR"]] < 0] <- 0
    dt_ave_sr[["Average FSR"]][dt_ave_sr[["Average FSR"]] < 0] <- 0
    dt_ave_sr[["Average DSR"]][dt_ave_sr[["Average DSR"]] < 0] <- 0
    dt_ave_sr[["Average M2SR"]][dt_ave_sr[["Average M2SR"]] < 0] <- 0

    list_dt_average_sr[[i]] <- dt_ave_sr
  }
  return(list_dt_average_sr)
}

#' edit_sr_server
#'
#' @description The function to create the server module for the stutter ratio
edit_sr_server <- function(id, path_pack){
  moduleServer(
    id,
    function(input, output, session){
      rv_edit_sr <- reactiveValues()
      rv_edit_sr$fn_stutter_data <- NULL
      rv_edit_sr$dt_setting_est_sr <- NULL
      rv_edit_sr$dt_result_est_bsr <- NULL
      rv_edit_sr$dt_result_est_fsr <- NULL
      rv_edit_sr$dt_result_est_dsr <- NULL
      rv_edit_sr$dt_result_est_m2sr <- NULL
      rv_edit_sr$list_dt_average_sr <- NULL
      rv_edit_sr$markers_auto <- character(0)
      rv_edit_sr$current_mk <- NULL

      shinyFileChoose(input, "load_stutter_data", roots = c(documents = paste0(path_pack, "/extdata/stutter")), session = session, filetypes = c("RData"))

      observeEvent(input$load_stutter_data, {
        req(input$load_stutter_data)
        path_stutter_data <- parseFilePaths(roots = c(documents = paste0(path_pack, "/extdata/stutter")), input$load_stutter_data)$datapath

        if(length(path_stutter_data) > 0){
          rv_edit_sr$fn_stutter_data <- fn_stutter_data <- basename(path_stutter_data)

          rv_edit_sr$dt_setting_est_sr <- NULL
          rv_edit_sr$dt_result_est_bsr <- NULL
          rv_edit_sr$dt_result_est_fsr <- NULL
          rv_edit_sr$dt_result_est_dsr <- NULL
          rv_edit_sr$dt_result_est_m2sr <- NULL
          rv_edit_sr$list_dt_average_sr <- NULL
          rv_edit_sr$markers_auto <- character(0)
          rv_edit_sr$current_mk <- NULL

          load(path_stutter_data)
          rv_edit_sr$dt_setting_est_sr <- dt_setting_est_sr
          rv_edit_sr$dt_result_est_bsr <- dt_result_est_bsr
          rv_edit_sr$dt_result_est_fsr <- dt_result_est_fsr
          rv_edit_sr$dt_result_est_dsr <- dt_result_est_dsr
          rv_edit_sr$dt_result_est_m2sr <- dt_result_est_m2sr
          rv_edit_sr$list_dt_average_sr <- list_dt_average_sr
          rv_edit_sr$markers_auto <- names(list_dt_average_sr)

          select_mk <- rv_edit_sr$markers_auto[1]

          dt_setting_est_sr <- rv_edit_sr$dt_setting_est_sr
          index_mk <- which(dt_setting_est_sr$Marker == select_mk)
          meth_bsr <- dt_setting_est_sr[["BSR model"]][index_mk]
          meth_fsr <- dt_setting_est_sr[["FSR model"]][index_mk]
          meth_dsr <- dt_setting_est_sr[["DSR model"]][index_mk]
          meth_m2sr <- dt_setting_est_sr[["M2SR model"]][index_mk]
          sr_types <- c("BSR", "FSR", "DSR", "M2SR")

          dt_ave_sr <- rv_edit_sr$list_dt_average_sr[[select_mk]]

          output$dt_edit_sr <- renderRHandsontable({
            rhandsontable(dt_ave_sr, rowHeaders = FALSE, height = 600) %>%
              hot_table(fillHandle = FALSE, rowHeights = 30) %>%
              hot_col(c("Allele"), type = "text", readOnly = TRUE, width = 50) %>%
              hot_col(c("Corrected allele (BSR)", "Corrected allele (FSR)", "Corrected allele (DSR)", "Corrected allele (M2SR)"), type = "numeric", format = "0.000", readOnly = TRUE, width = 175) %>%
              hot_col(c("Average BSR", "Average FSR", "Average DSR", "Average M2SR"), type = "numeric", format = "0.0000", width = 120)
          })

          rv_edit_sr$current_mk <- select_mk

          enable("add_al")
          enable("act_add_al")
          enable("act_del_al")
          enable("act_save_as_edit_sr")
          enable("act_save_edit_sr")
        }else{
          disable("add_al")
          disable("act_add_al")
          disable("act_del_al")
          disable("act_save_as_edit_sr")
          disable("act_save_edit_sr")
        }
      })

      output$fn_stutter_data <- renderUI({
        disabled(textInput(session$ns("fn_stutter_data"), label = "", value = rv_edit_sr$fn_stutter_data))
      })

      output$select_mk <- renderUI({
        markers_auto <- rv_edit_sr$markers_auto
        if(length(markers_auto) != 0){
          selectInput(session$ns("select_mk"), label = NULL, choices = markers_auto, selected = markers_auto[1], selectize = FALSE)
        }else{
          disabled(selectInput(session$ns("select_mk"), label = NULL, choices = markers_auto, selectize = FALSE))
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

          dt_ave_sr <- rv_edit_sr$list_dt_average_sr[[select_mk]]

          output$dt_edit_sr <- renderRHandsontable({
            rhandsontable(dt_ave_sr, rowHeaders = FALSE, height = 600) %>%
              hot_table(fillHandle = FALSE, rowHeights = 30) %>%
              hot_col(c("Allele"), type = "text", readOnly = TRUE, width = 50) %>%
              hot_col(c("Corrected allele (BSR)", "Corrected allele (FSR)", "Corrected allele (DSR)", "Corrected allele (M2SR)"), type = "numeric", format = "0.000", readOnly = TRUE, width = 175) %>%
              hot_col(c("Average BSR", "Average FSR", "Average DSR", "Average M2SR"), type = "numeric", format = "0.0000", width = 120)
          })
        }

        rv_edit_sr$current_mk <- select_mk
      })

      ## Add an allele ##

      observeEvent(input$act_add_al, {
        add_al <- input$add_al
        if(!is.na(add_al)){
          current_mk <- rv_edit_sr$current_mk
          dt_ave_sr <- rv_edit_sr$list_dt_average_sr[[current_mk]]
          dt_ave_sr <- rbind(dt_ave_sr, list(as.character(add_al), NA, 0, NA, 0, NA, 0, NA, 0))
          tmp_al <- as.numeric(dt_ave_sr$Allele)
          tmp_al_sort <- sort(tmp_al)
          dt_ave_sr <- dt_ave_sr[match(tmp_al_sort, tmp_al)]
          rv_edit_sr$list_dt_average_sr[[current_mk]] <- dt_ave_sr

          output$dt_edit_sr <- renderRHandsontable({
            rhandsontable(dt_ave_sr, rowHeaders = FALSE, height = 600) %>%
              hot_table(fillHandle = FALSE, rowHeights = 30) %>%
              hot_col(c("Allele"), type = "text", readOnly = TRUE, width = 50) %>%
              hot_col(c("Corrected allele (BSR)", "Corrected allele (FSR)", "Corrected allele (DSR)", "Corrected allele (M2SR)"), type = "numeric", format = "0.000", readOnly = TRUE, width = 175) %>%
              hot_col(c("Average BSR", "Average FSR", "Average DSR", "Average M2SR"), type = "numeric", format = "0.0000", width = 120)
          })
        }
      })

      ## Delete an allele ##

      output$select_del_al <- renderUI({
        current_mk <- rv_edit_sr$current_mk
        if(!is.null(current_mk)){
          dt_ave_sr <- rv_edit_sr$list_dt_average_sr[[current_mk]]
          current_alleles <- dt_ave_sr$Allele
          selectInput(session$ns("del_al"), label = NULL, choices = current_alleles, selected = current_alleles[1], selectize = FALSE)
        }else{
          disabled(selectInput(session$ns("del_al"), label = NULL, choices = numeric(0), selectize = FALSE))
        }
      })

      observeEvent(input$act_del_al, {
        del_al <- input$del_al
        current_mk <- rv_edit_sr$current_mk
        dt_ave_sr <- rv_edit_sr$list_dt_average_sr[[current_mk]]
        current_alleles <- dt_ave_sr$Allele
        dt_ave_sr <- dt_ave_sr[- which(current_alleles == del_al)]
        rv_edit_sr$list_dt_average_sr[[current_mk]] <- dt_ave_sr

        output$dt_edit_sr <- renderRHandsontable({
          rhandsontable(dt_ave_sr, rowHeaders = FALSE, height = 600) %>%
            hot_table(fillHandle = FALSE, rowHeights = 30) %>%
            hot_col(c("Allele"), type = "text", readOnly = TRUE, width = 50) %>%
            hot_col(c("Corrected allele (BSR)", "Corrected allele (FSR)", "Corrected allele (DSR)", "Corrected allele (M2SR)"), type = "numeric", format = "0.000", readOnly = TRUE, width = 175) %>%
            hot_col(c("Average BSR", "Average FSR", "Average DSR", "Average M2SR"), type = "numeric", format = "0.0000", width = 120)
        })
      })

      ## Save as ##

      observeEvent(input$act_save_as_edit_sr, {
        showModal(modalDialog(
          title = "Save as",
          textInput(session$ns("fn_edit_sr"), label = "Enter the file name"),
          br(),
          "Note: Blank or negative values in SR cells will be automatically replaced with zero.",
          footer = tagList(
            actionButton(session$ns("act_save_as_edit_sr_2"), "Save"),
            modalButton("Close")
          )
        ))
      })

      observeEvent(input$fn_edit_sr, {
        fn_edit_sr <- input$fn_edit_sr
        fn_stutter_all <- list.files(paste0(path_pack, "/extdata/stutter"))
        if(is.element(paste0(fn_edit_sr, ".RData"), fn_stutter_all)){
          showFeedbackDanger(inputId = "fn_edit_sr", text = "A file with the same name already exists.")
          disable("act_save_as_edit_sr_2")
        }else{
          hideFeedback("fn_edit_sr")
          enable("act_save_as_edit_sr_2")
        }
      })

      observeEvent(input$act_save_as_edit_sr_2, {
        # Update edited information
        if(!is.null(rv_edit_sr$current_mk)){
          dt_edit_sr <- hot_to_r(input$dt_edit_sr)
          rv_edit_sr$list_dt_average_sr[[rv_edit_sr$current_mk]] <- dt_edit_sr
        }

        dt_setting_est_sr <- rv_edit_sr$dt_setting_est_sr
        dt_result_est_bsr <- rv_edit_sr$dt_result_est_bsr
        dt_result_est_fsr <- rv_edit_sr$dt_result_est_fsr
        dt_result_est_dsr <- rv_edit_sr$dt_result_est_dsr
        dt_result_est_m2sr <- rv_edit_sr$dt_result_est_m2sr
        list_dt_average_sr <- rv_edit_sr$list_dt_average_sr

        list_dt_average_sr <- replace_ave_sr(list_dt_average_sr)

        save(list = c("dt_setting_est_sr",
                      "dt_result_est_bsr", "dt_result_est_fsr", "dt_result_est_dsr", "dt_result_est_m2sr",
                      "list_dt_average_sr"),
             file = paste0(path_pack, "/extdata/stutter/", input$fn_edit_sr, ".RData"))

        removeModal()
      })

      ## Save ##

      observeEvent(input$act_save_edit_sr, {
        showModal(modalDialog(
          title = "Save",
          paste0("The existing file '", rv_edit_sr$fn_stutter_data, "' will be overwritten."),
          br(),
          "Press OK to continue.",
          br(),
          br(),
          "Note: Blank or negative values in SR cells will be automatically replaced with zero.",
          footer = tagList(
            actionButton(session$ns("act_save_edit_sr_2"), "OK"),
            modalButton("Close")
          )
        ))
      })

      observeEvent(input$act_save_edit_sr_2, {
        # Update edited information
        if(!is.null(rv_edit_sr$current_mk)){
          dt_edit_sr <- hot_to_r(input$dt_edit_sr)
          rv_edit_sr$list_dt_average_sr[[rv_edit_sr$current_mk]] <- dt_edit_sr
        }

        dt_setting_est_sr <- rv_edit_sr$dt_setting_est_sr
        dt_result_est_bsr <- rv_edit_sr$dt_result_est_bsr
        dt_result_est_fsr <- rv_edit_sr$dt_result_est_fsr
        dt_result_est_dsr <- rv_edit_sr$dt_result_est_dsr
        dt_result_est_m2sr <- rv_edit_sr$dt_result_est_m2sr
        list_dt_average_sr <- rv_edit_sr$list_dt_average_sr

        list_dt_average_sr <- replace_ave_sr(list_dt_average_sr)

        save(list = c("dt_setting_est_sr",
                      "dt_result_est_bsr", "dt_result_est_fsr", "dt_result_est_dsr", "dt_result_est_m2sr",
                      "list_dt_average_sr"),
             file = paste0(path_pack, "/extdata/stutter/", rv_edit_sr$fn_stutter_data))
      })
    }
  )
}
