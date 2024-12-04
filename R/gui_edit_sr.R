#' edit_sr_ui
#'
#' @description The function to create the ui module for the stutter ratio
edit_sr_ui <- function(id){
  ns <- NS(id)

  tabPanel("Edit average stutter ratio",
           titlePanel("Edit average stutter ratio"),
           br(),
           fluidRow(
             column(1,
                    shinyFilesButton(ns("load_stutter_file"), label = "Load stutter file", title = "Select a file", multiple = FALSE, class = "btn btn-primary"),
                    br(),
                    br(),
                    actionButton(ns("act_save_as_edit_sr"), label = "Save as", class = "btn btn-primary"),
                    br(),
                    br(),
                    actionButton(ns("act_save_edit_sr"), label = "Save", class = "btn btn-primary")
             ),
             column(2,
                    wellPanel(
                      h4("Add an allele"),
                      numericInput(ns("add_al"), label = "Allele", value = NULL),
                      checkboxInput(ns("ck_auto_calc_sr"), label = "Calculate stutter ratios", value = FALSE, width = NULL),
                      disabled(numericInput(ns("alcor_bsr"), label = "Corrected allele (BSR)", value = NULL)),
                      disabled(numericInput(ns("alcor_fsr"), label = "Corrected allele (FSR)", value = NULL)),
                      disabled(numericInput(ns("alcor_dsr"), label = "Corrected allele (DSR)", value = NULL)),
                      disabled(numericInput(ns("alcor_m2sr"), label = "Corrected allele (M2SR)", value = NULL)),
                      actionButton(ns("act_add_al"), label = "Add")
                    ),
                    wellPanel(
                      h4("Delete an allele"),
                      uiOutput(ns("select_del_al")),
                      actionButton(ns("act_del_al"), label = "Delete")
                    )
             ),
             column(9,
                    textOutput(ns("fn_stutter_file")),
                    br(),
                    uiOutput(ns("select_mk")),
                    br(),
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

#' extrapolate_sr
#'
#' @description The function to extrapolate the stutter ratio of an allele
extrapolate_sr <- function(alcor, dt_result_est_sr, mk){
  if(!is.na(alcor)){
    index_mk <- which(dt_result_est_sr$Marker == mk)
    meth <- dt_result_est_sr$Model[index_mk]
    if(is.element(meth, c("Allele", "LUS", "Multi-seq"))){
      slope <- dt_result_est_sr$Slope[index_mk]
      intercept <- dt_result_est_sr$Intercept[index_mk]
      sr <- slope * alcor + intercept
      if(sr < 0){
        sr <- 0
      }
    }else if(is.element(meth, c("Uniform", "Uniform (Multiple loci)"))){
      sr <- dt_result_est_sr$Myu[index_mk]
      if(sr < 0){
        sr <- 0
      }
    }else{
      sr <- NA
    }
  }else{
    sr <- NA
  }
  return(sr)
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
      rv_edit_sr$fn_stutter_file <- NULL
      rv_edit_sr$dt_setting_est_sr <- NULL
      rv_edit_sr$dt_result_est_bsr <- NULL
      rv_edit_sr$dt_result_est_fsr <- NULL
      rv_edit_sr$dt_result_est_dsr <- NULL
      rv_edit_sr$dt_result_est_m2sr <- NULL
      rv_edit_sr$list_dt_average_sr <- NULL
      rv_edit_sr$list_trace_bsr <- NULL
      rv_edit_sr$list_trace_fsr <- NULL
      rv_edit_sr$list_trace_dsr <- NULL
      rv_edit_sr$list_trace_m2sr <- NULL
      rv_edit_sr$markers_auto <- character(0)
      rv_edit_sr$current_mk <- NULL

      shinyFileChoose(input, "load_stutter_file", roots = c(documents = paste0(path_pack, "/extdata/stutter_file")), session = session, filetypes = c("RData"))

      observeEvent(input$load_stutter_file, {
        req(input$load_stutter_file)
        path_stutter_file <- parseFilePaths(roots = c(documents = paste0(path_pack, "/extdata/stutter_file")), input$load_stutter_file)$datapath

        if(length(path_stutter_file) > 0){
          rv_edit_sr$fn_stutter_file <- fn_stutter_file <- basename(path_stutter_file)
          output$fn_stutter_file <- renderText({paste0("File name: ", fn_stutter_file)})

          rv_edit_sr$dt_setting_est_sr <- NULL
          rv_edit_sr$dt_result_est_bsr <- NULL
          rv_edit_sr$dt_result_est_fsr <- NULL
          rv_edit_sr$dt_result_est_dsr <- NULL
          rv_edit_sr$dt_result_est_m2sr <- NULL
          rv_edit_sr$list_dt_average_sr <- NULL
          rv_edit_sr$list_trace_bsr <- NULL
          rv_edit_sr$list_trace_fsr <- NULL
          rv_edit_sr$list_trace_dsr <- NULL
          rv_edit_sr$list_trace_m2sr <- NULL
          rv_edit_sr$markers_auto <- character(0)
          rv_edit_sr$current_mk <- NULL

          load(path_stutter_file)
          rv_edit_sr$dt_setting_est_sr <- dt_setting_est_sr
          rv_edit_sr$dt_result_est_bsr <- dt_result_est_bsr
          rv_edit_sr$dt_result_est_fsr <- dt_result_est_fsr
          rv_edit_sr$dt_result_est_dsr <- dt_result_est_dsr
          rv_edit_sr$dt_result_est_m2sr <- dt_result_est_m2sr
          rv_edit_sr$list_dt_average_sr <- list_dt_average_sr
          rv_edit_sr$list_trace_bsr <- list_trace_bsr
          rv_edit_sr$list_trace_fsr <- list_trace_fsr
          rv_edit_sr$list_trace_dsr <- list_trace_dsr
          rv_edit_sr$list_trace_m2sr <- list_trace_m2sr
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
        }
      })

      output$select_mk <- renderUI({
        markers_auto <- rv_edit_sr$markers_auto
        if(length(markers_auto) != 0){
          selectInput(session$ns("select_mk"), label = "Marker", choices = markers_auto, selected = markers_auto[1], selectize = FALSE)
        }else{
          NULL
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

      observeEvent(input$act_add_al, {
        add_al <- input$add_al
        if(!is.na(add_al)){
          current_mk <- rv_edit_sr$current_mk
          dt_ave_sr <- rv_edit_sr$list_dt_average_sr[[current_mk]]

          alcor_bsr <- input$alcor_bsr
          alcor_fsr <- input$alcor_fsr
          alcor_dsr <- input$alcor_dsr
          alcor_m2sr <- input$alcor_m2sr
          if(input$ck_auto_calc_sr){
            bsr <- extrapolate_sr(alcor_bsr, rv_edit_sr$dt_result_est_bsr, current_mk)
            fsr <- extrapolate_sr(alcor_fsr, rv_edit_sr$dt_result_est_fsr, current_mk)
            dsr <- extrapolate_sr(alcor_dsr, rv_edit_sr$dt_result_est_dsr, current_mk)
            m2sr <- extrapolate_sr(alcor_m2sr, rv_edit_sr$dt_result_est_m2sr, current_mk)
          }else{
            bsr <- NA
            fsr <- NA
            dsr <- NA
            m2sr <- NA
          }

          dt_ave_sr <- rbind(dt_ave_sr, list(as.character(add_al), alcor_bsr, bsr, alcor_fsr, fsr, alcor_dsr, dsr, alcor_m2sr, m2sr))
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
          selectInput(session$ns("del_al"), label = "Allele", choices = current_alleles, selected = current_alleles[1], selectize = FALSE)
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

      ## Save edited data ##

      observeEvent(input$act_save_as_edit_sr, {
        showModal(modalDialog(
          title = "Save stutter data",
          textInput(session$ns("fn_edit_sr"), label = "Enter the file name"),
          footer = tagList(
            actionButton(session$ns("act_save_as_edit_sr_2"), "Save"),
            modalButton("Close")
          )
        ))
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
        list_trace_bsr <- rv_edit_sr$list_trace_bsr
        list_trace_fsr <- rv_edit_sr$list_trace_fsr
        list_trace_dsr <- rv_edit_sr$list_trace_dsr
        list_trace_m2sr <- rv_edit_sr$list_trace_m2sr

        list_dt_average_sr <- replace_ave_sr(list_dt_average_sr)

        save(list = c("dt_setting_est_sr",
                      "dt_result_est_bsr", "dt_result_est_fsr", "dt_result_est_dsr", "dt_result_est_m2sr",
                      "list_trace_bsr", "list_trace_fsr", "list_trace_dsr", "list_trace_m2sr",
                      "list_dt_average_sr"),
             file = paste0(path_pack, "/extdata/stutter_file/", input$fn_edit_sr, ".RData"))

        removeModal()
      })

      observeEvent(input$act_save_edit_sr, {
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
        list_trace_bsr <- rv_edit_sr$list_trace_bsr
        list_trace_fsr <- rv_edit_sr$list_trace_fsr
        list_trace_dsr <- rv_edit_sr$list_trace_dsr
        list_trace_m2sr <- rv_edit_sr$list_trace_m2sr

        list_dt_average_sr <- replace_ave_sr(list_dt_average_sr)

        save(list = c("dt_setting_est_sr",
                      "dt_result_est_bsr", "dt_result_est_fsr", "dt_result_est_dsr", "dt_result_est_m2sr",
                      "list_trace_bsr", "list_trace_fsr", "list_trace_dsr", "list_trace_m2sr",
                      "list_dt_average_sr"),
             file = paste0(path_pack, "/extdata/stutter_file/", rv_edit_sr$fn_stutter_file))
      })
    }
  )
}
