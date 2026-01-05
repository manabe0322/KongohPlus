#' view_sr_ui
#'
#' @description The function to create the ui module for viewing stutter ratios
view_sr_ui <- function(id){
  ns <- NS(id)

  tabPanel(
    "View the stutter ratios",
    tags$style("
             .edit-modal .modal-dialog {
               width: 90% !important;   /* Expand up to 90% of the screen width */
               max-width: 2000px;       /* Set maximum width */
             }
             .edit-modal .modal-body {
               max-height: 80vh;        /* Limit height to 80% of the viewport */
               overflow-y: auto;        /* Enable vertical scrolling if needed */
             }
           "),
    titlePanel("View the stutter ratios"),
    br(),
    fluidRow(
      column(3,
             wellPanel(
               h4("Stutter data"),
               fluidRow(
                 column(3, br(), shinyFilesButton(ns("load_stutter_data"), label = "Load", title = "Select a file", multiple = FALSE)),
                 column(9, uiOutput(ns("fn_stutter_data")))
               ),
               selectInput(ns("select_mk"), label = h4("Marker"), choices = character(0), selected = NULL, multiple = FALSE, selectize = FALSE),
               selectInput(ns("select_st_type"), label = h4("Stutter type"), choices = character(0), selected = NULL, multiple = FALSE, selectize = FALSE)
             ),
             wellPanel(
               h4("Warnings"),
               htmlOutput(ns("warning_message")),
               br()
             )
      ),
      column(4,
             wellPanel(
               h4("Average stutter ratio"),
               dataTableOutput(ns("dt_average_sr_display"))
             ),
             wellPanel(
               h4("View all MLE attempts"),
               uiOutput(ns("buttons_all_mle"))
             )
      ),
      column(5,
             wellPanel(
               h4("Experimental data"),
               uiOutput(ns("ui_plot_sr")),
               br(),
               dataTableOutput(ns("dt_sr_data_display"))
             )
      )
    )
  )
}

#' create_dt_all_mle_sr
#'
#' @description The function to create the data.table of all MLE attempts for the stutter ratios
#' @param sr_type The type of the stutter ratios: BSR, FSR, DSR, or M2SR
#' @param dt_setting_est_sr The data.table for the stutter settings
#' @param opt_select_mk_sr The list of the optimization results for the arbitrary type of the stutter ratios
#' @param min_myu The minimum myu
#' @param max_myu The maximum myu
#' @param min_slope The minimum slope
#' @param max_slope The maximum slope
#' @param min_intercept The minimum intercept
#' @param max_intercept The maximum intercept
#' @param min_sigma The minimum sigma
#' @param max_sigma The maximum sigma
#' @param select_mk The selected marker
create_dt_all_mle_sr <- function(sr_type, dt_setting_est_sr, opt_select_mk_sr, min_myu, max_myu, min_slope, max_slope, min_intercept, max_intercept, min_sigma, max_sigma, select_mk){

  dt_all_mle <- NULL
  count_mle <- NULL
  name_converged <- NULL

  if(!is.null(select_mk)){

    model_sr <- dt_setting_est_sr[[paste0(sr_type, " model")]][dt_setting_est_sr[["Marker"]] == select_mk]

    if(!is.null(opt_select_mk_sr)){

      if(is.element(model_sr, c("Linear (Allele)", "Linear (Corrected Allele)"))){
        name_param <- c("Slope", "Intercept", "Sigma")
      }else{
        name_param <- c("Myu", "Sigma")
      }

      n_param <- length(name_param)
      count_mle <- length(opt_select_mk_sr)

      mat_all_mle <- matrix(NA, n_param + 1, count_mle)
      all_scores <- rep(NA, count_mle)

      for(i in 1:count_mle){
        opt_one <- opt_select_mk_sr[[i]]
        all_scores[i] <- score <- opt_one$value
        est_params <- opt_one$par

        vec_params <- rep(NA, n_param)
        if(n_param == 3){
          vec_params[1] <- est_params[1] * (max_slope - min_slope) + min_slope
          vec_params[2] <- est_params[2] * (max_intercept - min_intercept) + min_intercept
          vec_params[3] <- est_params[3] * (max_sigma - min_sigma) + min_sigma
        }else{
          vec_params[1] <- est_params[1] * (max_myu - min_myu) + min_myu
          vec_params[2] <- est_params[2] * (max_sigma - min_sigma) + min_sigma
        }

        mat_all_mle[1, i] <- - score / log(10)
        mat_all_mle[2:(n_param + 1), i] <- vec_params
      }

      dt_all_mle <- as.data.table(mat_all_mle)
      names(dt_all_mle) <- paste0("MLE ", 1:count_mle)

      dt_all_mle[, Parameter := c("Log10 (likelihood)", name_param)]
      setcolorder(dt_all_mle, c("Parameter", setdiff(names(dt_all_mle), "Parameter")))

      all_scores_round <- round(all_scores, 3)
      index_converged <- which(all_scores_round == min(all_scores_round))
      name_converged <- paste0("MLE ", index_converged)
    }
  }

  return(list(dt_all_mle = dt_all_mle, count_mle = count_mle, name_converged = name_converged))
}

#' view_sr_server
#'
#' @description The function to create the server module for viewing stutter ratios
view_sr_server <- function(id, path_pack){
  moduleServer(
    id,
    function(input, output, session){
      rv_view_sr <- reactiveValues()
      rv_view_sr$warning_message <- NULL
      rv_view_sr$fn_stutter_data <- NULL
      rv_view_sr$markers_auto <- character(0)
      rv_view_sr$list_dt_average_sr <- NULL
      rv_view_sr$dt_setting_est_sr <- NULL
      rv_view_sr$dt_sr_data <- NULL
      rv_view_sr$name_kit_stutter <- NULL

      rv_view_sr$opt_list_all_mk_bsr <- NULL
      rv_view_sr$opt_list_all_mk_fsr <- NULL
      rv_view_sr$opt_list_all_mk_dsr <- NULL
      rv_view_sr$opt_list_all_mk_m2sr <- NULL
      rv_view_sr$min_myu <- NULL
      rv_view_sr$max_myu <- NULL
      rv_view_sr$min_slope <- NULL
      rv_view_sr$max_slope <- NULL
      rv_view_sr$min_intercept <- NULL
      rv_view_sr$max_intercept <- NULL
      rv_view_sr$min_sigma <- NULL
      rv_view_sr$max_sigma <- NULL

      rv_view_sr$dt_average_sr <- NULL
      rv_view_sr$dt_st_type <- NULL
      rv_view_sr$meth_all_mk <- NULL
      rv_view_sr$meth <- NULL
      rv_view_sr$dt_st_mk <- NULL
      rv_view_sr$dt_kit_mk <- NULL
      rv_view_sr$df_select_sr <- NULL
      rv_view_sr$dt_sr_data_display <- NULL
      rv_view_sr$opt_select_mk_bsr <- NULL
      rv_view_sr$opt_select_mk_fsr <- NULL
      rv_view_sr$opt_select_mk_dsr <- NULL
      rv_view_sr$opt_select_mk_m2sr <- NULL

      shinyFileChoose(input, "load_stutter_data", roots = c(documents = paste0(path_pack, "/extdata/stutter")), session = session, filetypes = c("RData"))

      observeEvent(input$load_stutter_data, {
        req(input$load_stutter_data)
        path_stutter_data <- parseFilePaths(roots = c(documents = paste0(path_pack, "/extdata/stutter")), input$load_stutter_data)$datapath

        if(length(path_stutter_data) > 0){
          rv_view_sr$fn_stutter_data <- basename(path_stutter_data)

          load(path_stutter_data)

          rv_view_sr$warning_message <- data_est_sr$warning_message
          rv_view_sr$markers_auto <- names(data_est_sr$list_dt_average_sr)
          rv_view_sr$list_dt_average_sr <- data_est_sr$list_dt_average_sr
          rv_view_sr$dt_setting_est_sr <- data_est_sr$dt_setting_est_sr
          rv_view_sr$dt_sr_data <- data_est_sr$dt_sr_data
          rv_view_sr$name_kit_stutter <- data_est_sr$name_kit_stutter
          rv_view_sr$opt_list_all_mk_bsr <- data_est_sr$opt_list_all_mk_bsr
          rv_view_sr$opt_list_all_mk_fsr <- data_est_sr$opt_list_all_mk_fsr
          rv_view_sr$opt_list_all_mk_dsr <- data_est_sr$opt_list_all_mk_dsr
          rv_view_sr$opt_list_all_mk_m2sr <- data_est_sr$opt_list_all_mk_m2sr
          rv_view_sr$min_myu <- data_est_sr$min_myu
          rv_view_sr$max_myu <- data_est_sr$max_myu
          rv_view_sr$min_slope <- data_est_sr$min_slope
          rv_view_sr$max_slope <- data_est_sr$max_slope
          rv_view_sr$min_intercept <- data_est_sr$min_intercept
          rv_view_sr$max_intercept <- data_est_sr$max_intercept
          rv_view_sr$min_sigma <- data_est_sr$min_sigma
          rv_view_sr$max_sigma <- data_est_sr$max_sigma
        }
      })

      output$fn_stutter_data <- renderUI({
        disabled(textInput(session$ns("fn_stutter_data"), label = "", value = rv_view_sr$fn_stutter_data))
      })

      observeEvent(rv_view_sr$markers_auto, {
        markers_auto <- rv_view_sr$markers_auto
        if(length(markers_auto) == 0){
          updateSelectInput(session, "select_mk", choices = character(0), selected = character(0))
        }else{
          prev_mk <- isolate(input$select_mk)
          if(!is.null(prev_mk) && prev_mk %in% markers_auto){
            current_mk <- prev_mk
          }else{
            current_mk <- markers_auto[1]
          }
          updateSelectInput(session, "select_mk", choices = markers_auto, selected = current_mk)
        }
      })

      observeEvent(rv_view_sr$fn_stutter_data, {
        fn_stutter_data <- rv_view_sr$fn_stutter_data
        if(length(fn_stutter_data) == 0){
          updateSelectInput(session, "select_st_type", choices = character(0), selected = character(0))
        }else{
          prev_st_type <- isolate(input$select_st_type)
          if(!is.null(prev_st_type) && prev_st_type %in% c("Back stutter", "Forward stutter", "Double-back stutter", "Minus 2-nt stutter")){
            current_st_type <- prev_st_type
          }else{
            current_st_type <- "Back stutter"
          }
          updateSelectInput(session, "select_st_type", choices = c("Back stutter", "Forward stutter", "Double-back stutter", "Minus 2-nt stutter"), selected = current_st_type)
        }
      })

      observeEvent(list(input$select_mk, input$select_st_type), {
        req(rv_view_sr$markers_auto)
        req(rv_view_sr$list_dt_average_sr)
        req(rv_view_sr$dt_setting_est_sr)
        req(rv_view_sr$dt_sr_data)
        req(rv_view_sr$name_kit_stutter)
        req(rv_view_sr$opt_list_all_mk_bsr)
        req(rv_view_sr$opt_list_all_mk_fsr)
        req(rv_view_sr$opt_list_all_mk_dsr)
        req(rv_view_sr$opt_list_all_mk_m2sr)

        select_mk <- input$select_mk
        select_st_type <- input$select_st_type

        markers_auto <- rv_view_sr$markers_auto
        list_dt_average_sr <- rv_view_sr$list_dt_average_sr
        dt_setting_est_sr <- rv_view_sr$dt_setting_est_sr
        dt_sr_data <- rv_view_sr$dt_sr_data
        name_kit_stutter <- rv_view_sr$name_kit_stutter
        opt_list_all_mk_bsr <- rv_view_sr$opt_list_all_mk_bsr
        opt_list_all_mk_fsr <- rv_view_sr$opt_list_all_mk_fsr
        opt_list_all_mk_dsr <- rv_view_sr$opt_list_all_mk_dsr
        opt_list_all_mk_m2sr <- rv_view_sr$opt_list_all_mk_m2sr

        rv_view_sr$dt_average_sr <- NULL
        rv_view_sr$dt_st_type <- NULL
        rv_view_sr$meth_all_mk <- NULL
        rv_view_sr$meth <- NULL
        rv_view_sr$dt_st_mk <- NULL
        rv_view_sr$dt_kit_mk <- NULL
        rv_view_sr$df_select_sr <- NULL
        rv_view_sr$dt_sr_data_display <- NULL
        rv_view_sr$opt_select_mk_bsr <- NULL
        rv_view_sr$opt_select_mk_fsr <- NULL
        rv_view_sr$opt_select_mk_dsr <- NULL
        rv_view_sr$opt_select_mk_m2sr <- NULL

        if(!is.null(select_mk) && !is.null(select_st_type)){
          rv_view_sr$dt_average_sr <- list_dt_average_sr[[select_mk]]

          if(select_st_type == "Back stutter"){
            dt_st_type <- dt_sr_data[dt_sr_data$Bool_BS]
            meth_all_mk <- dt_setting_est_sr[["BSR model"]]
            meth <- meth_all_mk[dt_setting_est_sr$Marker == select_mk]
          }else if(select_st_type == "Forward stutter"){
            dt_st_type <- dt_sr_data[dt_sr_data$Bool_FS]
            meth_all_mk <- dt_setting_est_sr[["FSR model"]]
            meth <- meth_all_mk[dt_setting_est_sr$Marker == select_mk]
          }else if(select_st_type == "Double-back stutter"){
            dt_st_type <- dt_sr_data[dt_sr_data$Bool_DS]
            meth_all_mk <- dt_setting_est_sr[["DSR model"]]
            meth <- meth_all_mk[dt_setting_est_sr$Marker == select_mk]
          }else if(select_st_type == "Minus 2-nt stutter"){
            dt_st_type <- dt_sr_data[dt_sr_data$Bool_M2S]
            meth_all_mk <- dt_setting_est_sr[["M2SR model"]]
            meth <- meth_all_mk[dt_setting_est_sr$Marker == select_mk]
          }else{
            dt_st_type <- dt_sr_data[rep(FALSE, nrow(dt_sr_data))]
            meth_all_mk <- NULL
            meth <- NULL
          }
          rv_view_sr$dt_st_type <- dt_st_type
          rv_view_sr$meth_all_mk <- meth_all_mk
          rv_view_sr$meth <- meth

          if(nrow(dt_st_type) != 0){
            dt_st_mk <- dt_st_type[Marker == select_mk]
            dt_st_mk <- dt_st_mk[Allele != "Multiple"]
            dt_st_mk <- dt_st_mk[!dt_st_mk$Bool_SharedAllele]
            dt_st_mk <- dt_st_mk[!dt_st_mk$Bool_AlleleDropout]
            dt_st_mk <- dt_st_mk[dt_st_mk$Bool_InfoAllele]
            rv_view_sr$dt_st_mk <- dt_st_mk
          }

          index_mk <- which(markers_auto == select_mk)
          rv_view_sr$opt_select_mk_bsr <- opt_list_all_mk_bsr[[index_mk]]
          rv_view_sr$opt_select_mk_fsr <- opt_list_all_mk_fsr[[index_mk]]
          rv_view_sr$opt_select_mk_dsr <- opt_list_all_mk_dsr[[index_mk]]
          rv_view_sr$opt_select_mk_m2sr <- opt_list_all_mk_m2sr[[index_mk]]

          dt_kit <- fread(paste0(path_pack, "/extdata/kit/", name_kit_stutter, ".csv"))
          dt_kit <- try(change_data_type_kit(dt_kit))
          if(inherits(dt_kit, "try-error")){
            showModal(modalDialog(title = "Error", "The file format of the selected kit is invalid.", easyClose = TRUE, footer = NULL))
          }else{
            rv_view_sr$dt_kit_mk <- dt_kit[Marker == select_mk]
          }
        }
      })

      output$warning_message <- renderUI({
        req(rv_view_sr$warning_message)
        HTML(rv_view_sr$warning_message)
      })

      output$dt_average_sr_display <- renderDataTable({
        req(rv_view_sr$dt_average_sr)
        dt_average_sr <- rv_view_sr$dt_average_sr

        dt_average_sr_display <- dt_average_sr[, !c("Corrected allele (BSR)", "Corrected allele (FSR)", "Corrected allele (DSR)", "Corrected allele (M2SR)"), with = FALSE]

        datatable(
          dt_average_sr_display,
          selection = "none",
          options = list(dom = "t", iDisplayLength = nrow(dt_average_sr_display), ordering = FALSE, autoWidth = FALSE, scrollY = "500px", fixedHeader = TRUE),
          rownames = FALSE
        ) %>%
          formatSignif(columns = c("Average BSR", "Average FSR", "Average DSR", "Average M2SR"), digits = 4) %>%
          formatStyle(names(dt_average_sr_display), lineHeight = "80%")
      })

      output$buttons_all_mle <- renderUI({
        req(rv_view_sr$opt_list_all_mk_bsr)
        req(rv_view_sr$opt_list_all_mk_fsr)
        req(rv_view_sr$opt_list_all_mk_dsr)
        req(rv_view_sr$opt_list_all_mk_m2sr)

        tagList(
          actionButton(session$ns("act_all_mle_bsr"), label = "BSR"),
          actionButton(session$ns("act_all_mle_fsr"), label = "FSR"),
          actionButton(session$ns("act_all_mle_dsr"), label = "DSR"),
          actionButton(session$ns("act_all_mle_m2sr"), label = "M2SR")
        )
      })

      observeEvent(input$act_all_mle_bsr, {
        showModal(modalDialog(
          title = "All MLE attempts (BSR)",
          dataTableOutput(session$ns("dt_all_mle_bsr")),
          easyClose = TRUE,
          footer = NULL
        ) %>% tagAppendAttributes(class = "edit-modal")
        )
      })

      output$dt_all_mle_bsr <- renderDataTable({
        req(rv_view_sr$dt_setting_est_sr)
        req(rv_view_sr$opt_select_mk_bsr)
        req(rv_view_sr$min_myu)
        req(rv_view_sr$max_myu)
        req(rv_view_sr$min_slope)
        req(rv_view_sr$max_slope)
        req(rv_view_sr$min_intercept)
        req(rv_view_sr$max_intercept)
        req(rv_view_sr$min_sigma)
        req(rv_view_sr$max_sigma)
        req(input$select_mk)

        dt_setting_est_sr <- rv_view_sr$dt_setting_est_sr
        opt_select_mk_sr <- rv_view_sr$opt_select_mk_bsr
        min_myu <- rv_view_sr$min_myu
        max_myu <- rv_view_sr$max_myu
        min_slope <- rv_view_sr$min_slope
        max_slope <- rv_view_sr$max_slope
        min_intercept <- rv_view_sr$min_intercept
        max_intercept <- rv_view_sr$max_intercept
        min_sigma <- rv_view_sr$min_sigma
        max_sigma <- rv_view_sr$max_sigma

        select_mk <- input$select_mk

        tmp <- create_dt_all_mle_sr("BSR", dt_setting_est_sr, opt_select_mk_sr, min_myu, max_myu, min_slope, max_slope, min_intercept, max_intercept, min_sigma, max_sigma, select_mk)
        dt_all_mle <- tmp$dt_all_mle
        count_mle <- tmp$count_mle
        name_converged <- tmp$name_converged

        if(all(!is.null(dt_all_mle), !is.null(count_mle), !is.null(name_converged))){
          datatable(
            dt_all_mle,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_all_mle), ordering = FALSE, autowidth = TRUE),
            rownames = FALSE
          ) %>%
            formatSignif(columns = paste0("MLE ", 1:count_mle), digits = 6) %>%
            formatStyle(names(dt_all_mle), lineHeight = "80%") %>%
            formatStyle(name_converged, backgroundColor = "lightblue")
        }
      })

      observeEvent(input$act_all_mle_fsr, {
        showModal(modalDialog(
          title = "All MLE attempts (FSR)",
          dataTableOutput(session$ns("dt_all_mle_fsr")),
          easyClose = TRUE,
          footer = NULL
        ) %>% tagAppendAttributes(class = "edit-modal")
        )
      })

      output$dt_all_mle_fsr <- renderDataTable({
        req(rv_view_sr$dt_setting_est_sr)
        req(rv_view_sr$opt_select_mk_fsr)
        req(rv_view_sr$min_myu)
        req(rv_view_sr$max_myu)
        req(rv_view_sr$min_slope)
        req(rv_view_sr$max_slope)
        req(rv_view_sr$min_intercept)
        req(rv_view_sr$max_intercept)
        req(rv_view_sr$min_sigma)
        req(rv_view_sr$max_sigma)
        req(input$select_mk)

        dt_setting_est_sr <- rv_view_sr$dt_setting_est_sr
        opt_select_mk_sr <- rv_view_sr$opt_select_mk_fsr
        min_myu <- rv_view_sr$min_myu
        max_myu <- rv_view_sr$max_myu
        min_slope <- rv_view_sr$min_slope
        max_slope <- rv_view_sr$max_slope
        min_intercept <- rv_view_sr$min_intercept
        max_intercept <- rv_view_sr$max_intercept
        min_sigma <- rv_view_sr$min_sigma
        max_sigma <- rv_view_sr$max_sigma

        select_mk <- input$select_mk

        tmp <- create_dt_all_mle_sr("FSR", dt_setting_est_sr, opt_select_mk_sr, min_myu, max_myu, min_slope, max_slope, min_intercept, max_intercept, min_sigma, max_sigma, select_mk)
        dt_all_mle <- tmp$dt_all_mle
        count_mle <- tmp$count_mle
        name_converged <- tmp$name_converged

        if(all(!is.null(dt_all_mle), !is.null(count_mle), !is.null(name_converged))){
          datatable(
            dt_all_mle,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_all_mle), ordering = FALSE, autowidth = TRUE),
            rownames = FALSE
          ) %>%
            formatSignif(columns = paste0("MLE ", 1:count_mle), digits = 6) %>%
            formatStyle(names(dt_all_mle), lineHeight = "80%") %>%
            formatStyle(name_converged, backgroundColor = "lightblue")
        }
      })

      observeEvent(input$act_all_mle_dsr, {
        showModal(modalDialog(
          title = "All MLE attempts (DSR)",
          dataTableOutput(session$ns("dt_all_mle_dsr")),
          easyClose = TRUE,
          footer = NULL
        ) %>% tagAppendAttributes(class = "edit-modal")
        )
      })

      output$dt_all_mle_dsr <- renderDataTable({
        req(rv_view_sr$dt_setting_est_sr)
        req(rv_view_sr$opt_select_mk_dsr)
        req(rv_view_sr$min_myu)
        req(rv_view_sr$max_myu)
        req(rv_view_sr$min_slope)
        req(rv_view_sr$max_slope)
        req(rv_view_sr$min_intercept)
        req(rv_view_sr$max_intercept)
        req(rv_view_sr$min_sigma)
        req(rv_view_sr$max_sigma)
        req(input$select_mk)

        dt_setting_est_sr <- rv_view_sr$dt_setting_est_sr
        opt_select_mk_sr <- rv_view_sr$opt_select_mk_dsr
        min_myu <- rv_view_sr$min_myu
        max_myu <- rv_view_sr$max_myu
        min_slope <- rv_view_sr$min_slope
        max_slope <- rv_view_sr$max_slope
        min_intercept <- rv_view_sr$min_intercept
        max_intercept <- rv_view_sr$max_intercept
        min_sigma <- rv_view_sr$min_sigma
        max_sigma <- rv_view_sr$max_sigma

        select_mk <- input$select_mk

        tmp <- create_dt_all_mle_sr("DSR", dt_setting_est_sr, opt_select_mk_sr, min_myu, max_myu, min_slope, max_slope, min_intercept, max_intercept, min_sigma, max_sigma, select_mk)
        dt_all_mle <- tmp$dt_all_mle
        count_mle <- tmp$count_mle
        name_converged <- tmp$name_converged

        if(all(!is.null(dt_all_mle), !is.null(count_mle), !is.null(name_converged))){
          datatable(
            dt_all_mle,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_all_mle), ordering = FALSE, autowidth = TRUE),
            rownames = FALSE
          ) %>%
            formatSignif(columns = paste0("MLE ", 1:count_mle), digits = 6) %>%
            formatStyle(names(dt_all_mle), lineHeight = "80%") %>%
            formatStyle(name_converged, backgroundColor = "lightblue")
        }
      })

      observeEvent(input$act_all_mle_m2sr, {
        showModal(modalDialog(
          title = "All MLE attempts (M2SR)",
          dataTableOutput(session$ns("dt_all_mle_m2sr")),
          easyClose = TRUE,
          footer = NULL
        ) %>% tagAppendAttributes(class = "edit-modal")
        )
      })

      output$dt_all_mle_m2sr <- renderDataTable({
        req(rv_view_sr$dt_setting_est_sr)
        req(rv_view_sr$opt_select_mk_m2sr)
        req(rv_view_sr$min_myu)
        req(rv_view_sr$max_myu)
        req(rv_view_sr$min_slope)
        req(rv_view_sr$max_slope)
        req(rv_view_sr$min_intercept)
        req(rv_view_sr$max_intercept)
        req(rv_view_sr$min_sigma)
        req(rv_view_sr$max_sigma)
        req(input$select_mk)

        dt_setting_est_sr <- rv_view_sr$dt_setting_est_sr
        opt_select_mk_sr <- rv_view_sr$opt_select_mk_m2sr
        min_myu <- rv_view_sr$min_myu
        max_myu <- rv_view_sr$max_myu
        min_slope <- rv_view_sr$min_slope
        max_slope <- rv_view_sr$max_slope
        min_intercept <- rv_view_sr$min_intercept
        max_intercept <- rv_view_sr$max_intercept
        min_sigma <- rv_view_sr$min_sigma
        max_sigma <- rv_view_sr$max_sigma
        select_mk <- input$select_mk

        tmp <- create_dt_all_mle_sr("M2SR", dt_setting_est_sr, opt_select_mk_sr, min_myu, max_myu, min_slope, max_slope, min_intercept, max_intercept, min_sigma, max_sigma, select_mk)
        dt_all_mle <- tmp$dt_all_mle
        count_mle <- tmp$count_mle
        name_converged <- tmp$name_converged

        if(all(!is.null(dt_all_mle), !is.null(count_mle), !is.null(name_converged))){
          datatable(
            dt_all_mle,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_all_mle), ordering = FALSE, autowidth = TRUE),
            rownames = FALSE
          ) %>%
            formatSignif(columns = paste0("MLE ", 1:count_mle), digits = 6) %>%
            formatStyle(names(dt_all_mle), lineHeight = "80%") %>%
            formatStyle(name_converged, backgroundColor = "lightblue")
        }else{
          return(NULL)
        }
      })

      output$ui_plot_sr <- renderUI({
        req(rv_view_sr$dt_st_type)
        req(rv_view_sr$meth_all_mk)
        req(rv_view_sr$meth)
        req(rv_view_sr$dt_st_mk)
        req(rv_view_sr$dt_kit_mk)
        req(rv_view_sr$markers_auto)

        dt_st_type <- rv_view_sr$dt_st_type
        meth_all_mk <- rv_view_sr$meth_all_mk
        meth <- rv_view_sr$meth
        dt_st_mk <- rv_view_sr$dt_st_mk
        dt_kit_mk <- rv_view_sr$dt_kit_mk
        markers_auto <- rv_view_sr$markers_auto

        plotOutput(session$ns("plot_sr"), width = "100%", height = "300px")
      })

      output$plot_sr <- renderPlot({
        req(rv_view_sr$dt_st_type)
        req(rv_view_sr$meth_all_mk)
        req(rv_view_sr$meth)
        req(rv_view_sr$dt_st_mk)
        req(rv_view_sr$dt_kit_mk)
        req(rv_view_sr$markers_auto)

        dt_st_type <- rv_view_sr$dt_st_type
        meth_all_mk <- rv_view_sr$meth_all_mk
        meth <- rv_view_sr$meth
        dt_st_mk <- rv_view_sr$dt_st_mk
        dt_kit_mk <- rv_view_sr$dt_kit_mk
        markers_auto <- rv_view_sr$markers_auto

        select_mk <- input$select_mk
        select_st_type <- input$select_st_type
        df_select_sr <- rv_view_sr$df_select_sr

        x_al <- numeric(0)
        y_sr <- numeric(0)

        if(nrow(dt_st_mk) != 0){
          # Extract (corrected) allele repeat numbers
          if(meth == "Linear (Corrected Allele)"){
            alcor <- dt_kit_mk[["Corrected Allele"]]
            x_al <- alcor[match(dt_st_mk$Allele, dt_kit_mk[["Allele"]])]
          }else{
            x_al <- as.numeric(dt_st_mk$Allele)
          }

          y_sr <- dt_st_mk$SR
        }

        x_al_multi_other_mk <- numeric(0)
        y_sr_multi_other_mk <- numeric(0)

        if(meth == "Uniform (Multiple loci)"){

          for(i in 1:length(markers_auto)){
            if(markers_auto[i] != select_mk && meth_all_mk[i] == "Uniform (Multiple loci)"){
              dt_st_other_mk <- dt_st_type[Marker == markers_auto[i]]
              dt_st_other_mk <- dt_st_other_mk[Allele != "Multiple"]
              dt_st_other_mk <- dt_st_other_mk[!dt_st_other_mk$Bool_SharedAllele]
              dt_st_other_mk <- dt_st_other_mk[!dt_st_other_mk$Bool_AlleleDropout]
              dt_st_other_mk <- dt_st_other_mk[dt_st_other_mk$Bool_InfoAllele]

              x_al_multi_other_mk <- c(x_al_multi_other_mk, as.numeric(dt_st_other_mk$Allele))
              y_sr_multi_other_mk <- c(y_sr_multi_other_mk, dt_st_other_mk$SR)
            }
          }
        }

        # Label for x-axis
        if(meth == "Linear (Corrected Allele)"){
          x_lab <- "Corrected Allele"
        }else{
          x_lab <- "Allele"
        }

        # Label for y-axis
        if(select_st_type == "Back stutter"){
          ylab = "BSR"
        }else if(select_st_type == "Forward stutter"){
          ylab = "FSR"
        }else if(select_st_type == "Double-back stutter"){
          ylab = "DSR"
        }else if(select_st_type == "Minus 2-nt stutter"){
          ylab = "M2SR"
        }else{
          ylab = "SR"
        }

        df_target_mk <- data.frame(Allele = x_al, SR = y_sr)
        color_values <- c("Selected markers" = "#00afcc")

        if(length(x_al_multi_other_mk) == 0){
          base_plot <- ggplot() +
            geom_point(data = df_target_mk, aes(Allele, SR, color = "Selected markers"))

        }else{
          df_other_mk <- data.frame(Allele = x_al_multi_other_mk, SR = y_sr_multi_other_mk)

          base_plot <- ggplot() +
            geom_point(data = df_other_mk, aes(Allele, SR, color = "Other markers")) +
            geom_point(data = df_target_mk, aes(Allele, SR, color = "Selected markers"))

          color_values <- c(color_values, "Other markers" = "#c0c0c0")
        }

        if(!is.null(df_select_sr)){
          base_plot <- base_plot + geom_point(data = df_select_sr, aes(Allele, SR, color = "Selected data"))

          color_values <- c(color_values, "Selected data" = "#e9bc00")
        }

        base_plot +
          labs(x = x_lab, y = ylab, color = NULL) +
          scale_color_manual(values = color_values) +
          theme_minimal() +
          theme(
            plot.title = element_text(size = 16),
            legend.text = element_text(size = 14),
            axis.title.x = element_text(size = 14, margin = margin(t = 20)),
            axis.title.y = element_text(size = 14, margin = margin(r = 20)),
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14)
          )
      })

      output$dt_sr_data_display <- renderDataTable({
        req(rv_view_sr$dt_st_mk)

        dt_st_mk <- rv_view_sr$dt_st_mk

        dt_sr_data_display <- dt_st_mk[, c("SampleName", "Allele", "SR", "Peak height (stutter)", "Peak height (allele)"), with = FALSE]
        dt_sr_data_display[["SampleName"]] <- factor(dt_sr_data_display[["SampleName"]], levels = unique(dt_sr_data_display[["SampleName"]]), labels = unique(dt_sr_data_display[["SampleName"]]))
        tmp_alleles <- as.character(sort(as.numeric(unique(dt_sr_data_display[["Allele"]]))))
        dt_sr_data_display[["Allele"]] <- factor(dt_sr_data_display[["Allele"]], levels = tmp_alleles, labels = tmp_alleles)
        rv_view_sr$dt_sr_data_display <- dt_sr_data_display

        datatable(
          dt_sr_data_display,
          filter = "top",
          selection = list(mode = "single", target = "row"),
          options = list(
            iDisplayLength = 10,
            autoWidth = FALSE,
            scrollX = TRUE,
            columnDefs = list(
              list(
                targets = 0,
                width = '250px'
              ),
              list(
                targets = 1,
                width = '75px'
              ),
              list(
                targets = c(2, 3, 4),
                searchable = FALSE
              )
            )
          ),
          rownames = FALSE
        ) %>%
          formatSignif(columns = c("SR"), digits = 4) %>%
          formatStyle(names(dt_sr_data_display), lineHeight = "80%")
      })

      observeEvent(input$dt_sr_data_display_rows_selected, {
        index_select <- input$dt_sr_data_display_rows_selected
        dt_kit_mk <- rv_view_sr$dt_kit_mk

        if(length(index_select) == 1 && !is.null(dt_kit_mk)){
          dt_sr_data_display <- rv_view_sr$dt_sr_data_display
          al <- as.character(dt_sr_data_display[["Allele"]][index_select])
          sr <- dt_sr_data_display[["SR"]][index_select]

          meth <- rv_view_sr$meth
          if(meth == "Linear (Corrected Allele)"){

            alcor <- dt_kit_mk[["Corrected Allele"]][dt_kit_mk[["Allele"]] == al]
          }else{
            alcor <- as.numeric(al)
          }

          rv_view_sr$df_select_sr <- data.frame(Allele = alcor, SR = sr)
        }else{
          rv_view_sr$df_select_sr <- NULL
        }
      })
    }
  )
}
