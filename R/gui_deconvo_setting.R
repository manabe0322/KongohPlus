#' deconvo_setting_ui
#'
#' @description The function to create the ui module for setting conditions of deconvolution
deconvo_setting_ui <- function(id){
  ns <- NS(id)

  tabPanel("Setting",
           useWaiter(),
           tags$style(HTML("
                         .handsontable .htCore td {
                           vertical-align: middle !important;
                           text-align: center !important;
                         }"
           )),
           br(),
           fluidRow(
             column(4,
                    wellPanel(
                      h4("Number of contributors"),
                      uiOutput(ns("ui_noc_from")),
                      uiOutput(ns("ui_noc_to")),
                      br(),
                      h4("Theta"),
                      uiOutput(ns("ui_theta")),
                      br(),
                      h4("Analysis method"),
                      fluidRow(
                        column(3, br(), shinyFilesButton(ns("load_anameth_deconvo"), label = "Load", title = "Select a file", multiple = FALSE)),
                        column(9, uiOutput(ns("anameth_deconvo")))
                      ),
                      br(),
                      actionButton(ns("act_advance_deconvo"), label = "Advanced setting"),
                      br(),
                      br(),
                      actionButton(ns("act_analyze_deconvo"), label = "Analyze", class = "btn btn-primary btn-lg")
                    )
             ),
             column(4,
                    h4("Crime stain profile"),
                    dataTableOutput(ns("dt_info_csp")),
                    br(),
                    h4("Allele frequency"),
                    dataTableOutput(ns("dt_info_af")),
                    br()
             ),
             column(4,
                    h4(textOutput(ns("title_ignore_mk"))),
                    rHandsontableOutput(ns("dt_ignore_mk"))
             )
           )
  )
}

#' deconvo_setting_server
#'
#' @description The function to create the server module for setting conditions of deconvolution
deconvo_setting_server <- function(id, path_pack, rv_input){
  moduleServer(
    id,
    function(input, output, session){
      rv_deconvo_setting <- reactiveValues()
      rv_deconvo_setting$dt_csp <- NULL
      rv_deconvo_setting$dt_af <- NULL
      rv_deconvo_setting$fn_csp <- NULL
      rv_deconvo_setting$fn_af <- NULL
      rv_deconvo_setting$hnc_from <- 1
      rv_deconvo_setting$hnc_to <- 1
      rv_deconvo_setting$theta <- 0
      rv_deconvo_setting$fn_anameth <- ""
      rv_deconvo_setting$dt_anameth <- NULL
      rv_deconvo_setting$min_temp <- 1
      rv_deconvo_setting$max_temp <- 32000
      rv_deconvo_setting$min_rep_ae <- 0.2
      rv_deconvo_setting$max_rep_ae <- 5
      rv_deconvo_setting$min_deg <- 0.0001
      rv_deconvo_setting$stop_iter <- 500
      rv_deconvo_setting$max_iter <- 5000
      rv_deconvo_setting$num_converged <- 3
      rv_deconvo_setting$max_mle_runs <- 10
      rv_deconvo_setting$bool_correct_noc <- TRUE
      rv_deconvo_setting$bool_correct_theta <- TRUE
      rv_deconvo_setting$bool_correct_advanced_temp <- TRUE
      rv_deconvo_setting$bool_correct_advanced_rep_ae <- TRUE
      rv_deconvo_setting$bool_correct_advanced_deg <- TRUE
      rv_deconvo_setting$bool_correct_advanced_1 <- TRUE
      rv_deconvo_setting$bool_correct_advanced_2 <- TRUE
      rv_deconvo_setting$dt_ignore_mk <- NULL

      rv_deconvo <- reactiveValues()
      rv_deconvo$data_deconvo <- NULL

      waiter_ui <- Waiter$new(html = tagList(spin_3k(),
                                             h3("Checking input files...")),
                              color = transparent(.5))

      observe({
        req(rv_input)
        req(rv_deconvo)
        data_deconvo <- rv_deconvo$data_deconvo

        if(is.null(data_deconvo)){
          rv_deconvo_setting$dt_csp <- rv_input$dt_csp
          rv_deconvo_setting$dt_af <- rv_input$dt_af
          rv_deconvo_setting$fn_csp <- rv_input$fn_csp
          rv_deconvo_setting$fn_af <- rv_input$fn_af
        }else{
          rv_deconvo_setting$dt_csp <- data_deconvo$dt_csp
          rv_deconvo_setting$dt_af <- data_deconvo$dt_af
          rv_deconvo_setting$fn_csp <- data_deconvo$fn_csp
          rv_deconvo_setting$fn_af <- data_deconvo$fn_af
          hnc_all <- data_deconvo$hnc_all
          rv_deconvo_setting$hnc_to <- max(hnc_all)
          rv_deconvo_setting$hnc_from <- min(hnc_all)
          rv_deconvo_setting$theta <- data_deconvo$theta
          rv_deconvo_setting$fn_anameth <- data_deconvo$fn_anameth
          rv_deconvo_setting$dt_ignore_mk <- data_deconvo$dt_ignore_mk
        }
      })

      output$dt_info_csp <- renderDataTable(server = FALSE, {
        dt_csp <- rv_deconvo_setting$dt_csp
        dt_info_csp <- data.table("Sample Name" = unique(dt_csp$SampleName))

        datatable(
          dt_info_csp,
          selection = "none",
          options = list(dom = "t", iDisplayLength = nrow(dt_info_csp), ordering = FALSE, autowidth = TRUE),
          rownames = FALSE,
          escape = FALSE
        )
      })

      output$dt_info_af <- renderDataTable(server = FALSE, {
        dt_info_af <- data.table("File Name" = rv_deconvo_setting$fn_af)

        datatable(
          dt_info_af,
          selection = "none",
          options = list(dom = "t", iDisplayLength = nrow(dt_info_af), ordering = FALSE, autowidth = TRUE),
          rownames = FALSE,
          escape = FALSE
        )
      })

      # --- Set the number of contributors --- #

      output$ui_noc_from <- renderUI({
        numericInput(session$ns("noc_from"), label = "From", value = rv_deconvo_setting$hnc_from, min = 1, max = 4, step = 1)
      })

      output$ui_noc_to <- renderUI({
        numericInput(session$ns("noc_to"), label = "To", value = rv_deconvo_setting$hnc_to, min = 1, max = 4, step = 1)
      })

      observeEvent(input$noc_from | input$noc_to, {
        noc_from <- input$noc_from
        noc_to <- input$noc_to

        hideFeedback("noc_from")
        hideFeedback("noc_to")
        bool_correct_noc <- TRUE

        if(!is.integer(noc_from) || noc_from < 1 || noc_from > 4){
          showFeedbackDanger(inputId = "noc_from", text = "Only 1, 2, 3, and 4 is allowed.")
          bool_correct_noc <- FALSE
        }

        if(!is.integer(noc_to) || noc_to < 1 || noc_to > 4){
          showFeedbackDanger(inputId = "noc_to", text = "Only 1, 2, 3, and 4 is allowed.")
          bool_correct_noc <- FALSE
        }

        if(is.integer(noc_from) && is.integer(noc_to) && noc_from > noc_to){
          showFeedbackDanger(inputId = "noc_from", text = "The From number must not exceed the To number.")
          showFeedbackDanger(inputId = "noc_to", text = "The From number must not exceed the To number.")
          bool_correct_noc <- FALSE
        }

        rv_deconvo_setting$bool_correct_noc <- bool_correct_noc

        if(all(bool_correct_noc, rv_deconvo_setting$bool_correct_theta, rv_deconvo_setting$fn_anameth != "")){
          enable("act_analyze_deconvo")
        }else{
          disable("act_analyze_deconvo")
        }
      })

      # --- Set the theta value --- #

      output$ui_theta <- renderUI({
        numericInput(session$ns("theta"), label = "", value = rv_deconvo_setting$theta, min = 0, max = 1, step = 0.01)
      })

      observeEvent(input$theta, {
        theta <- input$theta

        if(!is.numeric(theta) || theta < 0 || theta > 1){
          showFeedbackDanger(inputId = "theta", text = "The value must be between 0 and 1.")
          bool_correct_theta <- FALSE
        }else{
          hideFeedback("theta")
          bool_correct_theta <- TRUE
        }

        if(all(bool_correct_theta, rv_deconvo_setting$bool_correct_noc, rv_deconvo_setting$fn_anameth != "")){
          enable("act_analyze_deconvo")
        }else{
          disable("act_analyze_deconvo")
        }

        rv_deconvo_setting$bool_correct_theta <- bool_correct_theta
      })

      # --- Set analysis method --- #

      shinyFileChoose(input, "load_anameth_deconvo", roots = c(documents = paste0(path_pack, "/extdata/anameth")), session = session, filetypes = c("csv"))

      observeEvent(input$load_anameth_deconvo, {
        req(input$load_anameth_deconvo)
        path_anameth <- parseFilePaths(roots = c(documents = paste0(path_pack, "/extdata/anameth")), input$load_anameth_deconvo)$datapath
        if(length(path_anameth) > 0){
          rv_deconvo_setting$fn_anameth <- basename(path_anameth)
          if(all(rv_deconvo_setting$bool_correct_noc, rv_deconvo_setting$bool_correct_theta, rv_deconvo_setting$fn_anameth != "")){
            enable("act_analyze_deconvo")
          }else{
            disable("act_analyze_deconvo")
          }
        }
      })

      output$anameth_deconvo <- renderUI({
        disabled(textInput(session$ns("anameth_deconvo"), label = "", value = rv_deconvo_setting$fn_anameth))
      })

      observeEvent(input$anameth_deconvo, {
        fn_anameth <- input$anameth_deconvo
        if(fn_anameth != ""){
          rv_deconvo_setting$dt_anameth <- dt_anameth <- fread(paste0(path_pack, "/extdata/anameth/", fn_anameth))
          markers_auto <- dt_anameth$Subcategory[dt_anameth$Category == "AT"]

          output$title_ignore_mk <- renderText({"Ignored marker"})

          rv_deconvo_setting$dt_ignore_mk <- data.table(Ignore = rep(FALSE, length(markers_auto)),
                                                        Marker = markers_auto)
        }
      })

      output$dt_ignore_mk <- renderRHandsontable({
        dt_ignore_mk <- rv_deconvo_setting$dt_ignore_mk

        if(!is.null(dt_ignore_mk)){
          rhandsontable(dt_ignore_mk, rowHeaders = FALSE) %>%
            hot_table(fillHandle = FALSE, rowHeights = 30) %>%
            hot_col("Ignore", type = "checkbox", width = 100) %>%
            hot_col("Marker", type = "text", readOnly = TRUE, width = 100)
        }
      })

      # --- Advanced setting --- #

      observeEvent(input$act_advance_deconvo, {
        showModal(modalDialog(
          title = "MLE conditions",
          fluidRow(
            column(6,
                   numericInput(session$ns("stop_iter"), label = "Minimum number of iterations", value = rv_deconvo_setting$stop_iter),
                   numericInput(session$ns("num_converged"), label = "Minimum number of MLE attempts", value = rv_deconvo_setting$num_converged),
                   numericInput(session$ns("min_temp"), label = "Minimum template amount", value = rv_deconvo_setting$min_temp),
                   numericInput(session$ns("min_rep_ae"), label = "Minimum replicate amplification efficiency", value = rv_deconvo_setting$min_rep_ae),
                   numericInput(session$ns("min_deg"), label = "Minimum degradation parameter", value = rv_deconvo_setting$min_deg)
            ),
            column(6,
                   numericInput(session$ns("max_iter"), label = "Maximum number of iterations", value = rv_deconvo_setting$max_iter),
                   numericInput(session$ns("max_mle_runs"), label = "Maximum number of MLE attempts", value = rv_deconvo_setting$max_mle_runs),
                   numericInput(session$ns("max_temp"), label = "Maximum template amount", value = rv_deconvo_setting$max_temp),
                   numericInput(session$ns("max_rep_ae"), label = "Maximum replicate amplification efficiency", value = rv_deconvo_setting$max_rep_ae)
            )
          ),
          footer = tagList(
            actionButton(session$ns("act_advance_deconvo_save"), "Save"),
            modalButton("Cancel")
          ),
          size = "l"
        ))
      })

      observeEvent(input$min_temp | input$max_temp, {
        min_temp <- input$min_temp
        max_temp <- input$max_temp

        hideFeedback("min_temp")
        hideFeedback("max_temp")
        bool_correct_advanced_temp <- TRUE

        if(!is.numeric(min_temp) || min_temp <= 0){
          showFeedbackDanger(inputId = "min_temp", text = "Value must be positive.")
          bool_correct_advanced_temp <- FALSE
        }

        if(!is.numeric(max_temp) || max_temp <= 0){
          showFeedbackDanger(inputId = "max_temp", text = "Value must be positive.")
          bool_correct_advanced_temp <- FALSE
        }

        if(is.numeric(min_temp) && is.numeric(max_temp) && max_temp < min_temp){
          showFeedbackDanger(inputId = "min_temp", text = "Minimum value cannot be greater than maximum value.")
          showFeedbackDanger(inputId = "max_temp", text = "Minimum value cannot be greater than maximum value.")
          bool_correct_advanced_temp <- FALSE
        }

        rv_deconvo_setting$bool_correct_advanced_temp <- bool_correct_advanced_temp
      })

      observeEvent(input$min_rep_ae | input$max_rep_ae, {
        min_rep_ae <- input$min_rep_ae
        max_rep_ae <- input$max_rep_ae

        hideFeedback("min_rep_ae")
        hideFeedback("max_rep_ae")
        bool_correct_advanced_rep_ae <- TRUE

        if(!is.numeric(min_rep_ae) || min_rep_ae <= 0){
          showFeedbackDanger(inputId = "min_rep_ae", text = "Value must be positive.")
          bool_correct_advanced_rep_ae <- FALSE
        }

        if(!is.numeric(max_rep_ae) || max_rep_ae <= 0){
          showFeedbackDanger(inputId = "max_rep_ae", text = "Value must be positive.")
          bool_correct_advanced_rep_ae <- FALSE
        }

        if(is.numeric(min_rep_ae) && is.numeric(max_rep_ae) && max_rep_ae < min_rep_ae){
          showFeedbackDanger(inputId = "min_rep_ae", text = "Minimum value cannot be greater than maximum value.")
          showFeedbackDanger(inputId = "max_rep_ae", text = "Minimum value cannot be greater than maximum value.")
          bool_correct_advanced_rep_ae <- FALSE
        }

        rv_deconvo_setting$bool_correct_advanced_rep_ae <- bool_correct_advanced_rep_ae
      })

      observeEvent(input$min_deg, {
        min_deg <- input$min_deg

        hideFeedback("min_deg")
        bool_correct_advanced_deg <- TRUE

        if(!is.numeric(min_deg) || min_deg <= 0 || min_deg > 1){
          showFeedbackDanger(inputId = "min_deg", text = "Value must be between 0 and 1.")
          bool_correct_advanced_deg <- FALSE
        }

        rv_deconvo_setting$bool_correct_advanced_deg <- bool_correct_advanced_deg
      })

      observeEvent(input$stop_iter | input$max_iter, {
        stop_iter <- input$stop_iter
        max_iter <- input$max_iter

        hideFeedback("stop_iter")
        hideFeedback("max_iter")
        bool_correct_advanced_1 <- TRUE

        if(!is.integer(stop_iter) || stop_iter <= 0){
          showFeedbackDanger(inputId = "stop_iter", text = "Value must be an integer greater than zero.")
          bool_correct_advanced_1 <- FALSE
        }

        if(!is.integer(max_iter) || max_iter <= 0){
          showFeedbackDanger(inputId = "max_iter", text = "Value must be an integer greater than zero.")
          bool_correct_advanced_1 <- FALSE
        }

        if(is.integer(stop_iter) && is.integer(max_iter) && stop_iter > max_iter){
          showFeedbackDanger(inputId = "stop_iter", text = "Minimum value cannot be greater than maximum value.")
          showFeedbackDanger(inputId = "max_iter", text = "Minimum value cannot be greater than maximum value.")
          bool_correct_advanced_1 <- FALSE
        }

        rv_deconvo_setting$bool_correct_advanced_1 <- bool_correct_advanced_1
      })

      observeEvent(input$num_converged | input$max_mle_runs, {
        num_converged <- input$num_converged
        max_mle_runs <- input$max_mle_runs

        hideFeedback("num_converged")
        hideFeedback("max_mle_runs")
        bool_correct_advanced_2 <- TRUE

        if(!is.integer(num_converged) || num_converged <= 0){
          showFeedbackDanger(inputId = "num_converged", text = "Value must be an integer greater than zero.")
          bool_correct_advanced_2 <- FALSE
        }

        if(!is.integer(max_mle_runs) || max_mle_runs <= 0){
          showFeedbackDanger(inputId = "max_mle_runs", text = "Value must be an integer greater than zero.")
          bool_correct_advanced_2 <- FALSE
        }

        if(is.integer(num_converged) && is.integer(max_mle_runs) && num_converged > max_mle_runs){
          showFeedbackDanger(inputId = "num_converged", text = "Minimum value cannot be greater than maximum value.")
          showFeedbackDanger(inputId = "max_mle_runs", text = "Minimum value cannot be greater than maximum value.")
          bool_correct_advanced_2 <- FALSE
        }

        rv_deconvo_setting$bool_correct_advanced_2 <- bool_correct_advanced_2
      })

      observe({
        req(rv_deconvo_setting)

        if(all(rv_deconvo_setting$bool_correct_advanced_temp,
               rv_deconvo_setting$bool_correct_advanced_rep_ae,
               rv_deconvo_setting$bool_correct_advanced_deg,
               rv_deconvo_setting$bool_correct_advanced_1,
               rv_deconvo_setting$bool_correct_advanced_2)){
          enable("act_advance_deconvo_save")
        }else{
          disable("act_advance_deconvo_save")
        }
      })

      observeEvent(input$act_advance_deconvo_save, {
        rv_deconvo_setting$min_temp <- input$min_temp
        rv_deconvo_setting$max_temp <- input$max_temp
        rv_deconvo_setting$min_rep_ae <- input$min_rep_ae
        rv_deconvo_setting$max_rep_ae <- input$max_rep_ae
        rv_deconvo_setting$min_deg <- input$min_deg
        rv_deconvo_setting$stop_iter <- input$stop_iter
        rv_deconvo_setting$max_iter <- input$max_iter
        rv_deconvo_setting$num_converged <- input$num_converged
        rv_deconvo_setting$max_mle_runs <- input$max_mle_runs
        removeModal()
      })

      # --- Perform deconvolution --- #

      observeEvent(input$act_analyze_deconvo, {
        waiter_ui$show()

        # Fix setting
        rv_deconvo_setting$hnc_from <- hnc_from <- input$noc_from
        rv_deconvo_setting$hnc_to <- hnc_to <- input$noc_to
        rv_deconvo_setting$theta <- theta <- input$theta

        fn_csp <- rv_deconvo_setting$fn_csp
        dt_csp <- rv_deconvo_setting$dt_csp
        fn_af <- rv_deconvo_setting$fn_af
        dt_af <- rv_deconvo_setting$dt_af
        fn_anameth <- rv_deconvo_setting$fn_anameth
        dt_anameth <- rv_deconvo_setting$dt_anameth
        min_temp <- rv_deconvo_setting$min_temp
        max_temp <- rv_deconvo_setting$max_temp
        min_rep_ae <- rv_deconvo_setting$min_rep_ae
        max_rep_ae <- rv_deconvo_setting$max_rep_ae
        min_deg <- rv_deconvo_setting$min_deg
        stop_iter <- rv_deconvo_setting$stop_iter
        max_iter <- rv_deconvo_setting$max_iter
        num_converged <- rv_deconvo_setting$num_converged
        max_mle_runs <- rv_deconvo_setting$max_mle_runs

        # Extract information on the kit
        dt_kit <- fread(paste0(path_pack, "/extdata/kit/", dt_anameth$Name[dt_anameth$Category == "Kit"], ".csv"))
        dt_kit <- try(change_data_type_kit(dt_kit))
        if(inherits(dt_kit, "try-error")){
          showModal(modalDialog(title = "Error", "The file format of the selected kit is invalid. Please update a file in the correct format.", easyClose = TRUE, footer = NULL))
        }else{
          dt_kit_auto <- extract_kit_auto(dt_kit)
          markers_auto <- unique(dt_kit_auto[, Marker])

          rv_deconvo_setting$dt_ignore_mk <- dt_ignore_mk <- hot_to_r(input$dt_ignore_mk)

          # Check input files
          error_message <- check_main_file(dt_csp, NULL, dt_af, markers_auto, method = "deconvo")
          if(error_message != ""){
            waiter_ui$hide()
            showModal(modalDialog(title = "Error", HTML(error_message), easyClose = TRUE, footer = NULL))
          }else{
            # Check registration of the stutter ratio
            list_dt_average_sr <- create_list_dt_average_sr(dt_anameth, markers_auto)
            message_sr <- check_csp_sr(dt_csp, dt_ignore_mk, dt_anameth, list_dt_average_sr)

            if(message_sr == "error"){
              waiter_ui$hide()
              showModal(modalDialog(title = "Error", "All markers are ignored!", easyClose = TRUE, footer = NULL))
            }else{
              hnc_all <- hnc_from:hnc_to
              data_deconvo <- perform_deconvo_all_noc(dt_csp, dt_af, dt_kit_auto, dt_ignore_mk, dt_anameth, list_dt_average_sr, hnc_all, theta,
                                                      min_temp, max_temp, min_rep_ae, max_rep_ae, min_deg, stop_iter, max_iter, num_converged, max_mle_runs,
                                                      waiter_ui)
              data_deconvo$hnc_all <- hnc_all
              data_deconvo$theta <- theta
              data_deconvo$dt_csp <- dt_csp
              data_deconvo$dt_af <- dt_af
              data_deconvo$fn_csp <- fn_csp
              data_deconvo$fn_af <- fn_af
              data_deconvo$fn_anameth <- fn_anameth
              data_deconvo$dt_anameth <- dt_anameth
              data_deconvo$dt_kit <- dt_kit
              data_deconvo$dt_ignore_mk <- dt_ignore_mk
              data_deconvo$list_dt_average_sr <- list_dt_average_sr
              data_deconvo$message_sr <- message_sr
              data_deconvo$min_temp <- min_temp
              data_deconvo$max_temp <- max_temp
              data_deconvo$min_rep_ae <- min_rep_ae
              data_deconvo$max_rep_ae <- max_rep_ae
              data_deconvo$min_deg <- min_deg
              data_deconvo$stop_iter <- stop_iter
              data_deconvo$max_iter <- max_iter
              data_deconvo$num_converged <- num_converged
              data_deconvo$max_mle_runs <- max_mle_runs
              rv_deconvo$data_deconvo <- data_deconvo

              waiter_ui$hide()
            }
          }
        }
      })

      return(rv_deconvo)
    }
  )
}
