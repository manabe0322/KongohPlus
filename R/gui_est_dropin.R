#' est_dropin_ui
#'
#' @description The function to create the ui module for estimating drop-in parameters
#' @param id The ui ID
est_dropin_ui <- function(id){
  ns <- NS(id)

  tabPanel("6. Estimate the drop-in parameters",
           titlePanel("6. Estimate the drop-in parameters"),
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
           br(),
           fluidRow(
             column(
               3,
               wellPanel(
                 fileInput(ns("ntc_data"), label = h4("Non template control data"), accept = ".csv", buttonLabel = "Load"),
                 br(),
                 h4("Kit"),
                 fluidRow(
                   column(3, br(), shinyFilesButton(ns("load_kit"), label = "Load", title = "Kit", multiple = FALSE)),
                   column(9, uiOutput(ns("name_kit_dropin")))
                 ),
                 br(),
                 numericInput(ns("at_dropin"), label = h4("Analytical threshold"), value = 35),
                 br(),
                 actionButton(ns("act_advance_dropin"), "Advanced setting"),
                 br(),
                 br(),
                 br(),
                 uiOutput(ns("ui_analyze_est_dropin")),
                 br(),
                 br(),
                 disabled(actionButton(ns("act_saveas_est_dropin_1"), label = "Save as", class = "btn btn-primary btn-lg"))
               )
             ),
             column(
               4,
               wellPanel(
                 h4("Estimated parameters"),
                 dataTableOutput(ns("dt_est_dropin")),
                 br(),
                 uiOutput(ns("button_all_mle"))
               ),
               wellPanel(
                 h4("Warnings"),
                 htmlOutput(ns("warning_message")),
                 br()
               )
             ),
             column(
               5,
               wellPanel(
                 h4("Experimental drop-in data"),
                 uiOutput(ns("ui_plot_dropin")),
                 br(),
                 dataTableOutput(ns("dt_dropin_data"))
               )
             )
           )
  )
}

#' est_dropin_server
#'
#' @description The function to create the server module for estimating drop-in parameters
#' @param id The server ID
#' @param path_pack The path of the package KongohPlus
est_dropin_server <- function(id, path_pack){
  moduleServer(
    id,
    function(input, output, session){
      rv_est_dropin_setting <- reactiveValues()
      rv_est_dropin_setting$dt_ntc <- NULL
      rv_est_dropin_setting$fn_ntc <- NULL
      rv_est_dropin_setting$dt_kit <- NULL
      rv_est_dropin_setting$name_kit_dropin <- NULL
      rv_est_dropin_setting$min_shape_di <- 0.01
      rv_est_dropin_setting$max_shape_di <- 100
      rv_est_dropin_setting$min_scale_di <- 0.01
      rv_est_dropin_setting$max_scale_di <- 100
      rv_est_dropin_setting$stop_iter <- 500
      rv_est_dropin_setting$max_iter <- 5000
      rv_est_dropin_setting$num_converged <- 3
      rv_est_dropin_setting$max_mle_runs <- 10
      rv_est_dropin_setting$max_diff_like <- 0.005
      rv_est_dropin_setting$bool_correct_mle_shape_di <- TRUE
      rv_est_dropin_setting$bool_correct_mle_scale_di <- TRUE
      rv_est_dropin_setting$bool_correct_advanced_1 <- TRUE
      rv_est_dropin_setting$bool_correct_advanced_2 <- TRUE
      rv_est_dropin_setting$bool_correct_advanced_3 <- TRUE
      rv_est_dropin_setting$data_est_dropin <- NULL

      rv_est_dropin_result <- reactiveValues()
      rv_est_dropin_result$warning_message <- NULL
      rv_est_dropin_result$dt_est_dropin <- NULL
      rv_est_dropin_result$dt_dropin_data <- NULL
      rv_est_dropin_result$opt_list <- NULL
      rv_est_dropin_result$index_converged <- NULL
      rv_est_dropin_result$min_shape_di <- NULL
      rv_est_dropin_result$max_shape_di <- NULL
      rv_est_dropin_result$min_scale_di <- NULL
      rv_est_dropin_result$max_scale_di <- NULL

      # Load non template control data
      observeEvent(input$ntc_data, {
        file_input <- input$ntc_data
        if(!is.null(file_input)){
          dt_ntc <- fread(file_input$datapath)
          dt_ntc <- try(change_data_type_csp(dt_ntc))
          if(!inherits(dt_ntc, "try-error")){
            rv_est_dropin_setting$dt_ntc <- dt_ntc
            rv_est_dropin_setting$fn_ntc <- file_input$name
          }else{
            rv_est_dropin_setting$dt_ntc <- NULL
            rv_est_dropin_setting$fn_ntc <- NULL
            showModal(modalDialog(title = "Error", "The selected file format is invalid. Please upload a file in the correct format.", easyClose = TRUE, footer = NULL))
          }
        }
      })

      # Select a kit
      shinyFileChoose(input, "load_kit", roots = c(documents = paste0(path_pack, "/extdata/kit")), session = session, filetypes = c("csv"))

      observeEvent(input$load_kit, {
        req(input$load_kit)
        path_kit <- parseFilePaths(roots = c(documents = paste0(path_pack, "/extdata/kit")), input$load_kit)$datapath
        if(length(path_kit) > 0){
          dt_kit <- fread(path_kit)
          dt_kit <- try(change_data_type_kit(dt_kit))
          if(inherits(dt_kit, "try-error")){
            rv_est_dropin_setting$dt_kit <- NULL
            showModal(modalDialog(title = "Error", "The file format of the selected kit is invalid. Please update a file in the correct format.", easyClose = TRUE, footer = NULL))
          }else{
            rv_est_dropin_setting$dt_kit <- dt_kit
            rv_est_dropin_setting$name_kit_dropin <- str_replace(basename(path_kit), ".{4}$", "")
          }
        }
      })

      output$name_kit_dropin <- renderUI({
        disabled(textInput(session$ns("name_kit_dropin"), label = "", value = rv_est_dropin_setting$name_kit_dropin))
      })

      # --- Advanced setting --- #

      observeEvent(input$act_advance_dropin, {
        showModal(modalDialog(
          title = "MLE conditions",
          fluidRow(
            column(6,
                   numericInput(session$ns("stop_iter"), label = "Minimum number of iterations", value = rv_est_dropin_setting$stop_iter),
                   numericInput(session$ns("max_iter"), label = "Maximum number of iterations", value = rv_est_dropin_setting$max_iter),
                   numericInput(session$ns("num_converged"), label = "Minimum number of MLE attempts", rv_est_dropin_setting$num_converged),
                   numericInput(session$ns("max_mle_runs"), label = "Maximum number of MLE attempts", value = rv_est_dropin_setting$max_mle_runs),
                   numericInput(session$ns("max_diff_like"), label = "Reproducibility threshold (log10-likelihood difference)", value = rv_est_dropin_setting$max_diff_like)
            ),
            column(6,
                   numericInput(session$ns("min_shape_di"), label = "Minimum shape parameter", value = rv_est_dropin_setting$min_shape_di),
                   numericInput(session$ns("min_scale_di"), label = "Minimum scale parameter", value = rv_est_dropin_setting$min_scale_di),
                   numericInput(session$ns("max_shape_di"), label = "Maximum shape parameter", value = rv_est_dropin_setting$max_shape_di),
                   numericInput(session$ns("max_scale_di"), label = "Maximum scale parameter", value = rv_est_dropin_setting$max_scale_di)
            )
          ),
          footer = tagList(
            actionButton(session$ns("act_advance_dropin_save"), "Save"),
            modalButton("Cancel")
          ),
          size = "l"
        ))
      })

      observeEvent(input$min_shape_di | input$max_shape_di, {
        min_shape_di <- input$min_shape_di
        max_shape_di <- input$max_shape_di

        hideFeedback("min_shape_di")
        hideFeedback("max_shape_di")
        rv_est_dropin_setting$bool_correct_mle_shape_di <- TRUE

        if(!is.numeric(min_shape_di) || min_shape_di <= 0){
          showFeedbackDanger(inputId = "min_shape_di", text = "Value must be positive.")
          rv_est_dropin_setting$bool_correct_mle_shape_di <- FALSE
        }

        if(!is.numeric(max_shape_di) || max_shape_di <= 0){
          showFeedbackDanger(inputId = "max_shape_di", text = "Value must be positive.")
          rv_est_dropin_setting$bool_correct_mle_shape_di <- FALSE
        }

        if(is.numeric(min_shape_di) && is.numeric(max_shape_di) && max_shape_di < min_shape_di){
          showFeedbackDanger(inputId = "min_shape_di", text = "Minimum value cannot be greater than maximum value.")
          showFeedbackDanger(inputId = "max_shape_di", text = "Minimum value cannot be greater than maximum value.")
          rv_est_dropin_setting$bool_correct_mle_shape_di <- FALSE
        }
      })

      observeEvent(input$min_scale_di | input$max_scale_di, {
        min_scale_di <- input$min_scale_di
        max_scale_di <- input$max_scale_di

        hideFeedback("min_scale_di")
        hideFeedback("max_scale_di")
        rv_est_dropin_setting$bool_correct_mle_scale_di <- TRUE

        if(!is.numeric(min_scale_di) || min_scale_di <= 0){
          showFeedbackDanger(inputId = "min_scale_di", text = "Value must be positive.")
          rv_est_dropin_setting$bool_correct_mle_scale_di <- FALSE
        }

        if(!is.numeric(max_scale_di) || max_scale_di <= 0){
          showFeedbackDanger(inputId = "max_scale_di", text = "Value must be positive.")
          rv_est_dropin_setting$bool_correct_mle_scale_di <- FALSE
        }

        if(is.numeric(min_scale_di) && is.numeric(max_scale_di) && max_scale_di < min_scale_di){
          showFeedbackDanger(inputId = "min_scale_di", text = "Minimum value cannot be greater than maximum value.")
          showFeedbackDanger(inputId = "max_scale_di", text = "Minimum value cannot be greater than maximum value.")
          rv_est_dropin_setting$bool_correct_mle_scale_di <- FALSE
        }
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

        rv_est_dropin_setting$bool_correct_advanced_1 <- bool_correct_advanced_1
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

        rv_est_dropin_setting$bool_correct_advanced_2 <- bool_correct_advanced_2
      })

      observeEvent(input$max_diff_like, {
        max_diff_like <- input$max_diff_like

        hideFeedback("max_diff_like")
        bool_correct_advanced_3 <- TRUE

        if(!is.numeric(max_diff_like) || max_diff_like < 0){
          showFeedbackDanger(inputId = "max_diff_like", text = "Value must be a numeric value greater than zero.")
          bool_correct_advanced_3 <- FALSE
        }

        rv_est_dropin_setting$bool_correct_advanced_3 <- bool_correct_advanced_3
      })

      observe({
        req(rv_est_dropin_setting)

        if(all(rv_est_dropin_setting$bool_correct_mle_shape_di,
               rv_est_dropin_setting$bool_correct_mle_scale_di,
               rv_est_dropin_setting$bool_correct_advanced_1,
               rv_est_dropin_setting$bool_correct_advanced_2,
               rv_est_dropin_setting$bool_correct_advanced_3)){
          enable("act_advance_dropin_save")
        }else{
          disable("act_advance_dropin_save")
        }
      })

      observeEvent(input$act_advance_dropin_save, {
        rv_est_dropin_setting$min_shape_di <- input$min_shape_di
        rv_est_dropin_setting$max_shape_di <- input$max_shape_di
        rv_est_dropin_setting$min_scale_di <- input$min_scale_di
        rv_est_dropin_setting$max_scale_di <- input$max_scale_di
        rv_est_dropin_setting$stop_iter <- input$stop_iter
        rv_est_dropin_setting$max_iter <- input$max_iter
        rv_est_dropin_setting$num_converged <- input$num_converged
        rv_est_dropin_setting$max_mle_runs <- input$max_mle_runs
        rv_est_dropin_setting$max_diff_like <- input$max_diff_like
        removeModal()
      })

      output$ui_analyze_est_dropin <- renderUI({
        dt_ntc <- rv_est_dropin_setting$dt_ntc
        dt_kit <- rv_est_dropin_setting$dt_kit
        if(any(is.null(dt_ntc), is.null(dt_kit))){
          disabled(actionButton(session$ns("act_analyze_est_dropin"), label = "Analyze", class = "btn btn-primary btn-lg"))
        }else{
          actionButton(session$ns("act_analyze_est_dropin"), label = "Analyze", class = "btn btn-primary btn-lg")
        }
      })

      # Start estimating dropin parameters
      observeEvent(input$act_analyze_est_dropin, {
        req(rv_est_dropin_setting$dt_ntc)
        req(rv_est_dropin_setting$fn_ntc)
        req(rv_est_dropin_setting$name_kit_dropin)
        req(rv_est_dropin_setting$min_shape_di)
        req(rv_est_dropin_setting$max_shape_di)
        req(rv_est_dropin_setting$min_scale_di)
        req(rv_est_dropin_setting$max_scale_di)
        req(rv_est_dropin_setting$stop_iter)
        req(rv_est_dropin_setting$max_iter)
        req(rv_est_dropin_setting$num_converged)
        req(rv_est_dropin_setting$max_mle_runs)
        req(rv_est_dropin_setting$max_diff_like)
        req(rv_est_dropin_setting$dt_kit)

        dt_ntc <- rv_est_dropin_setting$dt_ntc
        fn_ntc <- rv_est_dropin_setting$fn_ntc
        name_kit_dropin <- rv_est_dropin_setting$name_kit_dropin
        min_shape_di <- rv_est_dropin_setting$min_shape_di
        max_shape_di <- rv_est_dropin_setting$max_shape_di
        min_scale_di <- rv_est_dropin_setting$min_scale_di
        max_scale_di <- rv_est_dropin_setting$max_scale_di
        stop_iter <- rv_est_dropin_setting$stop_iter
        max_iter <- rv_est_dropin_setting$max_iter
        num_converged <- rv_est_dropin_setting$num_converged
        max_mle_runs <- rv_est_dropin_setting$max_mle_runs
        max_diff_like <- rv_est_dropin_setting$max_diff_like
        dt_kit <- rv_est_dropin_setting$dt_kit
        at_dropin <- input$at_dropin

        waiter_ui <- Waiter$new(html = tagList(spin_3k(),
                                               h3("Checking the settings...")),
                                color = transparent(.5))
        waiter_ui$show()

        # Check input files and setting
        error_message <- check_est_dropin(dt_ntc, dt_kit)
        if(error_message != ""){
          waiter_ui$hide()
          showModal(modalDialog(title = "Error", HTML(error_message), easyClose = TRUE, footer = NULL))
        }else{
          data_est_dropin <- try(estimate_dropin_params(dt_ntc, dt_kit, at_dropin, min_shape_di, max_shape_di, min_scale_di, max_scale_di, stop_iter, max_iter, num_converged, max_mle_runs, max_diff_like, waiter_ui))
          if(inherits(data_est_dropin, "try-error")){
            waiter_ui$hide()
            showModal(modalDialog(title = "Error", "An unexpected error has occurred.", easyClose = TRUE, footer = NULL))
          }else{
            data_est_dropin$fn_ntc <- fn_ntc
            data_est_dropin$name_kit_dropin <- name_kit_dropin
            data_est_dropin$at_dropin <- at_dropin
            data_est_dropin$min_shape_di <- min_shape_di
            data_est_dropin$max_shape_di <- max_shape_di
            data_est_dropin$min_scale_di <- min_scale_di
            data_est_dropin$max_scale_di <- max_scale_di
            data_est_dropin$stop_iter <- stop_iter
            data_est_dropin$max_iter <- max_iter
            data_est_dropin$num_converged <- num_converged
            data_est_dropin$max_mle_runs <- max_mle_runs
            data_est_dropin$max_diff_like <- max_diff_like
            rv_est_dropin_setting$data_est_dropin <- data_est_dropin

            waiter_ui$hide()
          }
        }
      })

      observe({
        req(rv_est_dropin_setting)
        data_est_dropin <- rv_est_dropin_setting$data_est_dropin
        if(!is.null(data_est_dropin)){
          rv_est_dropin_result$warning_message <- data_est_dropin$warning_message
          rv_est_dropin_result$dt_est_dropin <- data_est_dropin$dt_est_dropin
          rv_est_dropin_result$dt_dropin_data <- data_est_dropin$dt_dropin_data
          rv_est_dropin_result$opt_list <- data_est_dropin$opt_list
          rv_est_dropin_result$index_converged <- data_est_dropin$index_converged
          rv_est_dropin_result$min_shape_di <- data_est_dropin$min_shape_di
          rv_est_dropin_result$max_shape_di <- data_est_dropin$max_shape_di
          rv_est_dropin_result$min_scale_di <- data_est_dropin$min_scale_di
          rv_est_dropin_result$max_scale_di <- data_est_dropin$max_scale_di

          enable("act_saveas_est_dropin_1")
        }else{
          rv_est_dropin_result$warning_message <- NULL
          rv_est_dropin_result$dt_est_dropin <- NULL
          rv_est_dropin_result$dt_dropin_data <- NULL
          rv_est_dropin_result$opt_list <- NULL
          rv_est_dropin_result$index_converged <- NULL
          rv_est_dropin_result$min_shape_di <- NULL
          rv_est_dropin_result$max_shape_di <- NULL
          rv_est_dropin_result$min_scale_di <- NULL
          rv_est_dropin_result$max_scale_di <- NULL

          disable("act_saveas_est_dropin_1")
        }
      })

      observeEvent(input$act_saveas_est_dropin_1, {
        req(rv_est_dropin_setting)
        data_est_dropin <- rv_est_dropin_setting$data_est_dropin

        if(!is.null(data_est_dropin)){
          showModal(modalDialog(
            title = "Save drop-in data",
            textInput(session$ns("fn_est_dropin"), label = "Enter the file name"),
            footer = tagList(
              disabled(actionButton(session$ns("act_saveas_est_dropin_2"), "Save")),
              modalButton("Close")
            )
          ))
        }

        else{
          showModal(modalDialog(title = "Error", "No drop-in data.", easyClose = TRUE, footer = NULL))
        }
      })

      observeEvent(input$fn_est_dropin, {
        fn_est_dropin <- input$fn_est_dropin
        dropin_file_names <- list.files(paste0(path_pack, "/extdata/dropin"))
        if(is.element(paste0(fn_est_dropin, ".RData"), dropin_file_names)){
          showFeedbackDanger(inputId = "fn_est_dropin", text = "A file with the same name already exists.")
          disable("act_saveas_est_dropin_2")
        }else if(nchar(fn_est_dropin) == 0){
          hideFeedback("fn_est_dropin")
          disable("act_saveas_est_dropin_2")
        }else{
          hideFeedback("fn_est_dropin")
          enable("act_saveas_est_dropin_2")
        }
      })

      observeEvent(input$act_saveas_est_dropin_2, {
        req(rv_est_dropin_setting)
        data_est_dropin <- rv_est_dropin_setting$data_est_dropin
        fn_est_dropin <- input$fn_est_dropin

        if(!is.null(data_est_dropin)){
          save(list = c("data_est_dropin"),
               file = paste0(path_pack, "/extdata/dropin/", fn_est_dropin, ".RData"))

          removeModal()
        }

        else{
          removeModal()
          showModal(modalDialog(title = "Error", "No drop-in data.", easyClose = TRUE, footer = NULL))
        }
      })

      output$dt_est_dropin <- renderDataTable({
        req(rv_est_dropin_result$dt_est_dropin)
        dt_est_dropin <- rv_est_dropin_result$dt_est_dropin

        datatable(
          dt_est_dropin,
          selection = "none",
          options = list(dom = "t", iDisplayLength = nrow(dt_est_dropin), ordering = FALSE, autowidth = FALSE),
          rownames = FALSE
        )
      })

      output$button_all_mle <- renderUI({
        req(rv_est_dropin_result$dt_est_dropin)
        actionButton(session$ns("act_all_mle"), label = "View all MLE attempts")
      })

      observeEvent(input$act_all_mle, {
        showModal(modalDialog(
          title = "All MLE attempts",
          dataTableOutput(session$ns("dt_all_mle")),
          easyClose = TRUE,
          footer = NULL
        ) %>% tagAppendAttributes(class = "edit-modal")
        )
      })

      output$dt_all_mle <- renderDataTable({
        req(rv_est_dropin_result$opt_list)
        req(rv_est_dropin_result$index_converged)
        req(rv_est_dropin_result$min_shape_di)
        req(rv_est_dropin_result$max_shape_di)
        req(rv_est_dropin_result$min_scale_di)
        req(rv_est_dropin_result$max_scale_di)

        opt_list <- rv_est_dropin_result$opt_list
        index_converged <- rv_est_dropin_result$index_converged
        min_shape_di <- rv_est_dropin_result$min_shape_di
        max_shape_di <- rv_est_dropin_result$max_shape_di
        min_scale_di <- rv_est_dropin_result$min_scale_di
        max_scale_di <- rv_est_dropin_result$max_scale_di

        n_param <- length(opt_list[[1]][["par"]])
        count_mle <- length(opt_list)

        mat_all_mle <- matrix(NA, n_param + 1, count_mle)
        all_scores <- rep(NA, count_mle)

        for(i in 1:count_mle){
          opt_one <- opt_list[[i]]
          all_scores[i] <- score <- opt_one$value
          est_params <- opt_one$par

          shape <- est_params[1] * (max_shape_di - min_shape_di) + min_shape_di
          scale <- est_params[2] * (max_scale_di - min_scale_di) + min_scale_di

          mat_all_mle[1, i] <- - score
          mat_all_mle[2:(n_param + 1), i] <- c(shape, scale)
        }

        dt_all_mle <- as.data.table(mat_all_mle)
        names(dt_all_mle) <- paste0("MLE ", 1:count_mle)

        name_param <- c("Shape parameter", "Scale parameter")
        dt_all_mle[, Parameter := c("Log10 (likelihood)", name_param)]
        setcolorder(dt_all_mle, c("Parameter", setdiff(names(dt_all_mle), "Parameter")))

        name_converged <- paste0("MLE ", index_converged)

        datatable(
          dt_all_mle,
          selection = "none",
          options = list(dom = "t", iDisplayLength = nrow(dt_all_mle), ordering = FALSE, autowidth = TRUE),
          rownames = FALSE
        ) %>%
          formatSignif(columns = paste0("MLE ", 1:count_mle), digits = 6) %>%
          formatStyle(names(dt_all_mle), lineHeight = "80%") %>%
          formatStyle(name_converged, backgroundColor = "lightblue")
      })

      output$warning_message <- renderUI({
        req(rv_est_dropin_result$warning_message)
        HTML(rv_est_dropin_result$warning_message)
      })

      output$ui_plot_dropin <- renderUI({
        req(rv_est_dropin_result$dt_dropin_data)
        plotOutput(session$ns("plot_dropin"), width = "100%", height = "300px")
      })

      output$plot_dropin <- renderPlot({
        req(rv_est_dropin_result$dt_est_dropin)
        req(rv_est_dropin_result$dt_dropin_data)
        dt_est_dropin <- rv_est_dropin_result$dt_est_dropin
        dt_dropin_data <- rv_est_dropin_result$dt_dropin_data

        heights <- dt_dropin_data$Height
        shape_di <- dt_est_dropin$Value[dt_est_dropin$Parameter == "Shape parameter"]
        scale_di <- dt_est_dropin$Value[dt_est_dropin$Parameter == "Scale parameter"]

        tb_ph <- table(heights)
        x_val <- 0:round(max(heights) * 1.5, 0)
        y_val <- rep(0, length(x_val))
        y_val[match(as.numeric(names(tb_ph)), x_val)] <- tb_ph

        ggplot() +
          geom_line(data = data.frame(X = x_val, Y = y_val, P = y_val * 0.01), aes(X, P, color = "Observed"), linewidth = 1) +
          stat_function(fun = dgamma, args = list(shape = shape_di, scale = scale_di), aes(color = "Estimated"), linewidth = 1) +
          xlim(0, max(x_val)) +
          labs(color = NULL, x = "Peak height (RFU)") +
          scale_y_continuous(
            name = "Probability density",
            sec.axis = sec_axis(~ . * 100, name = "Count")
          ) +
          scale_color_manual(values = c("Observed" = "#4169e1", "Estimated" = "#c71585")) +
          theme(
            plot.title = element_text(size = 14),
            legend.text = element_text(size = 12),
            axis.title.x = element_text(size = 12),
            axis.title.y = element_text(size = 12),
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12)
          )
      })

      output$dt_dropin_data <- renderDataTable({
        req(rv_est_dropin_result$dt_dropin_data)
        dt_dropin_data <- rv_est_dropin_result$dt_dropin_data

        datatable(
          dt_dropin_data,
          selection = "none",
          options = list(
            iDisplayLength = 10,
            autoWidth = FALSE
          ),
          rownames = FALSE
        ) %>%
          formatStyle(names(dt_dropin_data), lineHeight = "80%")
      })
    }
  )
}
