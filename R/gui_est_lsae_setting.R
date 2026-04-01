#' est_lsae_setting_ui
#'
#' @description The function to create the ui module for setting to estimate locus-specific amplification efficiencies
#' @param id The ui ID
est_lsae_setting_ui <- function(id){
  ns <- NS(id)

  tabPanel("Setting",
           br(),
           sidebarLayout(
             sidebarPanel(width = 3,
                          fileInput(ns("exp_data"), label = h4("Experimental data"), accept = ".csv", buttonLabel = "Load"),
                          fileInput(ns("gt_data"), label = h4("Genotype data"), accept = ".csv", buttonLabel = "Load"),
                          h4("Stutter ratio"),
                          fluidRow(
                            column(3, br(), shinyFilesButton(ns("load_sr"), label = "Load", title = "Stutter ratio", multiple = FALSE)),
                            column(9, uiOutput(ns("name_sr")))
                          ),
                          br(),
                          h4("Kit"),
                          uiOutput(ns("name_kit_lsae")),
                          br(),
                          actionButton(ns("act_advance_lsae"), "Advanced setting"),
                          br(),
                          br(),
                          br(),
                          uiOutput(ns("ui_analyze_est_lsae"))
             ),
             mainPanel(width = 9,
                       tags$style(HTML("
                         .handsontable .htCore td {
                           vertical-align: middle !important;
                           text-align: center !important;
                         }"
                       )),
                       rHandsontableOutput(ns("dt_setting_est_par")),
                       br()
             )
           )
  )
}

#' est_lsae_setting_server
#'
#' @description The function to create the server module for setting to estimate locus-specific amplification efficiencies
#' @param id The server ID
#' @param path_pack The path of the package KongohPlus
#' @param session_lsae_top The top session for estimating the LSAE
est_lsae_setting_server <- function(id, path_pack, session_lsae_top){
  moduleServer(
    id,
    function(input, output, session){

      # --- Define the initial reactive values --- #

      rv_est_lsae_setting <- reactiveValues()
      rv_est_lsae_setting$dt_exp_lsae <- NULL
      rv_est_lsae_setting$fn_exp_lsae <- NULL
      rv_est_lsae_setting$dt_gt_lsae <- NULL
      rv_est_lsae_setting$fn_gt_lsae <- NULL
      rv_est_lsae_setting$dt_setting_est_par <- NULL
      rv_est_lsae_setting$list_dt_average_sr <- NULL
      rv_est_lsae_setting$name_kit_lsae <- NULL
      rv_est_lsae_setting$name_sr <- NULL
      rv_est_lsae_setting$min_alpha <- 0.01
      rv_est_lsae_setting$max_alpha <- 100
      rv_est_lsae_setting$max_iter <- 500
      rv_est_lsae_setting$num_consistent <- 3
      rv_est_lsae_setting$max_mle_runs <- 10
      rv_est_lsae_setting$max_diff_like <- 0.005
      rv_est_lsae_setting$bool_correct_mle_alpha <- TRUE
      rv_est_lsae_setting$bool_correct_advanced_1 <- TRUE
      rv_est_lsae_setting$bool_correct_advanced_2 <- TRUE
      rv_est_lsae_setting$bool_correct_advanced_3 <- TRUE

      rv_est_lsae <- reactiveValues()
      rv_est_lsae$data_est_lsae <- NULL

      # --- Load experimental data --- #

      observeEvent(input$exp_data, {
        file_input <- input$exp_data
        if(!is.null(file_input)){
          dt_exp_lsae <- fread(file_input$datapath)
          dt_exp_lsae <- try(change_data_type_csp(dt_exp_lsae))
          if(!inherits(dt_exp_lsae, "try-error")){
            rv_est_lsae_setting$dt_exp_lsae <- dt_exp_lsae
            rv_est_lsae_setting$fn_exp_lsae <- file_input$name
          }else{
            rv_est_lsae_setting$dt_exp_lsae <- NULL
            rv_est_lsae_setting$fn_exp_lsae <- NULL
            showModal(modalDialog(title = "Error", "The selected file format is invalid. Please upload a file in the correct format.", easyClose = TRUE, footer = NULL))
          }
        }
      })

      # --- Load genotype data --- #

      observeEvent(input$gt_data, {
        file_input <- input$gt_data
        if(!is.null(file_input)){
          dt_gt_lsae <- fread(file_input$datapath)
          dt_gt_lsae <- try(change_data_type_ref(dt_gt_lsae))
          if(!inherits(dt_gt_lsae, "try-error")){
            rv_est_lsae_setting$dt_gt_lsae <- dt_gt_lsae
            rv_est_lsae_setting$fn_gt_lsae <- file_input$name
          }else{
            rv_est_lsae_setting$dt_gt_lsae <- NULL
            rv_est_lsae_setting$fn_gt_lsae <- NULL
            showModal(modalDialog(title = "Error", "The selected file format is invalid. Please upload a file in the correct format.", easyClose = TRUE, footer = NULL))
          }
        }
      })

      # --- Load data for stutter ratios --- #

      shinyFileChoose(input, "load_sr", roots = c(documents = paste0(path_pack, "/extdata/stutter")), session = session, filetypes = c("RData"))

      observeEvent(input$load_sr, {
        req(input$load_sr)
        path_sr <- parseFilePaths(roots = c(documents = paste0(path_pack, "/extdata/stutter")), input$load_sr)$datapath
        if(length(path_sr) > 0){
          load(path_sr)
          rv_est_lsae_setting$list_dt_average_sr <- data_est_sr$list_dt_average_sr
          rv_est_lsae_setting$name_kit_lsae <- data_est_sr$name_kit_stutter
          rv_est_lsae_setting$name_sr <- str_replace(basename(path_sr), ".{6}$", "")
        }
      })

      output$name_sr <- renderUI({
        disabled(textInput(session$ns("name_sr"), label = "", value = rv_est_lsae_setting$name_sr))
      })

      # --- Prepare dataset for the kit --- #

      output$name_kit_lsae <- renderUI({
        disabled(textInput(session$ns("name_kit_lsae"), label = "", value = rv_est_lsae_setting$name_kit_lsae))
      })

      observeEvent(input$name_kit_lsae, ignoreInit = TRUE, {
        name_kit_lsae <- rv_est_lsae_setting$name_kit_lsae
        if(!is.null(name_kit_lsae)){
          dt_kit <- fread(paste0(path_pack, "/extdata/kit/", name_kit_lsae, ".csv"))
          dt_kit <- try(change_data_type_kit(dt_kit))
          if(inherits(dt_kit, "try-error")){
            rv_est_lsae_setting$dt_setting_est_par <- NULL
            showModal(modalDialog(title = "Error", "The file format of the selected kit is invalid. Please update a file in the correct format.", easyClose = TRUE, footer = NULL))
          }else{
            rv_est_lsae_setting$dt_setting_est_par <- create_init_dt_setting_est_par(dt_kit)
          }
        }
      })

      output$dt_setting_est_par <- renderRHandsontable({
        req(rv_est_lsae_setting$dt_setting_est_par)

        rhandsontable(rv_est_lsae_setting$dt_setting_est_par, rowHeaders = FALSE) %>%
          hot_table(fillHandle = FALSE, rowHeights = 30) %>%
          hot_col(c("Marker"), type = "text", readOnly = TRUE, width = 100) %>%
          hot_col(c("AT"), type = "numeric", format = "0", width = 50)
      })

      # --- Advanced setting --- #

      observeEvent(input$act_advance_lsae, {
        showModal(modalDialog(
          title = "MLE conditions",
          fluidRow(
            column(6,
                   numericInput(session$ns("max_iter"), label = "Maximum number of iterations", value = rv_est_lsae_setting$max_iter),
                   numericInput(session$ns("max_mle_runs"), label = "Maximum number of MLE attempts", value = rv_est_lsae_setting$max_mle_runs),
                   numericInput(session$ns("num_consistent"), label = "Minimum number of consistent runs", value = rv_est_lsae_setting$num_consistent),
                   numericInput(session$ns("max_diff_like"), label = "Reproducibility threshold (log10-likelihood difference)", value = rv_est_lsae_setting$max_diff_like)
            ),
            column(6,
                   numericInput(session$ns("min_alpha"), label = "Minimum alpha for LSAE", value = rv_est_lsae_setting$min_alpha),
                   numericInput(session$ns("max_alpha"), label = "Maximum alpha for LSAE", value = rv_est_lsae_setting$max_alpha)
            )
          ),
          footer = tagList(
            actionButton(session$ns("act_advance_lsae_save"), "Save"),
            modalButton("Cancel")
          ),
          size = "l"
        ))
      })

      observeEvent(input$min_alpha | input$max_alpha, {
        min_alpha <- input$min_alpha
        max_alpha <- input$max_alpha

        hideFeedback("min_alpha")
        hideFeedback("max_alpha")
        rv_est_lsae_setting$bool_correct_mle_alpha <- TRUE

        if(!is.numeric(min_alpha) || min_alpha <= 0){
          showFeedbackDanger(inputId = "min_alpha", text = "Value must be a positive number.")
          rv_est_lsae_setting$bool_correct_mle_alpha <- FALSE
        }

        if(!is.numeric(max_alpha) || max_alpha <= 0){
          showFeedbackDanger(inputId = "max_alpha", text = "Value must be a positive number.")
          rv_est_lsae_setting$bool_correct_mle_alpha <- FALSE
        }

        if(is.numeric(min_alpha) && is.numeric(max_alpha) && max_alpha < min_alpha){
          showFeedbackDanger(inputId = "min_alpha", text = "The minimum value cannot exceed the maximum value.")
          showFeedbackDanger(inputId = "max_alpha", text = "The minimum value cannot exceed the maximum value.")
          rv_est_lsae_setting$bool_correct_mle_alpha <- FALSE
        }
      })

      observeEvent(input$max_iter, {
        max_iter <- input$max_iter

        hideFeedback("max_iter")
        bool_correct_advanced_1 <- TRUE

        if(!is.integer(max_iter) || max_iter <= 0){
          showFeedbackDanger(inputId = "max_iter", text = "Value must be a positive integer.")
          bool_correct_advanced_1 <- FALSE
        }

        rv_est_lsae_setting$bool_correct_advanced_1 <- bool_correct_advanced_1
      })

      observeEvent(input$num_consistent | input$max_mle_runs, {
        num_consistent <- input$num_consistent
        max_mle_runs <- input$max_mle_runs

        hideFeedback("num_consistent")
        hideFeedback("max_mle_runs")
        bool_correct_advanced_2 <- TRUE

        if(!is.integer(num_consistent) || num_consistent <= 0){
          showFeedbackDanger(inputId = "num_consistent", text = "Value must be a positive integer.")
          bool_correct_advanced_2 <- FALSE
        }

        if(!is.integer(max_mle_runs) || max_mle_runs <= 0){
          showFeedbackDanger(inputId = "max_mle_runs", text = "Value must be a positive integer.")
          bool_correct_advanced_2 <- FALSE
        }

        if(is.integer(num_consistent) && is.integer(max_mle_runs) && num_consistent > max_mle_runs){
          showFeedbackDanger(inputId = "num_consistent", text = "The minimum number of consistent runs cannot exceed the number of MLE attempts.")
          showFeedbackDanger(inputId = "max_mle_runs", text = "The minimum number of consistent runs cannot exceed the number of MLE attempts.")
          bool_correct_advanced_2 <- FALSE
        }

        rv_est_lsae_setting$bool_correct_advanced_2 <- bool_correct_advanced_2
      })

      observeEvent(input$max_diff_like, {
        max_diff_like <- input$max_diff_like

        hideFeedback("max_diff_like")
        bool_correct_advanced_3 <- TRUE

        if(!is.numeric(max_diff_like) || max_diff_like < 0){
          showFeedbackDanger(inputId = "max_diff_like", text = "Value must be a positive number.")
          bool_correct_advanced_3 <- FALSE
        }

        rv_est_lsae_setting$bool_correct_advanced_3 <- bool_correct_advanced_3
      })

      observe({
        req(rv_est_lsae)

        if(all(rv_est_lsae_setting$bool_correct_mle_alpha,
               rv_est_lsae_setting$bool_correct_advanced_1,
               rv_est_lsae_setting$bool_correct_advanced_2,
               rv_est_lsae_setting$bool_correct_advanced_3)){
          enable("act_advance_lsae_save")
        }else{
          disable("act_advance_lsae_save")
        }
      })

      observeEvent(input$act_advance_lsae_save, {
        rv_est_lsae_setting$min_alpha <- input$min_alpha
        rv_est_lsae_setting$max_alpha <- input$max_alpha
        rv_est_lsae_setting$max_iter <- input$max_iter
        rv_est_lsae_setting$num_consistent <- input$num_consistent
        rv_est_lsae_setting$max_mle_runs <- input$max_mle_runs
        rv_est_lsae_setting$max_diff_like <- input$max_diff_like
        removeModal()
      })

      # --- Check whether the required data has been loaded --- #

      output$ui_analyze_est_lsae <- renderUI({
        dt_exp_lsae <- rv_est_lsae_setting$dt_exp_lsae
        dt_gt_lsae <- rv_est_lsae_setting$dt_gt_lsae
        dt_setting_est_par <- rv_est_lsae_setting$dt_setting_est_par
        if(any(is.null(dt_exp_lsae), is.null(dt_gt_lsae), is.null(dt_setting_est_par))){
          disabled(actionButton(session$ns("act_analyze_est_lsae"), label = "Analyze", class = "btn btn-primary btn-lg"))
        }else{
          actionButton(session$ns("act_analyze_est_lsae"), label = "Analyze", class = "btn btn-primary btn-lg")
        }
      })

      # --- Start estimating locus-specific amplification efficiencies --- #

      observeEvent(input$act_analyze_est_lsae, {
        waiter_ui <- Waiter$new(html = tagList(spin_3k(),
                                               h3("Checking the settings...")),
                                color = transparent(.5))
        waiter_ui$show()

        # Fix data
        fn_exp_lsae <- rv_est_lsae_setting$fn_exp_lsae
        dt_exp_lsae <- rv_est_lsae_setting$dt_exp_lsae
        fn_gt_lsae <- rv_est_lsae_setting$fn_gt_lsae
        dt_gt_lsae <- rv_est_lsae_setting$dt_gt_lsae
        dt_setting_est_par <- hot_to_r(input$dt_setting_est_par)
        list_dt_average_sr <- rv_est_lsae_setting$list_dt_average_sr
        name_kit_lsae <- rv_est_lsae_setting$name_kit_lsae
        name_sr <- rv_est_lsae_setting$name_sr
        min_alpha <- rv_est_lsae_setting$min_alpha
        max_alpha <- rv_est_lsae_setting$max_alpha
        max_iter <- rv_est_lsae_setting$max_iter
        num_consistent <- rv_est_lsae_setting$num_consistent
        max_mle_runs <- rv_est_lsae_setting$max_mle_runs
        max_diff_like <- rv_est_lsae_setting$max_diff_like

        # Check input files and setting
        error_message <- check_est_par(dt_exp_lsae, dt_gt_lsae, dt_setting_est_par)
        if(error_message != ""){
          waiter_ui$hide()
          showModal(modalDialog(title = "Error", HTML(error_message), easyClose = TRUE, footer = NULL))
        }else{
          data_est_lsae <- try(estimate_lsae(dt_exp_lsae, dt_gt_lsae, list_dt_average_sr, dt_setting_est_par,
                                             min_alpha, max_alpha,
                                             max_iter, num_consistent, max_mle_runs, max_diff_like,
                                             waiter_ui))
          if(inherits(data_est_lsae, "try-error")){
            waiter_ui$hide()
            showModal(modalDialog(title = "Error", "An unexpected error has occurred.", easyClose = TRUE, footer = NULL))
          }else{
            data_est_lsae$fn_exp_lsae <- fn_exp_lsae
            data_est_lsae$fn_gt_lsae <- fn_gt_lsae
            data_est_lsae$dt_setting_est_par <- dt_setting_est_par
            data_est_lsae$list_dt_average_sr <- list_dt_average_sr
            data_est_lsae$name_kit_lsae <- name_kit_lsae
            data_est_lsae$name_sr <- name_sr
            data_est_lsae$min_alpha <- min_alpha
            data_est_lsae$max_alpha <- max_alpha
            data_est_lsae$max_iter <- max_iter
            data_est_lsae$num_consistent <- num_consistent
            data_est_lsae$max_mle_runs <- max_mle_runs
            data_est_lsae$max_diff_like <- max_diff_like

            rv_est_lsae$data_est_lsae <- data_est_lsae

            waiter_ui$hide()

            updateTabsetPanel(session = session_lsae_top, inputId = "tabset_est_lsae", selected = "Result")
          }
        }
      })

      return(rv_est_lsae)
    }
  )
}
