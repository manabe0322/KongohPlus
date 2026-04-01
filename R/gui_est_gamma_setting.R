#' est_gamma_setting_ui
#'
#' @description The function to create the ui module for setting to estimate the gamma scale parameter
#' @param id The ui ID
est_gamma_setting_ui <- function(id){
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
                          h4("Locus-specific amplification efficiency"),
                          fluidRow(
                            column(3, br(), shinyFilesButton(ns("load_lsae"), label = "Load", title = "Locus-specific amplification efficiency", multiple = FALSE)),
                            column(9, uiOutput(ns("name_lsae")))
                          ),
                          br(),
                          h4("Kit"),
                          uiOutput(ns("name_kit_gamma")),
                          br(),
                          actionButton(ns("act_advance_gamma"), "Advanced setting"),
                          br(),
                          br(),
                          br(),
                          uiOutput(ns("ui_analyze_est_gamma"))
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

#' est_gamma_setting_server
#'
#' @description The function to create the server module for setting to estimate the gamma scale parameter
#' @param id The server ID
#' @param path_pack The path of the package KongohPlus
#' @param session_gamma_top The top session for estimating the gamma scale parameters
est_gamma_setting_server <- function(id, path_pack, session_gamma_top){
  moduleServer(
    id,
    function(input, output, session){

      # --- Define the initial reactive values --- #

      rv_est_gamma_setting <- reactiveValues()
      rv_est_gamma_setting$dt_exp_gamma <- NULL
      rv_est_gamma_setting$fn_exp_gamma <- NULL
      rv_est_gamma_setting$dt_gt_gamma <- NULL
      rv_est_gamma_setting$fn_gt_gamma <- NULL
      rv_est_gamma_setting$dt_setting_est_par <- NULL
      rv_est_gamma_setting$list_dt_average_sr <- NULL
      rv_est_gamma_setting$name_kit_stutter <- NULL
      rv_est_gamma_setting$name_sr <- NULL
      rv_est_gamma_setting$dt_expect_ae <- NULL
      rv_est_gamma_setting$name_kit_lsae <- NULL
      rv_est_gamma_setting$name_lsae <- NULL
      rv_est_gamma_setting$name_kit_gamma <- NULL
      rv_est_gamma_setting$min_scale_allele <- 0.01
      rv_est_gamma_setting$max_scale_allele <- 400
      rv_est_gamma_setting$min_scale_stutter <- 0.01
      rv_est_gamma_setting$max_scale_stutter <- 400
      rv_est_gamma_setting$max_iter <- 500
      rv_est_gamma_setting$num_consistent <- 3
      rv_est_gamma_setting$max_mle_runs <- 10
      rv_est_gamma_setting$max_diff_like <- 0.005
      rv_est_gamma_setting$bool_correct_mle_scale_allele <- TRUE
      rv_est_gamma_setting$bool_correct_mle_scale_stutter <- TRUE
      rv_est_gamma_setting$bool_correct_advanced_1 <- TRUE
      rv_est_gamma_setting$bool_correct_advanced_2 <- TRUE
      rv_est_gamma_setting$bool_correct_advanced_3 <- TRUE

      rv_est_gamma <- reactiveValues()
      rv_est_gamma$data_est_gamma <- NULL

      # --- Load experimental data --- #

      observeEvent(input$exp_data, {
        file_input <- input$exp_data
        if(!is.null(file_input)){
          dt_exp_gamma <- fread(file_input$datapath)
          dt_exp_gamma <- try(change_data_type_csp(dt_exp_gamma))
          if(!inherits(dt_exp_gamma, "try-error")){
            rv_est_gamma_setting$dt_exp_gamma <- dt_exp_gamma
            rv_est_gamma_setting$fn_exp_gamma <- file_input$name
          }else{
            rv_est_gamma_setting$dt_exp_gamma <- NULL
            rv_est_gamma_setting$fn_exp_gamma <- NULL
            showModal(modalDialog(title = "Error", "The selected file format is invalid. Please upload a file in the correct format.", easyClose = TRUE, footer = NULL))
          }
        }
      })

      # --- Load genotype data --- #

      observeEvent(input$gt_data, {
        file_input <- input$gt_data
        if(!is.null(file_input)){
          dt_gt_gamma <- fread(file_input$datapath)
          dt_gt_gamma <- try(change_data_type_ref(dt_gt_gamma))
          if(!inherits(dt_gt_gamma, "try-error")){
            rv_est_gamma_setting$dt_gt_gamma <- dt_gt_gamma
            rv_est_gamma_setting$fn_gt_gamma <- file_input$name
          }else{
            rv_est_gamma_setting$dt_gt_gamma <- NULL
            rv_est_gamma_setting$fn_gt_gamma <- NULL
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
          rv_est_gamma_setting$list_dt_average_sr <- data_est_sr$list_dt_average_sr
          rv_est_gamma_setting$name_kit_stutter <- data_est_sr$name_kit_stutter
          rv_est_gamma_setting$name_sr <- str_replace(basename(path_sr), ".{6}$", "")
        }
      })

      output$name_sr <- renderUI({
        disabled(textInput(session$ns("name_sr"), label = "", value = rv_est_gamma_setting$name_sr))
      })

      # --- Load data for LSAE --- #

      shinyFileChoose(input, "load_lsae", roots = c(documents = paste0(path_pack, "/extdata/lsae")), session = session, filetypes = c("RData"))

      observeEvent(input$load_lsae, {
        req(input$load_lsae)
        path_lsae <- parseFilePaths(roots = c(documents = paste0(path_pack, "/extdata/lsae")), input$load_lsae)$datapath
        if(length(path_lsae) > 0){
          load(path_lsae)
          rv_est_gamma_setting$dt_expect_ae <- data_est_lsae$dt_expect_ae
          rv_est_gamma_setting$name_kit_lsae <- data_est_lsae$name_kit_lsae
          rv_est_gamma_setting$name_lsae <- str_replace(basename(path_lsae), ".{6}$", "")
        }
      })

      output$name_lsae <- renderUI({
        disabled(textInput(session$ns("name_lsae"), label = "", value = rv_est_gamma_setting$name_lsae))
      })

      # --- Prepare dataset for the kit --- #

      observe({
        req(rv_est_gamma_setting$name_kit_stutter)
        req(rv_est_gamma_setting$name_kit_lsae)

        name_kit_stutter <- rv_est_gamma_setting$name_kit_stutter
        name_kit_lsae <- rv_est_gamma_setting$name_kit_lsae

        if(name_kit_stutter == name_kit_lsae){
          rv_est_gamma_setting$name_kit_gamma <- name_kit_stutter

          dt_kit <- fread(paste0(path_pack, "/extdata/kit/", name_kit_stutter, ".csv"))
          dt_kit <- try(change_data_type_kit(dt_kit))
          if(inherits(dt_kit, "try-error")){
            rv_est_gamma_setting$dt_setting_est_par <- NULL
            showModal(modalDialog(title = "Error", "The file format of the selected kit is invalid. Please update a file in the correct format.", easyClose = TRUE, footer = NULL))
          }else{
            rv_est_gamma_setting$dt_setting_est_par <- create_init_dt_setting_est_par(dt_kit)
          }
        }

        else{
          rv_est_gamma_setting$name_kit_gamma <- NULL

          showModal(modalDialog(title = "Error", "The kit name in dataset for the stutter ratios does not match the kit name in dataset for the locus-specific amplification efficiencies.", easyClose = TRUE, footer = NULL))
        }
      })

      output$name_kit_gamma <- renderUI({
        disabled(textInput(session$ns("name_kit_gamma"), label = "", value = rv_est_gamma_setting$name_kit_gamma))
      })

      output$dt_setting_est_par <- renderRHandsontable({
        req(rv_est_gamma_setting$dt_setting_est_par)

        rhandsontable(rv_est_gamma_setting$dt_setting_est_par, rowHeaders = FALSE) %>%
          hot_table(fillHandle = FALSE, rowHeights = 30) %>%
          hot_col(c("Marker"), type = "text", readOnly = TRUE, width = 100) %>%
          hot_col(c("AT"), type = "numeric", format = "0", width = 50)
      })

      # --- Advanced setting --- #

      observeEvent(input$act_advance_gamma, {
        showModal(modalDialog(
          title = "MLE conditions",
          fluidRow(
            column(6,
                   numericInput(session$ns("max_iter"), label = "Maximum number of iterations", value = rv_est_gamma_setting$max_iter),
                   numericInput(session$ns("max_mle_runs"), label = "Maximum number of MLE attempts", value = rv_est_gamma_setting$max_mle_runs),
                   numericInput(session$ns("num_consistent"), label = "Minimum number of consistent runs", value = rv_est_gamma_setting$num_consistent),
                   numericInput(session$ns("max_diff_like"), label = "Reproducibility threshold (log10-likelihood difference)", value = rv_est_gamma_setting$max_diff_like)
            ),
            column(6,
                   numericInput(session$ns("min_scale_allele"), label = "Minimum scale parameter for alleles", value = rv_est_gamma_setting$min_scale_allele),
                   numericInput(session$ns("max_scale_allele"), label = "Maximum scale parameter for alleles", value = rv_est_gamma_setting$max_scale_allele),
                   numericInput(session$ns("min_scale_stutter"), label = "Minimum scale parameter for stutters", value = rv_est_gamma_setting$min_scale_stutter),
                   numericInput(session$ns("max_scale_stutter"), label = "Maximum scale parameter for stutters", value = rv_est_gamma_setting$max_scale_stutter)
            )
          ),
          footer = tagList(
            actionButton(session$ns("act_advance_gamma_save"), "Save"),
            modalButton("Cancel")
          ),
          size = "l"
        ))
      })

      observeEvent(input$min_scale_allele | input$max_scale_allele, {
        min_scale_allele <- input$min_scale_allele
        max_scale_allele <- input$max_scale_allele

        hideFeedback("min_scale_allele")
        hideFeedback("max_scale_allele")
        rv_est_gamma_setting$bool_correct_mle_scale_allele <- TRUE

        if(!is.numeric(min_scale_allele) || min_scale_allele <= 0){
          showFeedbackDanger(inputId = "min_scale_allele", text = "Value must be a positive number.")
          rv_est_gamma_setting$bool_correct_mle_scale_allele <- FALSE
        }

        if(!is.numeric(max_scale_allele) || max_scale_allele <= 0){
          showFeedbackDanger(inputId = "max_scale_allele", text = "Value must be a positive number.")
          rv_est_gamma_setting$bool_correct_mle_scale_allele <- FALSE
        }

        if(is.numeric(min_scale_allele) && is.numeric(max_scale_allele) && max_scale_allele < min_scale_allele){
          showFeedbackDanger(inputId = "min_scale_allele", text = "The minimum value cannot exceed the maximum value.")
          showFeedbackDanger(inputId = "max_scale_allele", text = "The minimum value cannot exceed the maximum value.")
          rv_est_gamma_setting$bool_correct_mle_scale_allele <- FALSE
        }
      })

      observeEvent(input$min_scale_stutter | input$max_scale_stutter, {
        min_scale_stutter <- input$min_scale_stutter
        max_scale_stutter <- input$max_scale_stutter

        hideFeedback("min_scale_stutter")
        hideFeedback("max_scale_stutter")
        rv_est_gamma_setting$bool_correct_mle_scale_stutter <- TRUE

        if(!is.numeric(min_scale_stutter) || min_scale_stutter <= 0){
          showFeedbackDanger(inputId = "min_scale_stutter", text = "Value must be a positive number.")
          rv_est_gamma_setting$bool_correct_mle_scale_stutter <- FALSE
        }

        if(!is.numeric(max_scale_stutter) || max_scale_stutter <= 0){
          showFeedbackDanger(inputId = "max_scale_stutter", text = "Value must be a positive number.")
          rv_est_gamma_setting$bool_correct_mle_scale_stutter <- FALSE
        }

        if(is.numeric(min_scale_stutter) && is.numeric(max_scale_stutter) && max_scale_stutter < min_scale_stutter){
          showFeedbackDanger(inputId = "min_scale_stutter", text = "The minimum value cannot exceed the maximum value.")
          showFeedbackDanger(inputId = "max_scale_stutter", text = "The minimum value cannot exceed the maximum value.")
          rv_est_gamma_setting$bool_correct_mle_scale_stutter <- FALSE
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

        rv_est_gamma_setting$bool_correct_advanced_1 <- bool_correct_advanced_1
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

        rv_est_gamma_setting$bool_correct_advanced_2 <- bool_correct_advanced_2
      })

      observeEvent(input$max_diff_like, {
        max_diff_like <- input$max_diff_like

        hideFeedback("max_diff_like")
        bool_correct_advanced_3 <- TRUE

        if(!is.numeric(max_diff_like) || max_diff_like < 0){
          showFeedbackDanger(inputId = "max_diff_like", text = "Value must be a positive number.")
          bool_correct_advanced_3 <- FALSE
        }

        rv_est_gamma_setting$bool_correct_advanced_3 <- bool_correct_advanced_3
      })

      observe({
        req(rv_est_gamma)

        if(all(rv_est_gamma_setting$bool_correct_mle_scale_allele,
               rv_est_gamma_setting$bool_correct_mle_scale_stutter,
               rv_est_gamma_setting$bool_correct_advanced_1,
               rv_est_gamma_setting$bool_correct_advanced_2,
               rv_est_gamma_setting$bool_correct_advanced_3)){
          enable("act_advance_gamma_save")
        }else{
          disable("act_advance_gamma_save")
        }
      })

      observeEvent(input$act_advance_gamma_save, {
        rv_est_gamma_setting$min_scale_allele <- input$min_scale_allele
        rv_est_gamma_setting$max_scale_allele <- input$max_scale_allele
        rv_est_gamma_setting$min_scale_stutter <- input$min_scale_stutter
        rv_est_gamma_setting$max_scale_stutter <- input$max_scale_stutter
        rv_est_gamma_setting$max_iter <- input$max_iter
        rv_est_gamma_setting$num_consistent <- input$num_consistent
        rv_est_gamma_setting$max_mle_runs <- input$max_mle_runs
        rv_est_gamma_setting$max_diff_like <- input$max_diff_like
        removeModal()
      })

      output$ui_analyze_est_gamma <- renderUI({
        dt_exp_gamma <- rv_est_gamma_setting$dt_exp_gamma
        dt_gt_gamma <- rv_est_gamma_setting$dt_gt_gamma
        dt_setting_est_par <- rv_est_gamma_setting$dt_setting_est_par
        name_kit_gamma <- rv_est_gamma_setting$name_kit_gamma
        if(any(is.null(dt_exp_gamma), is.null(dt_gt_gamma), is.null(dt_setting_est_par), is.null(name_kit_gamma))){
          disabled(actionButton(session$ns("act_analyze_est_gamma"), label = "Analyze", class = "btn btn-primary btn-lg"))
        }else{
          actionButton(session$ns("act_analyze_est_gamma"), label = "Analyze", class = "btn btn-primary btn-lg")
        }
      })

      # Start estimating parameters
      observeEvent(input$act_analyze_est_gamma, {
        waiter_ui <- Waiter$new(html = tagList(spin_3k(),
                                               h3("Checking the settings...")),
                                color = transparent(.5))
        waiter_ui$show()

        # Fix data
        fn_exp_gamma <- rv_est_gamma_setting$fn_exp_gamma
        dt_exp_gamma <- rv_est_gamma_setting$dt_exp_gamma
        fn_gt_gamma <- rv_est_gamma_setting$fn_gt_gamma
        dt_gt_gamma <- rv_est_gamma_setting$dt_gt_gamma
        dt_setting_est_par <- hot_to_r(input$dt_setting_est_par)
        list_dt_average_sr <- rv_est_gamma_setting$list_dt_average_sr
        name_sr <- rv_est_gamma_setting$name_sr
        dt_expect_ae <- rv_est_gamma_setting$dt_expect_ae
        name_lsae <- rv_est_gamma_setting$name_lsae
        name_kit_gamma <- rv_est_gamma_setting$name_kit_gamma
        min_scale_allele <- rv_est_gamma_setting$min_scale_allele
        max_scale_allele <- rv_est_gamma_setting$max_scale_allele
        min_scale_stutter <- rv_est_gamma_setting$min_scale_stutter
        max_scale_stutter <- rv_est_gamma_setting$max_scale_stutter
        max_iter <- rv_est_gamma_setting$max_iter
        num_consistent <- rv_est_gamma_setting$num_consistent
        max_mle_runs <- rv_est_gamma_setting$max_mle_runs
        max_diff_like <- rv_est_gamma_setting$max_diff_like

        # Check input files and setting
        error_message <- check_est_par(dt_exp_gamma, dt_gt_gamma, dt_setting_est_par)
        if(error_message != ""){
          waiter_ui$hide()
          showModal(modalDialog(title = "Error", HTML(error_message), easyClose = TRUE, footer = NULL))
        }else{
          data_est_gamma <- try(estimate_gamma(dt_exp_gamma, dt_gt_gamma, list_dt_average_sr, dt_setting_est_par, dt_expect_ae,
                                               min_scale_allele, max_scale_allele, min_scale_stutter, max_scale_stutter,
                                               max_iter, num_consistent, max_mle_runs, max_diff_like,
                                               waiter_ui))
          if(inherits(data_est_gamma, "try-error")){
            waiter_ui$hide()
            showModal(modalDialog(title = "Error", "An unexpected error has occurred.", easyClose = TRUE, footer = NULL))
          }else{
            data_est_gamma$fn_exp_gamma <- fn_exp_gamma
            data_est_gamma$fn_gt_gamma <- fn_gt_gamma
            data_est_gamma$dt_setting_est_par <- dt_setting_est_par
            data_est_gamma$list_dt_average_sr <- list_dt_average_sr
            data_est_gamma$name_kit_gamma <- name_kit_gamma
            data_est_gamma$name_sr <- name_sr
            data_est_gamma$name_lsae <- name_lsae
            data_est_gamma$min_scale_allele <- min_scale_allele
            data_est_gamma$max_scale_allele <- max_scale_allele
            data_est_gamma$min_scale_stutter <- min_scale_stutter
            data_est_gamma$max_scale_stutter <- max_scale_stutter
            data_est_gamma$max_iter <- max_iter
            data_est_gamma$num_consistent <- num_consistent
            data_est_gamma$max_mle_runs <- max_mle_runs
            data_est_gamma$max_diff_like <- max_diff_like

            rv_est_gamma$data_est_gamma <- data_est_gamma

            waiter_ui$hide()

            updateTabsetPanel(session = session_gamma_top, inputId = "tabset_est_gamma", selected = "Result")
          }
        }
      })

      return(rv_est_gamma)
    }
  )
}
