#' est_par_setting_ui
#'
#' @description The function to create the ui module for setting to estimate parameters
est_par_setting_ui <- function(id){
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
                          uiOutput(ns("name_kit_param")),
                          br(),
                          actionButton(ns("act_advance_par"), "Advanced setting"),
                          br(),
                          br(),
                          br(),
                          actionButton(ns("act_analyze_est_par"), label = "Analyze", class = "btn btn-primary btn-lg")
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

#' est_par_setting_server
#'
#' @description The function to create the server module for setting to estimate parameters
est_par_setting_server <- function(id, path_pack, session_param_top){
  moduleServer(
    id,
    function(input, output, session){
      rv_est_par_setting <- reactiveValues()
      rv_est_par_setting$dt_exp_param <- NULL
      rv_est_par_setting$fn_exp_param <- NULL
      rv_est_par_setting$dt_gt_param <- NULL
      rv_est_par_setting$fn_gt_param <- NULL
      rv_est_par_setting$dt_setting_est_par <- NULL
      rv_est_par_setting$list_dt_average_sr <- NULL
      rv_est_par_setting$name_kit_param <- NULL
      rv_est_par_setting$name_sr <- NULL
      rv_est_par_setting$min_alpha <- 0.01
      rv_est_par_setting$max_alpha <- 100
      rv_est_par_setting$min_scale_allele <- 0.01
      rv_est_par_setting$max_scale_allele <- 400
      rv_est_par_setting$min_scale_stutter <- 0.01
      rv_est_par_setting$max_scale_stutter <- 400
      rv_est_par_setting$stop_iter <- 500
      rv_est_par_setting$max_iter <- 5000
      rv_est_par_setting$num_converged <- 3
      rv_est_par_setting$max_mle_runs <- 10
      rv_est_par_setting$bool_correct_mle_alpha <- TRUE
      rv_est_par_setting$bool_correct_mle_scale_allele <- TRUE
      rv_est_par_setting$bool_correct_mle_scale_stutter <- TRUE
      rv_est_par_setting$bool_correct_advanced_1 <- TRUE
      rv_est_par_setting$bool_correct_advanced_2 <- TRUE

      rv_est_par <- reactiveValues()
      rv_est_par$data_est_par <- NULL

      # Load experimental data
      observeEvent(input$exp_data, {
        file_input <- input$exp_data
        if(!is.null(file_input)){
          dt_exp_param <- fread(file_input$datapath)
          dt_exp_param <- change_data_type_exp(dt_exp_param)
          rv_est_par_setting$dt_exp_param <- dt_exp_param
          rv_est_par_setting$fn_exp_param <- file_input$name
        }
      })

      # Load genotype data
      observeEvent(input$gt_data, {
        file_input <- input$gt_data
        if(!is.null(file_input)){
          dt_gt_param <- fread(file_input$datapath)
          dt_gt_param <- change_data_type_gt(dt_gt_param)
          rv_est_par_setting$dt_gt_param <- dt_gt_param
          rv_est_par_setting$fn_gt_param <- file_input$name
        }
      })

      # Select a file for stutter ratios
      shinyFileChoose(input, "load_sr", roots = c(documents = paste0(path_pack, "/extdata/stutter")), session = session, filetypes = c("RData"))

      observeEvent(input$load_sr, {
        req(input$load_sr)
        path_sr <- parseFilePaths(roots = c(documents = paste0(path_pack, "/extdata/stutter")), input$load_sr)$datapath
        if(length(path_sr) > 0){
          load(path_sr)
          rv_est_par_setting$list_dt_average_sr <- data_est_sr$list_dt_average_sr
          rv_est_par_setting$name_kit_param <- data_est_sr$name_kit_stutter
          rv_est_par_setting$name_sr <- str_replace(basename(path_sr), ".{6}$", "")
        }
      })

      output$name_sr <- renderUI({
        disabled(textInput(session$ns("name_sr"), label = "", value = rv_est_par_setting$name_sr))
      })

      # Kit
      output$name_kit_param <- renderUI({
        disabled(textInput(session$ns("name_kit_param"), label = "", value = rv_est_par_setting$name_kit_param))
      })

      observeEvent(input$name_kit_param, ignoreInit = TRUE, {
        name_kit_param <- rv_est_par_setting$name_kit_param
        if(!is.null(name_kit_param)){
          dt_kit <- fread(paste0(path_pack, "/extdata/kit/", name_kit_param, ".csv"))
          dt_kit <- change_data_type_kit(dt_kit)
          rv_est_par_setting$dt_setting_est_par <- dt_setting_est_par <- create_init_dt_setting_est_par(dt_kit)

          output$dt_setting_est_par <- renderRHandsontable({
            rhandsontable(dt_setting_est_par, rowHeaders = FALSE) %>%
              hot_table(fillHandle = FALSE, rowHeights = 30) %>%
              hot_col(c("Marker"), type = "text", readOnly = TRUE, width = 100) %>%
              hot_col(c("AT"), type = "numeric", format = "0", width = 50)
          })
        }
      })

      #### Advanced setting ####

      observeEvent(input$act_advance_par, {
        showModal(modalDialog(
          title = "MLE conditions",
          fluidRow(
            column(6,
                   numericInput(session$ns("stop_iter"), label = "Minimum number of iterations", value = rv_est_par_setting$stop_iter),
                   numericInput(session$ns("num_converged"), label = "Minimum number of MLE attempts", value = rv_est_par_setting$num_converged),
                   numericInput(session$ns("min_alpha"), label = "Minimum alpha for LSAE", value = rv_est_par_setting$min_alpha),
                   numericInput(session$ns("min_scale_allele"), label = "Minimum scale parameter for alleles", value = rv_est_par_setting$min_scale_allele),
                   numericInput(session$ns("min_scale_stutter"), label = "Minimum scale parameter for stutters", value = rv_est_par_setting$min_scale_stutter)
            ),
            column(6,
                   numericInput(session$ns("max_iter"), label = "Maximum number of iterations", value = rv_est_par_setting$max_iter),
                   numericInput(session$ns("max_mle_runs"), label = "Maximum number of MLE attempts", value = rv_est_par_setting$max_mle_runs),
                   numericInput(session$ns("max_alpha"), label = "Maximum alpha for LSAE", value = rv_est_par_setting$max_alpha),
                   numericInput(session$ns("max_scale_allele"), label = "Maximum scale parameter for alleles", value = rv_est_par_setting$max_scale_allele),
                   numericInput(session$ns("max_scale_stutter"), label = "Maximum scale parameter for stutters", value = rv_est_par_setting$max_scale_stutter)
            )
          ),
          footer = tagList(
            actionButton(session$ns("act_advance_par_save"), "Save"),
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
        rv_est_par_setting$bool_correct_mle_alpha <- TRUE

        if(!is.numeric(min_alpha) || min_alpha <= 0){
          showFeedbackDanger(inputId = "min_alpha", text = "Value must be positive.")
          rv_est_par_setting$bool_correct_mle_alpha <- FALSE
        }

        if(!is.numeric(max_alpha) || max_alpha <= 0){
          showFeedbackDanger(inputId = "max_alpha", text = "Value must be positive.")
          rv_est_par_setting$bool_correct_mle_alpha <- FALSE
        }

        if(is.numeric(min_alpha) && is.numeric(max_alpha) && max_alpha < min_alpha){
          showFeedbackDanger(inputId = "min_alpha", text = "Minimum value cannot be greater than maximum value.")
          showFeedbackDanger(inputId = "max_alpha", text = "Minimum value cannot be greater than maximum value.")
          rv_est_par_setting$bool_correct_mle_alpha <- FALSE
        }
      })

      observeEvent(input$min_scale_allele | input$max_scale_allele, {
        min_scale_allele <- input$min_scale_allele
        max_scale_allele <- input$max_scale_allele

        hideFeedback("min_scale_allele")
        hideFeedback("max_scale_allele")
        rv_est_par_setting$bool_correct_mle_scale_allele <- TRUE

        if(!is.numeric(min_scale_allele) || min_scale_allele <= 0){
          showFeedbackDanger(inputId = "min_scale_allele", text = "Value must be positive.")
          rv_est_par_setting$bool_correct_mle_scale_allele <- FALSE
        }

        if(!is.numeric(max_scale_allele) || max_scale_allele <= 0){
          showFeedbackDanger(inputId = "max_scale_allele", text = "Value must be positive.")
          rv_est_par_setting$bool_correct_mle_scale_allele <- FALSE
        }

        if(is.numeric(min_scale_allele) && is.numeric(max_scale_allele) && max_scale_allele < min_scale_allele){
          showFeedbackDanger(inputId = "min_scale_allele", text = "Minimum value cannot be greater than maximum value.")
          showFeedbackDanger(inputId = "max_scale_allele", text = "Minimum value cannot be greater than maximum value.")
          rv_est_par_setting$bool_correct_mle_scale_allele <- FALSE
        }
      })

      observeEvent(input$min_scale_stutter | input$max_scale_stutter, {
        min_scale_stutter <- input$min_scale_stutter
        max_scale_stutter <- input$max_scale_stutter

        hideFeedback("min_scale_stutter")
        hideFeedback("max_scale_stutter")
        rv_est_par_setting$bool_correct_mle_scale_stutter <- TRUE

        if(!is.numeric(min_scale_stutter) || min_scale_stutter <= 0){
          showFeedbackDanger(inputId = "min_scale_stutter", text = "Value must be positive.")
          rv_est_par_setting$bool_correct_mle_scale_stutter <- FALSE
        }

        if(!is.numeric(max_scale_stutter) || max_scale_stutter <= 0){
          showFeedbackDanger(inputId = "max_scale_stutter", text = "Value must be positive.")
          rv_est_par_setting$bool_correct_mle_scale_stutter <- FALSE
        }

        if(is.numeric(min_scale_stutter) && is.numeric(max_scale_stutter) && max_scale_stutter < min_scale_stutter){
          showFeedbackDanger(inputId = "min_scale_stutter", text = "Minimum value cannot be greater than maximum value.")
          showFeedbackDanger(inputId = "max_scale_stutter", text = "Minimum value cannot be greater than maximum value.")
          rv_est_par_setting$bool_correct_mle_scale_stutter <- FALSE
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

        rv_est_par_setting$bool_correct_advanced_1 <- bool_correct_advanced_1
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

        rv_est_par_setting$bool_correct_advanced_2 <- bool_correct_advanced_2
      })

      observe({
        req(rv_est_par)

        if(all(rv_est_par_setting$bool_correct_mle_alpha,
               rv_est_par_setting$bool_correct_mle_scale_allele,
               rv_est_par_setting$bool_correct_mle_scale_stutter,
               rv_est_par_setting$bool_correct_advanced_1,
               rv_est_par_setting$bool_correct_advanced_2)){
          enable("act_advance_par_save")
        }else{
          disable("act_advance_par_save")
        }
      })

      observeEvent(input$act_advance_par_save, {
        rv_est_par_setting$min_alpha <- input$min_alpha
        rv_est_par_setting$max_alpha <- input$max_alpha
        rv_est_par_setting$min_scale_allele <- input$min_scale_allele
        rv_est_par_setting$max_scale_allele <- input$max_scale_allele
        rv_est_par_setting$min_scale_stutter <- input$min_scale_stutter
        rv_est_par_setting$max_scale_stutter <- input$max_scale_stutter
        rv_est_par_setting$stop_iter <- input$stop_iter
        rv_est_par_setting$max_iter <- input$max_iter
        rv_est_par_setting$num_converged <- input$num_converged
        rv_est_par_setting$max_mle_runs <- input$max_mle_runs
        removeModal()
      })

      # Start estimating parameters
      observeEvent(input$act_analyze_est_par, {
        waiter_show(html = spin_3k(), color = "white")

        # Fix data
        fn_exp_param <- rv_est_par_setting$fn_exp_param
        dt_exp_param <- rv_est_par_setting$dt_exp_param
        fn_gt_param <- rv_est_par_setting$fn_gt_param
        dt_gt_param <- rv_est_par_setting$dt_gt_param
        dt_setting_est_par <- hot_to_r(input$dt_setting_est_par)
        list_dt_average_sr <- rv_est_par_setting$list_dt_average_sr
        name_kit_param <- rv_est_par_setting$name_kit_param
        name_sr <- rv_est_par_setting$name_sr
        min_alpha <- rv_est_par_setting$min_alpha
        max_alpha <- rv_est_par_setting$max_alpha
        min_scale_allele <- rv_est_par_setting$min_scale_allele
        max_scale_allele <- rv_est_par_setting$max_scale_allele
        min_scale_stutter <- rv_est_par_setting$min_scale_stutter
        max_scale_stutter <- rv_est_par_setting$max_scale_stutter
        stop_iter <- rv_est_par_setting$stop_iter
        max_iter <- rv_est_par_setting$max_iter
        num_converged <- rv_est_par_setting$num_converged
        max_mle_runs <- rv_est_par_setting$max_mle_runs

        # Check input files and setting
        error_message <- check_est_par(dt_exp_param, dt_gt_param, dt_setting_est_par)
        if(error_message != ""){
          waiter_hide()
          showModal(modalDialog(title = "Error", HTML(error_message), easyClose = TRUE, footer = NULL))
        }else{
          data_est_par <- estimate_parameters(dt_exp_param, dt_gt_param, list_dt_average_sr, dt_setting_est_par,
                                              min_alpha, max_alpha,
                                              min_scale_allele, max_scale_allele, min_scale_stutter, max_scale_stutter,
                                              stop_iter, max_iter, num_converged, max_mle_runs)

          data_est_par$fn_exp_param <- fn_exp_param
          data_est_par$fn_gt_param <- fn_gt_param
          data_est_par$dt_setting_est_par <- dt_setting_est_par
          data_est_par$list_dt_average_sr <- list_dt_average_sr
          data_est_par$name_kit_param <- name_kit_param
          data_est_par$name_sr <- name_sr
          data_est_par$min_alpha <- min_alpha
          data_est_par$max_alpha <- max_alpha
          data_est_par$min_scale_allele <- min_scale_allele
          data_est_par$max_scale_allele <- max_scale_allele
          data_est_par$min_scale_stutter <- min_scale_stutter
          data_est_par$max_scale_stutter <- max_scale_stutter
          data_est_par$stop_iter <- stop_iter
          data_est_par$max_iter <- max_iter
          data_est_par$num_converged <- num_converged
          data_est_par$max_mle_runs <- max_mle_runs

          rv_est_par$data_est_par <- data_est_par

          waiter_hide()

          updateTabsetPanel(session = session_param_top, inputId = "tabset_est_par", selected = "Result (LSAE)")
        }
      })

      return(rv_est_par)
    }
  )
}
