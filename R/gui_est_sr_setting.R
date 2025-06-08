#' est_sr_setting_ui
#'
#' @description The function to create the ui module to create the setting tab for estimating average stutter ratios
est_sr_setting_ui <- function(id){
  ns <- NS(id)

  tabPanel("Setting",
           br(),
           sidebarLayout(
             sidebarPanel(width = 3,
                          fileInput(ns("exp_data"), label = h4("Experimental data"), accept = ".csv", buttonLabel = "Load"),
                          fileInput(ns("gt_data"), label = h4("Genotype data"), accept = ".csv", buttonLabel = "Load"),
                          h4("Kit"),
                          fluidRow(
                            column(3, br(), shinyFilesButton(ns("load_kit"), label = "Load", title = "Kit", multiple = FALSE)),
                            column(9, uiOutput(ns("name_kit_stutter")))
                          ),
                          br(),
                          actionButton(ns("act_advance_sr"), "Advanced setting"),
                          br(),
                          br(),
                          br(),
                          uiOutput(ns("act_analyze_est_sr"))
             ),
             mainPanel(width = 9,
                       tags$style(HTML("
                         .handsontable .htCore td {
                           vertical-align: middle !important;
                           text-align: center !important;
                         }"
                       )),
                       rHandsontableOutput(ns("dt_setting_est_sr")),
                       br()
             )
           )
  )
}

#' create_init_dt_setting_est_sr
#'
#' @description The function to create the initial data.table in the setting tab for estimating average stutter ratios
#' @param dt_kit The data.table of information on a kit
create_init_dt_setting_est_sr <- function(dt_kit){
  dt_kit_auto <- extract_kit_auto(dt_kit)
  markers <- unique(dt_kit_auto$Marker)

  dt_setting_est_sr <- data.table(Marker = markers)

  if(nrow(dt_setting_est_sr) == 21){
    dt_setting_est_sr$AT <- c(rep(35, 5), rep(40, 3), rep(25, 4), rep(25, 5), rep(45, 4))
    dt_setting_est_sr$Model_BSR <- c("Linear (Allele)", "Linear (Allele)", "Linear (Allele)", "Linear (Allele)", "Linear (Allele)",
                                     "Linear (Corrected Allele)", "Linear (Corrected Allele)", "Linear (Allele)",
                                     "Linear (Corrected Allele)", "Linear (Allele)", "Linear (Corrected Allele)", "Linear (Allele)",
                                     "Linear (Allele)", "Linear (Allele)", "Linear (Allele)", "Linear (Allele)", "Linear (Corrected Allele)",
                                     "Linear (Allele)", "Linear (Corrected Allele)", "Linear (Corrected Allele)", "Linear (Corrected Allele)")
    dt_setting_est_sr$Model_FSR <- "Uniform (Multiple loci)"
    dt_setting_est_sr$Model_FSR[13] <- "Linear (Allele)"
    dt_setting_est_sr$Model_DSR <- "Uniform (Multiple loci)"
    dt_setting_est_sr$Model_M2SR <- "Not consider"
    dt_setting_est_sr$Model_M2SR[17] <- "Uniform"
    dt_setting_est_sr$Model_M2SR[19] <- "Uniform"
  }else{
    dt_setting_est_sr$AT <- 50
    dt_setting_est_sr$Model_BSR <- "Linear (Allele)"
    dt_setting_est_sr$Model_FSR <- "Uniform (Multiple loci)"
    dt_setting_est_sr$Model_DSR <- "Uniform (Multiple loci)"
    dt_setting_est_sr$Model_M2SR <- "Not consider"
  }

  col_int <- c("AT")
  dt_setting_est_sr[, (col_int) := lapply(.SD, as.integer), .SDcols = col_int]

  names(dt_setting_est_sr) <- c("Marker", "AT", "BSR model", "FSR model", "DSR model", "M2SR model")

  return(dt_setting_est_sr)
}

#' est_sr_setting_server
#'
#' @description The function to create the server module to create the setting tab for estimating average stutter ratios
est_sr_setting_server <- function(id, path_pack, session_est_sr_top){
  moduleServer(
    id,
    function(input, output, session){
      rv_est_sr <- reactiveValues()
      rv_est_sr$dt_exp_sr <- NULL
      rv_est_sr$fn_exp_sr <- NULL
      rv_est_sr$dt_gt_sr <- NULL
      rv_est_sr$fn_gt_sr <- NULL
      rv_est_sr$dt_kit <- NULL
      rv_est_sr$dt_setting_est_sr <- NULL
      rv_est_sr$name_kit_stutter <- NULL
      rv_est_sr$min_slope <- 0
      rv_est_sr$max_slope <- 1
      rv_est_sr$min_intercept <- -1
      rv_est_sr$max_intercept <- 1
      rv_est_sr$min_myu <- 0.001
      rv_est_sr$max_myu <- 1
      rv_est_sr$min_sigma <- 1
      rv_est_sr$max_sigma <- 10000
      rv_est_sr$max_call <- 3000
      rv_est_sr$num_converged_solutions <- 3
      rv_est_sr$max_mle_runs <- 10
      rv_est_sr$bool_correct_mle_slope <- TRUE
      rv_est_sr$bool_correct_mle_intercept <- TRUE
      rv_est_sr$bool_correct_mle_myu <- TRUE
      rv_est_sr$bool_correct_mle_sigma <- TRUE
      rv_est_sr$bool_correct_advanced_1 <- TRUE
      rv_est_sr$bool_correct_advanced_2 <- TRUE
      rv_est_sr$data_est_sr <- NULL

      # Load experimental data
      observeEvent(input$exp_data, {
        file_input <- input$exp_data
        if(!is.null(file_input)){
          dt_exp_sr <- fread(file_input$datapath)
          dt_exp_sr <- change_data_type_exp(dt_exp_sr)
          rv_est_sr$dt_exp_sr <- dt_exp_sr
          rv_est_sr$fn_exp_sr <- file_input$name
        }
      })

      # Load genotype data
      observeEvent(input$gt_data, {
        file_input <- input$gt_data
        if(!is.null(file_input)){
          dt_gt_sr <- fread(file_input$datapath)
          dt_gt_sr <- change_data_type_gt(dt_gt_sr)
          rv_est_sr$dt_gt_sr <- dt_gt_sr
          rv_est_sr$fn_gt_sr <- file_input$name
        }
      })

      # Select a kit
      shinyFileChoose(input, "load_kit", roots = c(documents = paste0(path_pack, "/extdata/kit")), session = session, filetypes = c("csv"))

      observeEvent(input$load_kit, {
        req(input$load_kit)
        path_kit <- parseFilePaths(roots = c(documents = paste0(path_pack, "/extdata/kit")), input$load_kit)$datapath
        if(length(path_kit) > 0){
          dt_kit <- fread(path_kit)
          rv_est_sr$dt_kit <- dt_kit <- change_data_type_kit(dt_kit)
          rv_est_sr$dt_setting_est_sr <- create_init_dt_setting_est_sr(dt_kit)
          rv_est_sr$name_kit_stutter <- str_replace(basename(path_kit), ".{4}$", "")
        }
      })

      output$name_kit_stutter <- renderUI({
        disabled(textInput(session$ns("name_kit_stutter"), label = "", value = rv_est_sr$name_kit_stutter))
      })

      observeEvent(input$act_advance_sr, {
        showModal(modalDialog(
          title = "Advanced setting",
          h4("MLE conditions"),
          fluidRow(
            column(6,
                   numericInput(session$ns("min_slope"), label = "Minimum slope", value = rv_est_sr$min_slope),
                   numericInput(session$ns("min_intercept"), label = "Minimum intercept", value = rv_est_sr$min_intercept),
                   numericInput(session$ns("min_myu"), label = "Minimum myu", value = rv_est_sr$min_myu),
                   numericInput(session$ns("min_sigma"), label = "Minimum sigma", value = rv_est_sr$min_sigma),
                   numericInput(session$ns("max_call"), label = "Number of simulations", value = rv_est_sr$max_call),
                   numericInput(session$ns("max_mle_runs"), label = "The maximum number of MLE attempts with different initial values", value = rv_est_sr$max_mle_runs)
            ),
            column(6,
                   numericInput(session$ns("max_slope"), label = "Maximum slope", value = rv_est_sr$max_slope),
                   numericInput(session$ns("max_intercept"), label = "Maximum intercept", value = rv_est_sr$max_intercept),
                   numericInput(session$ns("max_myu"), label = "Maximum myu", value = rv_est_sr$max_myu),
                   numericInput(session$ns("max_sigma"), label = "Maximum sigma", value = rv_est_sr$max_sigma),
                   numericInput(session$ns("num_converged_solutions"), label = "The number of converged solutions required to terminate the MLE", value = rv_est_sr$num_converged_solutions)
            )
          ),
          footer = tagList(
            actionButton(session$ns("act_advance_sr_save"), "Save"),
            modalButton("Cancel")
          ),
          size = "l"
        ))
      })

      observeEvent(input$min_slope | input$max_slope, {
        min_slope <- input$min_slope
        max_slope <- input$max_slope

        hideFeedback("min_slope")
        hideFeedback("max_slope")
        rv_est_sr$bool_correct_mle_slope <- TRUE

        if(!is.numeric(min_slope) || min_slope < 0){
          showFeedbackDanger(inputId = "min_slope", text = "Value must be non-negative.")
          rv_est_sr$bool_correct_mle_slope <- FALSE
        }

        if(!is.numeric(max_slope) || max_slope < 0){
          showFeedbackDanger(inputId = "max_slope", text = "Value must be non-negative.")
          rv_est_sr$bool_correct_mle_slope <- FALSE
        }

        if(is.numeric(min_slope) && is.numeric(max_slope) && max_slope < min_slope){
          showFeedbackDanger(inputId = "min_slope", text = "Minimum value cannot be greater than maximum value.")
          showFeedbackDanger(inputId = "max_slope", text = "Minimum value cannot be greater than maximum value.")
          rv_est_sr$bool_correct_mle_slope <- FALSE
        }
      })

      observeEvent(input$min_intercept | input$max_intercept, {
        min_intercept <- input$min_intercept
        max_intercept <- input$max_intercept

        hideFeedback("min_intercept")
        hideFeedback("max_intercept")
        rv_est_sr$bool_correct_mle_intercept <- TRUE

        if(!is.numeric(min_intercept)){
          showFeedbackDanger(inputId = "min_intercept", text = "Only numeric values are allowed.")
          rv_est_sr$bool_correct_mle_intercept <- FALSE
        }

        if(!is.numeric(max_intercept)){
          showFeedbackDanger(inputId = "max_intercept", text = "Only numeric values are allowed.")
          rv_est_sr$bool_correct_mle_intercept <- FALSE
        }

        if(is.numeric(min_intercept) && is.numeric(max_intercept) && max_intercept < min_intercept){
          showFeedbackDanger(inputId = "min_intercept", text = "Minimum value cannot be greater than maximum value.")
          showFeedbackDanger(inputId = "max_intercept", text = "Minimum value cannot be greater than maximum value.")
          rv_est_sr$bool_correct_mle_intercept <- FALSE
        }
      })

      observeEvent(input$min_myu | input$max_myu, {
        min_myu <- input$min_myu
        max_myu <- input$max_myu

        hideFeedback("min_myu")
        hideFeedback("max_myu")
        rv_est_sr$bool_correct_mle_myu <- TRUE

        if(!is.numeric(min_myu) || min_myu <= 0){
          showFeedbackDanger(inputId = "min_myu", text = "Value must be positive.")
          rv_est_sr$bool_correct_mle_myu <- FALSE
        }

        if(!is.numeric(max_myu) || max_myu <= 0){
          showFeedbackDanger(inputId = "max_myu", text = "Value must be positive.")
          rv_est_sr$bool_correct_mle_myu <- FALSE
        }

        if(is.numeric(min_myu) && is.numeric(max_myu) && max_myu < min_myu){
          showFeedbackDanger(inputId = "min_myu", text = "Minimum value cannot be greater than maximum value.")
          showFeedbackDanger(inputId = "max_myu", text = "Minimum value cannot be greater than maximum value.")
          rv_est_sr$bool_correct_mle_myu <- FALSE
        }
      })

      observeEvent(input$min_sigma | input$max_sigma, {
        min_sigma <- input$min_sigma
        max_sigma <- input$max_sigma

        hideFeedback("min_sigma")
        hideFeedback("max_sigma")
        rv_est_sr$bool_correct_mle_sigma <- TRUE

        if(!is.numeric(min_sigma) || min_sigma <= 0){
          showFeedbackDanger(inputId = "min_sigma", text = "Value must be positive.")
          rv_est_sr$bool_correct_mle_sigma <- FALSE
        }

        if(!is.numeric(max_sigma) || max_sigma <= 0){
          showFeedbackDanger(inputId = "max_sigma", text = "Value must be positive.")
          rv_est_sr$bool_correct_mle_sigma <- FALSE
        }

        if(is.numeric(min_sigma) && is.numeric(max_sigma) && max_sigma < min_sigma){
          showFeedbackDanger(inputId = "min_sigma", text = "Minimum value cannot be greater than maximum value.")
          showFeedbackDanger(inputId = "max_sigma", text = "Minimum value cannot be greater than maximum value.")
          rv_est_sr$bool_correct_mle_sigma <- FALSE
        }
      })

      observeEvent(input$max_call, {
        max_call <- input$max_call

        hideFeedback("max_call")
        bool_correct_advanced_1 <- TRUE

        if(!is.integer(max_call) || max_call <= 0){
          showFeedbackDanger(inputId = "max_call", text = "Value must be an integer greater than zero.")
          bool_correct_advanced_1 <- FALSE
        }

        rv_est_sr$bool_correct_advanced_1 <- bool_correct_advanced_1
      })

      observeEvent(input$num_converged_solutions | input$max_mle_runs, {
        num_converged_solutions <- input$num_converged_solutions
        max_mle_runs <- input$max_mle_runs

        hideFeedback("num_converged_solutions")
        hideFeedback("max_mle_runs")
        bool_correct_advanced_2 <- TRUE

        if(!is.integer(num_converged_solutions) || num_converged_solutions <= 0){
          showFeedbackDanger(inputId = "num_converged_solutions", text = "Value must be an integer greater than zero.")
          bool_correct_advanced_2 <- FALSE
        }

        if(!is.integer(max_mle_runs) || max_mle_runs <= 0){
          showFeedbackDanger(inputId = "max_mle_runs", text = "Value must be an integer greater than zero.")
          bool_correct_advanced_2 <- FALSE
        }

        if(is.integer(num_converged_solutions) && is.integer(max_mle_runs) && num_converged_solutions > max_mle_runs){
          showFeedbackDanger(inputId = "num_converged_solutions", text = "The number of converged solutions must not exceed the maximum number of MLE attempts.")
          showFeedbackDanger(inputId = "max_mle_runs", text = "The number of converged solutions must not exceed the maximum number of MLE attempts.")
          bool_correct_advanced_2 <- FALSE
        }

        rv_est_sr$bool_correct_advanced_2 <- bool_correct_advanced_2
      })

      observe({
        req(rv_est_sr)

        if(all(rv_est_sr$bool_correct_mle_slope,
               rv_est_sr$bool_correct_mle_intercept,
               rv_est_sr$bool_correct_mle_myu,
               rv_est_sr$bool_correct_mle_sigma,
               rv_est_sr$bool_correct_advanced_1,
               rv_est_sr$bool_correct_advanced_2)){
          enable("act_advance_sr_save")
        }else{
          disable("act_advance_sr_save")
        }
      })

      observeEvent(input$act_advance_sr_save, {
        rv_est_sr$min_slope <- input$min_slope
        rv_est_sr$max_slope <- input$max_slope
        rv_est_sr$min_intercept <- input$min_intercept
        rv_est_sr$max_intercept <- input$max_intercept
        rv_est_sr$min_myu <- input$min_myu
        rv_est_sr$max_myu <- input$max_myu
        rv_est_sr$min_sigma <- input$min_sigma
        rv_est_sr$max_sigma <- input$max_sigma
        rv_est_sr$max_call <- input$max_call
        rv_est_sr$num_converged_solutions <- input$num_converged_solutions
        rv_est_sr$max_mle_runs <- input$max_mle_runs
        removeModal()
      })

      output$act_analyze_est_sr <- renderUI({
        dt_exp_sr <- rv_est_sr$dt_exp_sr
        dt_gt_sr <- rv_est_sr$dt_gt_sr
        dt_kit <- rv_est_sr$dt_kit
        if(any(is.null(dt_exp_sr), is.null(dt_gt_sr), is.null(dt_kit))){
          disabled(actionButton(session$ns("act_analyze_est_sr"), label = "Analyze", class = "btn btn-primary btn-lg"))
        }else{
          actionButton(session$ns("act_analyze_est_sr"), label = "Analyze", class = "btn btn-primary btn-lg")
        }
      })

      # Start estimating average stutter ratios
      observeEvent(input$act_analyze_est_sr, {
        dt_exp_sr <- rv_est_sr$dt_exp_sr
        dt_gt_sr <- rv_est_sr$dt_gt_sr
        dt_kit <- rv_est_sr$dt_kit
        name_kit_stutter <- rv_est_sr$name_kit_stutter
        rv_est_sr$dt_setting_est_sr <- dt_setting_est_sr <- hot_to_r(input$dt_setting_est_sr)
        min_slope <- rv_est_sr$min_slope
        max_slope <- rv_est_sr$max_slope
        min_intercept <- rv_est_sr$min_intercept
        max_intercept <- rv_est_sr$max_intercept
        min_myu <- rv_est_sr$min_myu
        max_myu <- rv_est_sr$max_myu
        min_sigma <- rv_est_sr$min_sigma
        max_sigma <- rv_est_sr$max_sigma
        max_call <- rv_est_sr$max_call
        num_converged_solutions <- rv_est_sr$num_converged_solutions
        max_mle_runs <- rv_est_sr$max_mle_runs

        error_message <- check_est_sr(dt_exp_sr, dt_gt_sr, dt_kit, dt_setting_est_sr)
        if(error_message != ""){
          showModal(modalDialog(title = "Error", HTML(error_message), easyClose = TRUE, footer = NULL))
        }else{
          rv_est_sr$data_est_sr <- estimate_average_sr(dt_exp_sr, dt_gt_sr, dt_kit, dt_setting_est_sr,
                                                       min_slope, max_slope, min_intercept, max_intercept, min_myu, max_myu, min_sigma, max_sigma, max_call, num_converged_solutions, max_mle_runs)

          updateTabsetPanel(session = session_est_sr_top, inputId = "tabset_est_sr", selected = "Result")
        }
      })

      output$dt_setting_est_sr <- renderRHandsontable({
        dt_setting_est_sr <- rv_est_sr$dt_setting_est_sr
        if(!is.null(dt_setting_est_sr)){
          rhandsontable(rv_est_sr$dt_setting_est_sr, rowHeaders = FALSE, overflow = "visible") %>%
            hot_table(fillHandle = FALSE, stretchH = "all", rowHeights = 30) %>%
            hot_col(c("Marker"), type = "text", readOnly = TRUE) %>%
            hot_col(c("AT"), type = "numeric", format = "0") %>%
            hot_col(c("BSR model", "FSR model", "DSR model", "M2SR model"), type = "dropdown", source = c("Linear (Allele)", "Linear (Corrected Allele)", "Uniform", "Uniform (Multiple loci)", "Not consider"))
        }
      })

      return(rv_est_sr)
    }
  )
}
