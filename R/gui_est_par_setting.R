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
                          uiOutput(ns("stutter_file")),
                          uiOutput(ns("kit_par")),
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

#' create_init_dt_setting_est_par
#'
#' @description The function to create the initial data.table in the setting tab for estimating average stutter ratios
#' @param dt_kit The data.table of information on a kit
create_init_dt_setting_est_par <- function(dt_kit){
  dt_kit_auto <- extract_kit_auto(dt_kit)
  markers <- unique(dt_kit_auto$Marker)

  dt_setting_est_par <- data.table(Marker = markers,
                                   AT = c(rep(35, 5), rep(40, 3), rep(25, 4), rep(25, 5), rep(45, 4)))
  col_int <- c("AT")
  dt_setting_est_par[, (col_int) := lapply(.SD, as.integer), .SDcols = col_int]

  return(dt_setting_est_par)
}

#' est_par_setting_server
#'
#' @description The function to create the server module for setting to estimate parameters
est_par_setting_server <- function(id, path_pack, session_param_top){
  moduleServer(
    id,
    function(input, output, session){
      rv_est_par <- reactiveValues()
      rv_est_par$dt_exp <- NULL
      rv_est_par$dt_gt <- NULL
      rv_est_par$kit_name <- NULL
      rv_est_par$warmup <- 2000
      rv_est_par$sampling <- 2000
      rv_est_par$chains <- 4
      rv_est_par$data_est_par <- NULL
      rv_est_par$list_dt_average_sr <- NULL
      rv_est_par$fn_stutter_file <- NULL

      # Load experimental data
      observeEvent(input$exp_data, {
        file_input <- input$exp_data
        if(!is.null(file_input)){
          dt_exp <- fread(file_input$datapath)
          dt_exp <- change_data_type_exp(dt_exp)
          rv_est_par$dt_exp <- dt_exp
        }
      })

      # Load genotype data
      observeEvent(input$gt_data, {
        file_input <- input$gt_data
        if(!is.null(file_input)){
          dt_gt <- fread(file_input$datapath)
          dt_gt <- change_data_type_gt(dt_gt)
          rv_est_par$dt_gt <- dt_gt
        }
      })

      # Select a stutter file
      output$stutter_file <- renderUI({
        ns <- session$ns
        stutter_files <- list.files(paste0(path_pack, "/extdata/stutter_file"))
        selectInput(ns("stutter_file"), label = h4("Stutter file"), choices = stutter_files, selected = stutter_files[1], multiple = FALSE, selectize = FALSE)
      })

      # Select a kit
      output$kit_par <- renderUI({
        kit_names <- list.files(paste0(path_pack, "/extdata/kit"))
        kit_names <- gsub(".csv", "", kit_names)
        selectInput(session$ns("kit_name"), label = h4("Kit"), choices = kit_names, selected = kit_names[1], multiple = FALSE, selectize = FALSE)
      })

      observeEvent(input$kit_name, ignoreInit = TRUE, {
        rv_est_par$kit_name <- kit_name <- input$kit_name
        dt_kit <- fread(paste0(path_pack, "/extdata/kit/", kit_name, ".csv"))
        dt_kit <- change_data_type_kit(dt_kit)
        dt_setting_est_par <- create_init_dt_setting_est_par(dt_kit)

        output$dt_setting_est_par <- renderRHandsontable({
          rhandsontable(dt_setting_est_par, rowHeaders = FALSE) %>%
            hot_table(fillHandle = FALSE, rowHeights = 30) %>%
            hot_col(c("Marker"), type = "text", readOnly = TRUE, width = 100) %>%
            hot_col(c("AT"), type = "numeric", format = "0", width = 50)
        })
      })

      #### Advanced setting ####

      observeEvent(input$act_advance_par, {
        showModal(modalDialog(
          title = "Advanced setting",
          h4("MCMC conditions"),
          numericInput(session$ns("warmup"), label = "Warmup iterations", value = rv_est_par$warmup),
          numericInput(session$ns("sampling"), label = "Sampling iterations", value = rv_est_par$sampling),
          numericInput(session$ns("chains"), label = "Number of chains", value = rv_est_par$chains),
          footer = tagList(
            actionButton(session$ns("act_advance_save"), "Save"),
            modalButton("Cancel")
          ),
          size = "l"
        ))
      })

      observeEvent(input$warmup, {
        warmup <- input$warmup
        if(!is.integer(warmup)){
          showFeedbackDanger(inputId = "warmup", text = "A positive integer is required!")
        }else if(warmup <= 0){
          showFeedbackDanger(inputId = "warmup", text = "A positive integer is required!")
        }else{
          hideFeedback("warmup")
        }
      })

      observeEvent(input$sampling, {
        sampling <- input$sampling
        if(!is.integer(sampling)){
          showFeedbackDanger(inputId = "sampling", text = "A positive integer is required!")
        }else if(sampling <= 0){
          showFeedbackDanger(inputId = "sampling", text = "A positive integer is required!")
        }else{
          hideFeedback("sampling")
        }
      })

      observeEvent(input$chains, {
        chains <- input$chains
        if(!is.integer(chains)){
          showFeedbackDanger(inputId = "chains", text = "A positive integer is required!")
        }else if(chains <= 0){
          showFeedbackDanger(inputId = "chains", text = "A positive integer is required!")
        }else{
          hideFeedback("chains")
        }
      })

      observeEvent(input$act_advance_save, {
        req(input$warmup)
        req(input$sampling)
        req(input$chains)

        isolate(warmup <- input$warmup)
        isolate(sampling <- input$sampling)
        isolate(chains <- input$chains)

        if(all(is.integer(warmup), is.integer(sampling), is.integer(chains))){
          if(all(warmup > 0, sampling > 0, chains > 0)){
            rv_est_par$warmup <- warmup
            rv_est_par$sampling <- sampling
            rv_est_par$chains <- chains
            removeModal()
          }
        }
      })

      # Start estimating parameters
      observeEvent(input$act_analyze_est_par, {
        # Fix data
        dt_exp <- rv_est_par$dt_exp
        dt_gt <- rv_est_par$dt_gt
        dt_setting_est_par <- hot_to_r(input$dt_setting_est_par)
        kit_name <- rv_est_par$kit_name
        warmup <- rv_est_par$warmup
        sampling <- rv_est_par$sampling
        chains <- rv_est_par$chains

        fn_stutter_file <- input$stutter_file
        load(paste0(path_pack, "/extdata/stutter_file/", fn_stutter_file))

        # Check input files and setting
        error_message <- check_est_par(dt_exp, dt_gt, dt_setting_est_par)
        if(error_message != ""){
          showModal(modalDialog(title = "Error", HTML(error_message), easyClose = TRUE, footer = NULL))
        }else{
          # Load the MCMC models
          model_ae <- stan_package_model(name = "model_ae", package = "KongohPlus")
          model_gamma <- stan_package_model(name = "model_gamma", package = "KongohPlus")

          rv_est_par$data_est_par <- estimate_parameters(dt_exp, dt_gt, list_dt_average_sr, dt_setting_est_par,
                                                         model_ae, model_gamma,
                                                         warmup, sampling, chains)
          rv_est_par$list_dt_average_sr <- list_dt_average_sr
          rv_est_par$kit_name <- kit_name
          rv_est_par$fn_stutter_file <- fn_stutter_file

          updateTabsetPanel(session = session_param_top, inputId = "tabset_est_par", selected = "Result")
        }
      })

      return(rv_est_par)
    }
  )
}
