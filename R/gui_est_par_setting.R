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
                          uiOutput(ns("name_kit")),
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
#' @description The function to create the initial data.table in the setting tab for estimating parameters
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
      rv_est_par$fn_exp <- NULL
      rv_est_par$dt_gt <- NULL
      rv_est_par$fn_gt <- NULL
      rv_est_par$dt_setting_est_par <- NULL
      rv_est_par$list_dt_average_sr <- NULL
      rv_est_par$name_kit <- NULL
      rv_est_par$name_sr <- NULL
      rv_est_par$n_conv <- 10
      rv_est_par$n_max_init <- 100
      rv_est_par$data_est_par <- NULL

      # Load experimental data
      observeEvent(input$exp_data, {
        file_input <- input$exp_data
        if(!is.null(file_input)){
          dt_exp <- fread(file_input$datapath)
          dt_exp <- change_data_type_exp(dt_exp)
          rv_est_par$dt_exp <- dt_exp
          rv_est_par$fn_exp <- file_input$name
        }
      })

      # Load genotype data
      observeEvent(input$gt_data, {
        file_input <- input$gt_data
        if(!is.null(file_input)){
          dt_gt <- fread(file_input$datapath)
          dt_gt <- change_data_type_gt(dt_gt)
          rv_est_par$dt_gt <- dt_gt
          rv_est_par$fn_gt <- file_input$name
        }
      })

      # Select a file for stutter ratios
      shinyFileChoose(input, "load_sr", roots = c(documents = paste0(path_pack, "/extdata/stutter")), session = session, filetypes = c("RData"))

      observeEvent(input$load_sr, {
        req(input$load_sr)
        path_sr <- parseFilePaths(roots = c(documents = paste0(path_pack, "/extdata/stutter")), input$load_sr)$datapath
        if(length(path_sr) > 0){
          load(path_sr)
          rv_est_par$list_dt_average_sr <- list_dt_average_sr
          rv_est_par$name_kit <- name_kit
          rv_est_par$name_sr <- str_replace(basename(path_sr), ".{6}$", "")
        }
      })

      output$name_sr <- renderUI({
        disabled(textInput(session$ns("name_sr"), label = "", value = rv_est_par$name_sr))
      })

      # Kit
      output$name_kit <- renderUI({
        disabled(textInput(session$ns("name_kit"), label = "", value = rv_est_par$name_kit))
      })

      observeEvent(input$name_kit, ignoreInit = TRUE, {
        name_kit <- rv_est_par$name_kit
        if(!is.null(name_kit)){
          dt_kit <- fread(paste0(path_pack, "/extdata/kit/", name_kit, ".csv"))
          dt_kit <- change_data_type_kit(dt_kit)
          rv_est_par$dt_setting_est_par <- dt_setting_est_par <- create_init_dt_setting_est_par(dt_kit)

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
          title = "Advanced setting",
          h4("MCMC conditions"),
          numericInput(session$ns("n_conv"), label = "Number of successful convergences", value = rv_est_par$n_conv),
          numericInput(session$ns("n_max_init"), label = "Maximum number of initializations", value = rv_est_par$n_max_init),
          footer = tagList(
            actionButton(session$ns("act_advance_par_save"), "Save"),
            modalButton("Cancel")
          ),
          size = "l"
        ))
      })

      observeEvent(input$n_conv, {
        n_conv <- input$n_conv
        if(!is.integer(n_conv)){
          showFeedbackDanger(inputId = "n_conv", text = "A positive integer is required!")
          disable("act_advance_par_save")
        }else if(n_conv <= 0){
          showFeedbackDanger(inputId = "n_conv", text = "A positive integer is required!")
          disable("act_advance_par_save")
        }else{
          hideFeedback("n_conv")
          enable("act_advance_par_save")
        }
      })

      observeEvent(input$n_max_init, {
        isolate(n_conv <- input$n_conv)
        n_max_init <- input$n_max_init
        if(!is.integer(n_max_init)){
          showFeedbackDanger(inputId = "n_max_init", text = "A positive integer is required!")
          disable("act_advance_par_save")
        }else if(n_max_init <= 0){
          showFeedbackDanger(inputId = "n_max_init", text = "A positive integer is required!")
          disable("act_advance_par_save")
        }else if(n_max_init < n_conv){
          showFeedbackDanger(inputId = "n_max_init", text = "The number must be greater than the number of successful convergences!")
          disable("act_advance_par_save")
        }else{
          hideFeedback("n_max_init")
          enable("act_advance_par_save")
        }
      })

      observeEvent(input$act_advance_par_save, {
        rv_est_par$n_conv <- input$n_conv
        rv_est_par$n_max_init <- input$n_max_init
        removeModal()
      })

      # Start estimating parameters
      observeEvent(input$act_analyze_est_par, {
        # Fix data
        dt_exp <- rv_est_par$dt_exp
        dt_gt <- rv_est_par$dt_gt
        rv_est_par$dt_setting_est_par <- dt_setting_est_par <- hot_to_r(input$dt_setting_est_par)
        list_dt_average_sr <- rv_est_par$list_dt_average_sr
        n_conv <- rv_est_par$n_conv
        n_max_init <- rv_est_par$n_max_init

        # Check input files and setting
        error_message <- check_est_par(dt_exp, dt_gt, dt_setting_est_par)
        if(error_message != ""){
          showModal(modalDialog(title = "Error", HTML(error_message), easyClose = TRUE, footer = NULL))
        }else{
          # Load the MCMC models
          model_ae <- stan_package_model(name = "model_ae", package = "KongohPlus")
          model_gamma <- stan_package_model(name = "model_gamma", package = "KongohPlus")

          rv_est_par$data_est_par <- estimate_parameters(dt_exp, dt_gt, list_dt_average_sr, dt_setting_est_par, model_ae, model_gamma, n_conv, n_max_init)

          updateTabsetPanel(session = session_param_top, inputId = "tabset_est_par", selected = "Result")
        }
      })

      return(rv_est_par)
    }
  )
}
