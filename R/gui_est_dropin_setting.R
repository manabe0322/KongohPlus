#' est_dropin_setting_ui
#'
#' @description The function to create the ui module for setting to estimate drop-in parameters
est_dropin_setting_ui <- function(id){
  ns <- NS(id)

  tabPanel("Setting",
           br(),
           sidebarLayout(
             sidebarPanel(width = 3,
                          fileInput(ns("ntc_data"), label = h4("Non template control data"), accept = ".csv", buttonLabel = "Load"),
                          br(),
                          h4("Kit"),
                          fluidRow(
                            column(3, br(), shinyFilesButton(ns("load_kit"), label = "Load", title = "Kit", multiple = FALSE)),
                            column(9, uiOutput(ns("name_kit")))
                          ),
                          br(),
                          actionButton(ns("act_advance_dropin"), "Advanced setting"),
                          br(),
                          br(),
                          br(),
                          actionButton(ns("act_analyze_est_dropin"), label = "Analyze", class = "btn btn-primary btn-lg")
             ),
             mainPanel(width = 9,
                       tags$style(HTML("
                         .handsontable .htCore td {
                           vertical-align: middle !important;
                           text-align: center !important;
                         }"
                       )),
                       rHandsontableOutput(ns("dt_setting_est_dropin")),
                       br()
             )
           )
  )
}

#' create_init_dt_setting_est_dropin
#'
#' @description The function to create the initial data.table in the setting tab for estimating dropin parameters
#' @param dt_kit The data.table of information on a kit
create_init_dt_setting_est_dropin <- function(dt_kit){
  dt_kit_auto <- extract_kit_auto(dt_kit)
  markers <- unique(dt_kit_auto$Marker)

  dt_setting_est_dropin <- data.table(Marker = markers,
                                      AT = rep(45, 21))
  col_int <- c("AT")
  dt_setting_est_dropin[, (col_int) := lapply(.SD, as.integer), .SDcols = col_int]

  return(dt_setting_est_dropin)
}

#' est_dropin_setting_server
#'
#' @description The function to create the server module for setting to estimate drop-in parameters
est_dropin_setting_server <- function(id, path_pack, session_dropin_top){
  moduleServer(
    id,
    function(input, output, session){
      rv_est_dropin <- reactiveValues()
      rv_est_dropin$dt_ntc <- NULL
      rv_est_dropin$fn_ntc <- NULL
      rv_est_dropin$dt_setting_est_dropin <- NULL
      rv_est_dropin$dt_kit <- NULL
      rv_est_dropin$name_kit <- NULL
      rv_est_dropin$warmup <- 2000
      rv_est_dropin$sampling <- 2000
      rv_est_dropin$chains <- 4
      rv_est_dropin$data_est_dropin <- NULL

      # Load non template control data
      observeEvent(input$ntc_data, {
        file_input <- input$ntc_data
        if(!is.null(file_input)){
          dt_ntc <- fread(file_input$datapath)
          dt_ntc <- change_data_type_exp(dt_ntc)
          rv_est_dropin$dt_ntc <- dt_ntc
          rv_est_dropin$fn_ntc <- file_input$name
        }
      })

      # Select a kit
      shinyFileChoose(input, "load_kit", roots = c(documents = paste0(path_pack, "/extdata/kit")), session = session, filetypes = c("csv"))

      observeEvent(input$load_kit, {
        req(input$load_kit)
        path_kit <- parseFilePaths(roots = c(documents = paste0(path_pack, "/extdata/kit")), input$load_kit)$datapath
        if(length(path_kit) > 0){
          dt_kit <- fread(path_kit)
          rv_est_dropin$dt_kit <- dt_kit <- change_data_type_kit(dt_kit)
          rv_est_dropin$dt_setting_est_dropin <- create_init_dt_setting_est_dropin(dt_kit)
          rv_est_dropin$name_kit <- str_replace(basename(path_kit), ".{4}$", "")
        }
      })

      output$name_kit <- renderUI({
        disabled(textInput(session$ns("name_kit"), label = "", value = rv_est_sr$name_kit))
      })

      output$dt_setting_est_dropin <- renderRHandsontable({
        dt_setting_est_dropin <- rv_est_dropin$dt_setting_est_dropin
        if(!is.null(dt_setting_est_dropin)){
          rhandsontable(dt_setting_est_dropin, rowHeaders = FALSE) %>%
            hot_table(fillHandle = FALSE, rowHeights = 30) %>%
            hot_col(c("Marker"), type = "text", readOnly = TRUE, width = 100) %>%
            hot_col(c("AT"), type = "numeric", format = "0", width = 50)
        }
      })

      #### Advanced setting ####

      observeEvent(input$act_advance_dropin, {
        showModal(modalDialog(
          title = "Advanced setting",
          h4("MCMC conditions"),
          numericInput(session$ns("warmup"), label = "Warmup iterations", value = rv_est_dropin$warmup),
          numericInput(session$ns("sampling"), label = "Sampling iterations", value = rv_est_dropin$sampling),
          numericInput(session$ns("chains"), label = "Number of chains", value = rv_est_dropin$chains),
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
            rv_est_dropin$warmup <- warmup
            rv_est_dropin$sampling <- sampling
            rv_est_dropin$chains <- chains
            removeModal()
          }
        }
      })

      # Start estimating dropin parameters
      observeEvent(input$act_analyze_est_dropin, {
        # Fix data
        dt_ntc <- rv_est_dropin$dt_ntc
        rv_est_dropin$dt_setting_est_dropin <- dt_setting_est_dropin <- hot_to_r(input$dt_setting_est_dropin)
        warmup <- rv_est_dropin$warmup
        sampling <- rv_est_dropin$sampling
        chains <- rv_est_dropin$chains

        # Check input files and setting
        error_message <- check_est_dropin(dt_ntc, dt_setting_est_dropin)
        if(error_message != ""){
          showModal(modalDialog(title = "Error", HTML(error_message), easyClose = TRUE, footer = NULL))
        }else{

        }
      })

      return(rv_est_dropin)
    }
  )
}
