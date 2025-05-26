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
                            column(9, uiOutput(ns("name_kit")))
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
      rv_est_sr$dt_exp <- NULL
      rv_est_sr$fn_exp <- NULL
      rv_est_sr$dt_gt <- NULL
      rv_est_sr$fn_gt <- NULL
      rv_est_sr$dt_kit <- NULL
      rv_est_sr$dt_setting_est_sr <- NULL
      rv_est_sr$name_kit <- NULL
      rv_est_sr$n_conv <- 10
      rv_est_sr$n_max_init <- 100
      rv_est_sr$data_est_sr <- NULL

      # Load experimental data
      observeEvent(input$exp_data, {
        file_input <- input$exp_data
        if(!is.null(file_input)){
          dt_exp <- fread(file_input$datapath)
          dt_exp <- change_data_type_exp(dt_exp)
          rv_est_sr$dt_exp <- dt_exp
          rv_est_sr$fn_exp <- file_input$name
        }
      })

      # Load genotype data
      observeEvent(input$gt_data, {
        file_input <- input$gt_data
        if(!is.null(file_input)){
          dt_gt <- fread(file_input$datapath)
          dt_gt <- change_data_type_gt(dt_gt)
          rv_est_sr$dt_gt <- dt_gt
          rv_est_sr$fn_gt <- file_input$name
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
          rv_est_sr$name_kit <- str_replace(basename(path_kit), ".{4}$", "")
        }
      })

      output$name_kit <- renderUI({
        disabled(textInput(session$ns("name_kit"), label = "", value = rv_est_sr$name_kit))
      })

      observeEvent(input$act_advance_sr, {
        showModal(modalDialog(
          title = "Advanced setting",
          h4("MLE conditions"),
          numericInput(session$ns("n_conv"), label = "Number of successful convergences", value = rv_est_sr$n_conv),
          numericInput(session$ns("n_max_init"), label = "Maximum number of initializations", value = rv_est_sr$n_max_init),
          footer = tagList(
            actionButton(session$ns("act_advance_sr_save"), "Save"),
            modalButton("Cancel")
          ),
          size = "l"
        ))
      })

      observeEvent(input$n_conv | input$n_max_init, {
        n_conv <- input$n_conv
        n_max_init <- input$n_max_init

        hideFeedback("n_conv")
        hideFeedback("n_max_init")
        bool_correct_mle <- TRUE

        if(!is.integer(n_conv) || n_conv <= 0){
          showFeedbackDanger(inputId = "n_conv", text = "A positive integer is required!")
          bool_correct_mle <- FALSE
        }

        if(!is.integer(n_max_init) || n_max_init <= 0){
          showFeedbackDanger(inputId = "n_max_init", text = "A positive integer is required!")
          bool_correct_mle <- FALSE
        }

        if(is.integer(n_conv) && is.integer(n_max_init) && n_conv > n_max_init){
          showFeedbackDanger(inputId = "n_conv", text = "The number of successful convergences must not exceed the maximum number of initializations!")
          showFeedbackDanger(inputId = "n_max_init", text = "The number of successful convergences must not exceed the maximum number of initializations!")
          bool_correct_mle <- FALSE
        }

        if(bool_correct_mle){
          enable("act_advance_sr_save")
        }else{
          disable("act_advance_sr_save")
        }
      })

      observeEvent(input$act_advance_sr_save, {
        rv_est_sr$n_conv <- input$n_conv
        rv_est_sr$n_max_init <- input$n_max_init
        removeModal()
      })

      output$act_analyze_est_sr <- renderUI({
        dt_exp <- rv_est_sr$dt_exp
        dt_gt <- rv_est_sr$dt_gt
        dt_kit <- rv_est_sr$dt_kit
        if(any(is.null(dt_exp), is.null(dt_gt), is.null(dt_kit))){
          disabled(actionButton(session$ns("act_analyze_est_sr"), label = "Analyze", class = "btn btn-primary btn-lg"))
        }else{
          actionButton(session$ns("act_analyze_est_sr"), label = "Analyze", class = "btn btn-primary btn-lg")
        }
      })

      # Start estimating average stutter ratios
      observeEvent(input$act_analyze_est_sr, {
        dt_exp <- rv_est_sr$dt_exp
        dt_gt <- rv_est_sr$dt_gt
        dt_kit <- rv_est_sr$dt_kit
        name_kit <- rv_est_sr$name_kit
        rv_est_sr$dt_setting_est_sr <- dt_setting_est_sr <- hot_to_r(input$dt_setting_est_sr)
        n_conv <- rv_est_sr$n_conv
        n_max_init <- rv_est_sr$n_max_init

        error_message <- check_est_sr(dt_exp, dt_gt, dt_kit, dt_setting_est_sr)
        if(error_message != ""){
          showModal(modalDialog(title = "Error", HTML(error_message), easyClose = TRUE, footer = NULL))
        }else{
          # Load the stan models
          model_sr_gen <- stan_package_model(name = "model_sr_gen", package = "KongohPlus")
          model_sr_uni <- stan_package_model(name = "model_sr_uni", package = "KongohPlus")

          rv_est_sr$data_est_sr <- estimate_average_sr(dt_exp, dt_gt, dt_kit, dt_setting_est_sr, model_sr_gen, model_sr_uni, n_conv, n_max_init)

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
