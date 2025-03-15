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
                          uiOutput(ns("seq_data")),
                          uiOutput(ns("kit_sr")),
                          br(),
                          actionButton(ns("act_advance_sr"), "Advanced setting"),
                          br(),
                          br(),
                          br(),
                          actionButton(ns("act_analyze_est_sr"), label = "Analyze", class = "btn btn-primary btn-lg")
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
    dt_setting_est_sr$Model_BSR <- c("Allele", "Allele", "Allele", "Allele", "Allele",
                                     "LUS", "LUS", "Allele",
                                     "LUS", "Allele", "LUS", "Allele",
                                     "Allele", "Allele", "Allele", "Allele", "Multi-seq",
                                     "Allele", "LUS", "LUS", "LUS")
    dt_setting_est_sr$Model_FSR <- "Uniform (Multiple loci)"
    dt_setting_est_sr$Model_FSR[13] <- "Allele"
    dt_setting_est_sr$Model_DSR <- "Uniform (Multiple loci)"
    dt_setting_est_sr$Model_DSR[13] <- "Allele"
    dt_setting_est_sr$Model_M2SR <- "Not consider"
    dt_setting_est_sr$Model_M2SR[17] <- "Uniform"
    dt_setting_est_sr$Model_M2SR[19] <- "Uniform"
    dt_setting_est_sr$X_Multi <- NA
    dt_setting_est_sr$X_Multi[17] <- 7
  }else{
    dt_setting_est_sr$AT <- 50
    dt_setting_est_sr$Model_BSR <- "Allele"
    dt_setting_est_sr$Model_FSR <- "Uniform (Multiple loci)"
    dt_setting_est_sr$Model_DSR <- "Uniform (Multiple loci)"
    dt_setting_est_sr$Model_M2SR <- "Not consider"
    dt_setting_est_sr$X_Multi <- NA
  }

  col_int <- c("AT", "X_Multi")
  dt_setting_est_sr[, (col_int) := lapply(.SD, as.integer), .SDcols = col_int]

  names(dt_setting_est_sr) <- c("Marker", "AT", "BSR model", "FSR model", "DSR model", "M2SR model", "Minimum uninterrupted stretch (for Multi-seq)")

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
      rv_est_sr$dt_gt <- NULL
      rv_est_sr$dt_seq <- NULL
      rv_est_sr$dt_setting_est_sr <- NULL
      rv_est_sr$kit_name <- NULL
      rv_est_sr$warmup <- 2000
      rv_est_sr$sampling <- 2000
      rv_est_sr$chains <- 4
      rv_est_sr$data_est_sr <- NULL

      # Load experimental data
      observeEvent(input$exp_data, {
        file_input <- input$exp_data
        if(!is.null(file_input)){
          dt_exp <- fread(file_input$datapath)
          dt_exp <- change_data_type_exp(dt_exp)
          rv_est_sr$dt_exp <- dt_exp
        }
      })

      # Load genotype data
      observeEvent(input$gt_data, {
        file_input <- input$gt_data
        if(!is.null(file_input)){
          dt_gt <- fread(file_input$datapath)
          dt_gt <- change_data_type_gt(dt_gt)
          rv_est_sr$dt_gt <- dt_gt
        }
      })

      # Select sequence data
      output$seq_data <- renderUI({
        ns <- session$ns
        seq_files <- list.files(paste0(path_pack, "/extdata/seq"))
        selectInput(ns("seq_data"), label = h4("Sequence data"), choices = seq_files, selected = seq_files[1], multiple = FALSE, selectize = FALSE)
      })

      observeEvent(input$seq_data, {
        file_input <- input$seq_data
        if(!is.null(file_input)){
          dt_seq <- fread(paste0(path_pack, "/extdata/seq/", file_input))
          dt_seq <- change_data_type_seq(dt_seq)
          rv_est_sr$dt_seq <- dt_seq
        }
      })

      # Select a kit
      output$kit_sr <- renderUI({
        kit_names <- list.files(paste0(path_pack, "/extdata/kit"))
        kit_names <- gsub(".csv", "", kit_names)
        selectInput(session$ns("kit_name"), label = h4("Kit"), choices = kit_names, selected = kit_names[1], multiple = FALSE, selectize = FALSE)
      })

      observeEvent(input$kit_name, ignoreInit = TRUE, {
        rv_est_sr$kit_name <- kit_name <- input$kit_name
        dt_kit <- fread(paste0(path_pack, "/extdata/kit/", kit_name, ".csv"))
        dt_kit <- change_data_type_kit(dt_kit)
        dt_setting_est_sr <- create_init_dt_setting_est_sr(dt_kit)

        output$dt_setting_est_sr <- renderRHandsontable({
          rhandsontable(dt_setting_est_sr, rowHeaders = FALSE) %>%
            hot_table(fillHandle = FALSE, rowHeights = 30) %>%
            hot_col(c("Marker"), type = "text", readOnly = TRUE, width = 100) %>%
            hot_col(c("AT"), type = "numeric", format = "0", width = 50) %>%
            hot_col(c("BSR model", "FSR model", "DSR model", "M2SR model"), type = "dropdown", source = c("Allele", "LUS", "Multi-seq", "Uniform", "Uniform (Multiple loci)", "Not consider"), width = 175) %>%
            hot_col(c("Minimum uninterrupted stretch (for Multi-seq)"), type = "numeric", format = "0", width = 450)
        })
      })

      observeEvent(input$act_advance_sr, {
        showModal(modalDialog(
          title = "Advanced setting",
          h4("MCMC conditions"),
          numericInput(session$ns("warmup"), label = "Warmup iterations", value = rv_est_sr$warmup),
          numericInput(session$ns("sampling"), label = "Sampling iterations", value = rv_est_sr$sampling),
          numericInput(session$ns("chains"), label = "Number of chains", value = rv_est_sr$chains),
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
            rv_est_sr$warmup <- warmup
            rv_est_sr$sampling <- sampling
            rv_est_sr$chains <- chains
            removeModal()
          }
        }
      })

      # Start estimating average stutter ratios
      observeEvent(input$act_analyze_est_sr, {
        dt_exp <- rv_est_sr$dt_exp
        dt_gt <- rv_est_sr$dt_gt
        dt_seq <- rv_est_sr$dt_seq
        rv_est_sr$dt_setting_est_sr <- dt_setting_est_sr <- hot_to_r(input$dt_setting_est_sr)
        kit_name <- rv_est_sr$kit_name
        warmup <- rv_est_sr$warmup
        sampling <- rv_est_sr$sampling
        chains <- rv_est_sr$chains

        error_message <- check_est_sr(dt_exp, dt_gt, dt_seq, dt_setting_est_sr)
        if(error_message != ""){
          showModal(modalDialog(title = "Error", HTML(error_message), easyClose = TRUE, footer = NULL))
        }else{
          # Load the MCMC models
          model_sr_gen <- stan_package_model(name = "model_sr_gen", package = "KongohPlus")
          model_sr_uni <- stan_package_model(name = "model_sr_uni", package = "KongohPlus")

          rv_est_sr$data_est_sr <- estimate_average_sr(dt_exp, dt_gt, dt_seq, dt_setting_est_sr, kit_name,
                                                       model_sr_gen, model_sr_uni,
                                                       warmup, sampling, chains)

          updateTabsetPanel(session = session_est_sr_top, inputId = "tabset_est_sr", selected = "Result")
        }
      })

      return(rv_est_sr)
    }
  )
}
