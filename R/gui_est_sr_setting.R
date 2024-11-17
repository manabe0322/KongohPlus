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
#' @param path_pack The package path of KongohPlus
create_init_dt_setting_est_sr <- function(dt_kit, path_pack){
  dt_kit_auto <- extract_kit_auto(dt_kit)
  markers <- unique(dt_kit_auto$Marker)
  dyes <- extract_dye(dt_kit_auto, markers)

  dt_setting_est_sr <- data.table(Marker = markers, Dye = dyes)
  dt_setting_est_sr$AT <- 50
  dt_setting_est_sr$Model_BSR <- "Allele"
  dt_setting_est_sr$Model_FSR <- "Uniform (Multiple loci)"
  dt_setting_est_sr$Model_DSR <- "Uniform (Multiple loci)"
  dt_setting_est_sr$Model_M2SR <- "Not consider"
  dt_setting_est_sr$X_Multi <- NA

  names(dt_setting_est_sr) <- c("Marker", "Dye", "AT", "BSR model", "FSR model", "DSR model", "M2SR model", "Min. number of repeats that affects stuttering (for Multi-seq)")

  return(dt_setting_est_sr)
}

#' est_sr_setting_server
#'
#' @description The function to create the server module to create the setting tab for estimating average stutter ratios
est_sr_setting_server <- function(id, path_pack){
  moduleServer(
    id,
    function(input, output, session){
      rv_est_sr <- reactiveValues()
      rv_est_sr$dt_exp <- NULL
      rv_est_sr$dt_gt <- NULL
      rv_est_sr$dt_seq <- NULL
      rv_est_sr$dt_kit <- NULL
      rv_est_sr$results_est_sr <- NULL

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

      # Select kit
      output$kit_sr <- renderUI({
        kit_names <- list.files(paste0(path_pack, "/extdata/kit"))
        kit_names <- gsub(".csv", "", kit_names)
        selectInput(session$ns("kit_name"), label = h4("Kit"), choices = kit_names, selected = kit_names[1], multiple = FALSE, selectize = FALSE)
      })

      observeEvent(input$kit_name, ignoreInit = TRUE, {
        dt_kit <- fread(paste0(path_pack, "/extdata/kit/", input$kit_name, ".csv"))
        rv_est_sr$dt_kit <- dt_kit <- change_data_type_kit(dt_kit)
        dt_setting_est_sr <- create_init_dt_setting_est_sr(dt_kit, path_pack)

        output$dt_setting_est_sr <- renderRHandsontable({
          rhandsontable(dt_setting_est_sr, rowHeaders = FALSE) %>%
            hot_table(fillHandle = FALSE, rowHeights = 30) %>%
            hot_col(c("Marker"), type = "text", readOnly = TRUE, width = 100) %>%
            hot_col(c("Dye"), type = "text", readOnly = TRUE, width = 50) %>%
            hot_col(c("AT"), type = "numeric", format = "0", width = 50) %>%
            hot_col(c("BSR model", "FSR model", "DSR model", "M2SR model"), type = "dropdown", source = c("Allele", "LUS", "Multi-seq", "Uniform", "Uniform (Multiple loci)", "Not consider"), width = 175) %>%
            hot_col(c("Min. number of repeats that affects stuttering (for Multi-seq)"), type = "numeric", format = "0", width = 450)
        })
      })

      # Advanced setting
      rv_advanced <- reactiveValues()
      rv_advanced$warmup <- 2000
      rv_advanced$sampling <- 2000
      rv_advanced$chains <- 4
      rv_advanced$pu_remove <- TRUE
      rv_advanced$pu_dist <- 1
      rv_advanced$at_sex <- 100

      observeEvent(input$act_advance_sr, {
        showModal(modalDialog(
          title = "Advanced setting",
          fluidRow(
            column(6,
                   h4("MCMC conditions"),
                   numericInput(session$ns("warmup"), label = "Warmup iterations", value = rv_advanced$warmup),
                   numericInput(session$ns("sampling"), label = "Sampling iterations", value = rv_advanced$sampling),
                   numericInput(session$ns("chains"), label = "Number of chains", value = rv_advanced$chains)
            ),
            column(6,
                   h4("Pull-up"),
                   checkboxInput(session$ns("pu_remove"), label = "Automatically removed", value = rv_advanced$pu_remove),
                   numericInput(session$ns("pu_dist"), label = "Maximum pull-up distance (bp)", value = rv_advanced$pu_dist),
                   numericInput(session$ns("at_sex"), label = "AT of the sex-chromosomal markers", value = rv_advanced$at_sex)
            )
          ),
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

      observeEvent(input$pu_dist, {
        pu_dist <- input$pu_dist
        if(!is.numeric(pu_dist)){
          showFeedbackDanger(inputId = "pu_dist", text = "A positive number is required!")
        }else if(pu_dist <= 0){
          showFeedbackDanger(inputId = "pu_dist", text = "A positive number is required!")
        }else{
          hideFeedback("pu_dist")
        }
      })

      observeEvent(input$at_sex, {
        at_sex <- input$at_sex
        if(!is.numeric(at_sex)){
          showFeedbackDanger(inputId = "at_sex", text = "A positive number is required!")
        }else if(at_sex <= 0){
          showFeedbackDanger(inputId = "at_sex", text = "A positive number is required!")
        }else{
          hideFeedback("at_sex")
        }
      })

      observeEvent(input$act_advance_save, {
        req(input$warmup)
        req(input$sampling)
        req(input$chains)
        req(input$pu_dist)
        req(input$at_sex)

        isolate(warmup <- input$warmup)
        isolate(sampling <- input$sampling)
        isolate(chains <- input$chains)
        isolate(pu_remove <- input$pu_remove)
        isolate(pu_dist <- input$pu_dist)
        isolate(at_sex <- input$at_sex)

        if(all(is.integer(warmup), is.integer(sampling), is.integer(chains), is.numeric(pu_dist), is.numeric(at_sex))){
          if(all(warmup > 0, sampling > 0, chains > 0, pu_dist > 0, at_sex > 0)){
            rv_advanced$warmup <- warmup
            rv_advanced$sampling <- sampling
            rv_advanced$chains <- chains
            rv_advanced$pu_remove <- pu_remove
            rv_advanced$pu_dist <- pu_dist
            rv_advanced$at_sex <- at_sex
            removeModal()
          }
        }
      })

      # Start estimating average stutter ratios
      observeEvent(input$act_analyze_est_sr, {
        dt_exp <- rv_est_sr$dt_exp
        dt_gt <- rv_est_sr$dt_gt
        dt_seq <- rv_est_sr$dt_seq
        dt_kit <- rv_est_sr$dt_kit
        dt_setting_est_sr <- hot_to_r(input$dt_setting_est_sr)
        kit_name <- input$kit_name
        warmup <- rv_advanced$warmup
        sampling <- rv_advanced$sampling
        chains <- rv_advanced$chains
        pu_remove <- rv_advanced$pu_remove
        pu_dist <- rv_advanced$pu_dist
        at_sex <- rv_advanced$at_sex

        error_message <- check_est_sr(dt_exp, dt_gt, dt_seq, dt_setting_est_sr)
        if(error_message != ""){
          showModal(modalDialog(title = "Error", HTML(error_message), easyClose = TRUE, footer = NULL))
        }else{
          model_sr_gen <- stan_package_model(name = "model_sr_gen", package = "KongohPlus")
          model_sr_uni <- stan_package_model(name = "model_sr_uni", package = "KongohPlus")

          estimate_average_sr(dt_exp, dt_gt, dt_seq, dt_kit, dt_setting_est_sr, kit_name)
        }
      })

      return(rv_est_sr)
    }
  )
}
