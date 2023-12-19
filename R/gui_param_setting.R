#' est_par_setting_ui
#'
#' @description The function to create the ui module for setting to estimate parameters
est_par_setting_ui <- function(id){
  ns <- NS(id)

  tabPanel("Setting",
           br(),

           tags$style(
             HTML(
               "td > div.form-group {
                  margin: 0px;
                  height: 30px;
                  display: flex;
                  align-items: center;
                }
                td:nth-child(1) {
                  width: 100px;
                }
                td:nth-child(2) {
                  width: 50px;
                }
                td:nth-child(3) {
                  width: 100px;
                }
                td:nth-child(4) {
                  width: 150px;
                }
                td:nth-child(5) {
                  width: 150px;
                }
                td:nth-child(6) {
                  width: 150px;
                }
                td:nth-child(7) {
                  width: 150px;
                }"
             )
           ),

           sidebarLayout(
             sidebarPanel(width = 3,
                          fileInput(ns("exp_data"), label = h4("Experimental data"), accept = ".csv"),
                          fileInput(ns("gt_data"), label = h4("Genotype data"), accept = ".csv"),
                          uiOutput(ns("seq_data")),
                          uiOutput(ns("kit_par")),
                          br(),
                          actionButton(ns("act_advanced"), "Advanced setting"),
                          br(),
                          br(),
                          br(),
                          actionButton(ns("act_analyze_est_par"), label = "Analyze", class = "btn btn-primary btn-lg")
             ),
             mainPanel(width = 9,
                       dataTableOutput(ns("dt_kit_set"))
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
      selectColumn <- function(len, col, choices, select_default, ...){
        inputs <- character(len)
        for(i in seq_len(len)){
          inputs[i] <- as.character(
            selectizeInput(paste0("select_", col, "_", i), label = NULL, choices = choices, selected = select_default[i], ...)
          )
        }
        return(inputs)
      }

      numericColumn <- function(len, col, ...){
        inputs <- character(len)
        for(i in seq_len(len)){
          inputs[i] <- as.character(
            numericInput(paste0("numeric_", col, "_", i), label = NULL, ...)
          )
        }
        return(inputs)
      }

      js <- function(dtid, col_numeric, cols_check){
        code <- vector("list", 1 + length(cols_check))
        code[[1]] <- c(
          sprintf("$('body').on('change', '[id^=numeric_%d_]', function() {", col_numeric),
          "  var id = this.getAttribute('id');",
          sprintf("  var i = parseInt(/numeric_%d_(\\d+)/.exec(id)[1]);", col_numeric),
          "  var value = $(this).val();",
          sprintf("  var info = [{row: i, col: %d, value: value}];", col_numeric),
          sprintf("  Shiny.setInputValue('%s', info);", session$ns(sprintf("%s_cell_edit:DT.cellInfo", dtid))),
          "});"
        )

        for(i in seq_along(cols_check)){
          col <- cols_check[i]
          code[[i + 1]] <- c(
            sprintf("$('body').on('change', '[id^=select_%d_]', function() {", col),
            "  var id = this.getAttribute('id');",
            sprintf("  var i = parseInt(/select_%d_(\\d+)/.exec(id)[1]);", col),
            "  var value = $(this).val();",
            sprintf("  var info = [{row: i, col: %d, value: value}];", col),
            sprintf("  Shiny.setInputValue('%s', info);", session$ns(sprintf("%s_cell_edit:DT.cellInfo", dtid))),
            "});"
          )
        }
        do.call(c, code)
      }

      rv_output <- reactiveValues()
      rv_output$dt_exp <- NULL
      rv_output$dt_gt <- NULL
      rv_output$dt_kit_set <- NULL
      rv_output$results_est_par <- NULL

      #### Experimental data ####

      observeEvent(input$exp_data, {
        file_input <- input$exp_data
        if(!is.null(file_input)){
          dt_exp <- fread(file_input$datapath)
          dt_exp <- change_data_type_exp(dt_exp)
          rv_output$dt_exp <- dt_exp
        }
      })

#      react_dt_exp <- reactive({
#        file_input <- input$exp_data
#        if(is.null(file_input)){
#          return(NULL)
#        }else{
#          dt_exp <- fread(file_input$datapath)
#          dt_exp <- change_data_type_exp(dt_exp)
#          return(dt_exp)
#        }
#      })

      #### Genotype data ####

      observeEvent(input$gt_data, {
        file_input <- input$gt_data
        if(!is.null(file_input)){
          dt_gt <- fread(file_input$datapath)
          dt_gt <- change_data_type_gt(dt_gt)
          rv_output$dt_gt <- dt_gt
        }
      })

#      react_dt_gt <- reactive({
#        file_input <- input$gt_data
#        if(is.null(file_input)){
#          return(NULL)
#        }else{
#          dt_gt <- fread(file_input$datapath)
#          dt_gt <- change_data_type_gt(dt_gt)
#          return(dt_gt)
#        }
#      })

      #### Sequence data ####

      output$seq_data <- renderUI({
        ns <- session$ns
        seq_files <- list.files(paste0(path_pack, "/extdata/sequence_data"))
        selectInput(ns("seq_data"), label = h4("Sequence data"), choices = seq_files, selected = seq_files[1], multiple = FALSE)
      })

      observeEvent(input$seq_data, {
        file_input <- input$seq_data
        if(!is.null(file_input)){
          dt_seq <- fread(paste0(path_pack, "/extdata/sequence_data/", file_input))
          dt_seq <- change_data_type_seq(dt_seq)
          rv_output$dt_seq <- dt_seq
        }
      })

#      react_dt_seq <- reactive({
#        file_input <- input$seq_data
#        if(is.null(file_input)){
#          return(NULL)
#        }else{
#          dt_seq <- fread(paste0(path_pack, "/extdata/sequence_data/", file_input))
#          dt_seq <- change_data_type_seq(dt_seq)
#          return(dt_seq)
#        }
#      })

      #### Kit-setting ####

      # Initial data.table
      react_dt_kit_set0 <- reactive({
        # Extract autosomal markers
        kit <- input$kit_par
        dt_kit_info <- fread(paste0(path_pack, "/extdata/kit/", kit, ".csv"))
        setkey(dt_kit_info, Sex_Chr)
        dt_kit_auto <- dt_kit_info[""]
        markers_par <- dt_kit_auto[, Marker]
        dyes_par <- dt_kit_auto[, Dye]

        # Make a data.table for kit-setting
        dt_kit_set0 <- data.table(Marker = markers_par, Dye = dyes_par)
        return(dt_kit_set0)
      })

      # Main data.table
      react_dt_kit_set1 <- reactive({
        dt_kit_set1 <- react_dt_kit_set0()
        dt_kit_set1$AT <- 50
        dt_kit_set1$Model_BSR <- "Locus-specific"
        dt_kit_set1$Model_FSR <- "Multiple loci together"
        dt_kit_set1$Model_DSR <- "Multiple loci together"
        dt_kit_set1$Model_M2SR <- "Not consider"
        return(dt_kit_set1)
      })

      # Displayed data.table
      react_dt_kit_set2 <- reactive({
        dt_kit_set2 <- react_dt_kit_set0()
        n_mk <- nrow(dt_kit_set2)
        dt_kit_set2$numeric_at <- numericColumn(n_mk, 3, value = 50, width = "75%")
        dt_kit_set2$select_bsr <- selectColumn(n_mk, 4, choices = c("Not consider", "Locus-specific", "Multiple loci together"), select_default = rep("Locus-specific", n_mk), width = "90%")
        dt_kit_set2$select_fsr <- selectColumn(n_mk, 5, choices = c("Not consider", "Locus-specific", "Multiple loci together"), select_default = rep("Multiple loci together", n_mk), width = "90%")
        dt_kit_set2$select_dsr <- selectColumn(n_mk, 6, choices = c("Not consider", "Locus-specific", "Multiple loci together"), select_default = rep("Multiple loci together", n_mk), width = "90%")
        dt_kit_set2$select_m2sr <- selectColumn(n_mk, 7, choices = c("Not consider", "Locus-specific", "Multiple loci together"), select_default = rep("Not consider", n_mk), width = "90%")
        return(dt_kit_set2)
      })

      # Select kit
      output$kit_par <- renderUI({
        kit_names <- list.files(paste0(path_pack, "/extdata/kit"))
        kit_names <- gsub(".csv", "", kit_names)
        selectInput(session$ns("kit_par"), label = h4("Kit"), choices = kit_names, selected = kit_names[1], multiple = FALSE)
      })

      # Output data.table
      observeEvent(input$kit_par, ignoreInit = TRUE, {
        rv_output$dt_kit_set <- react_dt_kit_set1()

        output$dt_kit_set <- renderDataTable(server = FALSE, {
          datatable(
            isolate(react_dt_kit_set2()),
            colnames = c("Marker", "Dye", "AT", "BSR model", "FSR model", "DSR model", "M2SR model"),
            selection = "none",
            options = list(iDisplayLength = 50, ordering = FALSE, autowidth = TRUE),
            rownames = FALSE,
            escape = FALSE,
            callback = JS(js("dt_kit_set", 3, c(4, 5, 6, 7)))
          )
        })
      })

      # Reflect changing to the main data.table
      observeEvent(input[["dt_kit_set_cell_edit"]], {
        info <- input[["dt_kit_set_cell_edit"]]
        dt_kit_set <- rv_output$dt_kit_set
        dt_kit_set[[info[1, "col"]]][info[1, "row"]] <- info[1, "value"]
        rv_output$dt_kit_set <- dt_kit_set
      })

      #### Advanced setting ####

      # Define reactive values
      rv_advanced <- reactiveValues()
      rv_advanced$warmup <- 2000
      rv_advanced$sampling <- 2000
      rv_advanced$chains <- 4
      rv_advanced$pu_auto <- TRUE
      rv_advanced$pu_dist <- 1
      rv_advanced$motif_lus <- 2
      rv_advanced$motif_multi <- 3

      # Make the setting pop-up window
      observeEvent(input$act_advanced, {
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
                   checkboxInput(session$ns("pu_auto"), label = "Automatically removed", value = rv_advanced$pu_auto),
                   numericInput(session$ns("pu_dist"), label = "Maximum pull-up distance (bp)", value = rv_advanced$pu_dist),
                   br(),
                   h4("Stutter model"),
                   numericInput(session$ns("motif_lus"), label = "Minimum number of repeat motifs for considering the LUS model", value = rv_advanced$motif_lus),
                   numericInput(session$ns("motif_multi"), label = "Minimum number of repeat motifs for considering the multi-sequence model", value = rv_advanced$motif_multi)
            )
          ),
          footer = tagList(
            actionButton(session$ns("act_advance_save"), "Save"),
            modalButton("Cancel")
          ),
          size = "l"
        ))

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

        observeEvent(input$motif_lus, {
          motif_lus <- input$motif_lus
          if(!is.integer(motif_lus)){
            showFeedbackDanger(inputId = "motif_lus", text = "A positive integer is required!")
          }else if(motif_lus <= 0){
            showFeedbackDanger(inputId = "motif_lus", text = "A positive integer is required!")
          }else{
            hideFeedback("motif_lus")
          }
        })

        observeEvent(input$motif_multi, {
          motif_multi <- input$motif_multi
          if(!is.integer(motif_multi)){
            showFeedbackDanger(inputId = "motif_multi", text = "A positive integer is required!")
          }else if(motif_multi <= 0){
            showFeedbackDanger(inputId = "motif_multi", text = "A positive integer is required!")
          }else{
            hideFeedback("motif_multi")
          }
        })
      })

      # Save the setting
      observeEvent(input$act_advance_save, {
        req(input$warmup)
        req(input$sampling)
        req(input$chains)
        req(input$pu_dist)
        req(input$motif_lus)
        req(input$motif_multi)

        isolate(warmup <- input$warmup)
        isolate(sampling <- input$sampling)
        isolate(chains <- input$chains)
        isolate(pu_auto <- input$pu_auto)
        isolate(pu_dist <- input$pu_dist)
        isolate(motif_lus <- input$motif_lus)
        isolate(motif_multi <- input$motif_multi)

        if(all(is.integer(warmup), is.integer(sampling), is.integer(chains), is.numeric(pu_dist), is.integer(motif_lus), is.integer(motif_multi))){
          rv_advanced$warmup <- warmup
          rv_advanced$sampling <- sampling
          rv_advanced$chains <- chains
          rv_advanced$pu_auto <- pu_auto
          rv_advanced$pu_dist <- pu_dist
          rv_advanced$motif_lus <- motif_lus
          rv_advanced$motif_multi <- motif_multi

          removeModal()
        }
      })

      #### Analyze ####

      observeEvent(input$act_analyze_est_par, {
        # Fix data.table
        dt_exp <- rv_output$dt_exp
        dt_gt <- rv_output$dt_gt
        dt_seq <- rv_output$dt_seq
        dt_kit_set <- rv_output$dt_kit_set

        # Fix parameters
        n_warmup <- rv_advanced$warmup
        n_sampling <- rv_advanced$sampling
        n_chains <- rv_advanced$chains
        pu_dist <- rv_advanced$pu_dist
        motif_lus <- rv_advanced$motif_lus
        motif_multi <- rv_advanced$motif_multi

        # Check input files and setting
        error_message <- check_est_par(dt_exp, dt_gt, dt_seq, dt_kit_set)
        if(error_message != ""){
          showModal(modalDialog(title = "Error", HTML(error_message), easyClose = TRUE, footer = NULL))
        }else{
          # Load the MCMC models
          model_ae <- readRDS(paste0(path_pack, "/model/model_ae.rds"))
          model_sr_gen <- readRDS(paste0(path_pack, "/model/model_sr_gen.rds"))
          model_sr_together <- readRDS(paste0(path_pack, "/model/model_sr_together.rds"))
          model_sr_uni <- readRDS(paste0(path_pack, "/model/model_sr_uni.rds"))

          results_est_par <- estimate_parameters(dt_exp, dt_gt, dt_seq, dt_kit_set, model_ae, model_sr_gen, model_sr_together, model_sr_uni, n_warmup, n_sampling, n_chains)
          rv_output$results_est_par <- results_est_par

          updateTabsetPanel(session = session_param_top, inputId = "tabset_est_par", selected = "Result")
        }
      })

      return(rv_output)
    }
  )
}
