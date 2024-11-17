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
                          fileInput(ns("exp_data"), label = h4("Experimental data"), accept = ".csv", buttonLabel = "Load"),
                          fileInput(ns("gt_data"), label = h4("Genotype data"), accept = ".csv", buttonLabel = "Load"),
                          uiOutput(ns("seq_data")),
                          uiOutput(ns("kit_par")),
                          br(),
                          actionButton(ns("act_advance_par"), "Advanced setting"),
                          br(),
                          br(),
                          br(),
                          actionButton(ns("act_analyze_est_par"), label = "Analyze", class = "btn btn-primary btn-lg")
             ),
             mainPanel(width = 9,
                       dataTableOutput(ns("dt_kit_set")),
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
      rv_output$dt_kit <- NULL
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

      #### Genotype data ####

      observeEvent(input$gt_data, {
        file_input <- input$gt_data
        if(!is.null(file_input)){
          dt_gt <- fread(file_input$datapath)
          dt_gt <- change_data_type_gt(dt_gt)
          rv_output$dt_gt <- dt_gt
        }
      })

      #### Sequence data ####

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
          rv_output$dt_seq <- dt_seq
        }
      })

      #### Kit-setting ####

      # Initial data.table
      react_dt_kit_set0 <- reactive({
        # Extract autosomal markers
        kit <- input$kit_par
        dt_kit <- fread(paste0(path_pack, "/extdata/kit/", kit, ".csv"))
        rv_output$dt_kit <- dt_kit <- change_data_type_kit(dt_kit)
        dt_kit_auto <- extract_kit_auto(dt_kit)
        markers_par <- unique(dt_kit_auto[, Marker])
        dyes_par <- extract_dye(dt_kit_auto, markers_par)

        # Make a data.table for kit-setting
        dt_kit_set0 <- data.table(Marker = markers_par, Dye = dyes_par)
        return(dt_kit_set0)
      })

      # Main data.table
      react_dt_kit_set1 <- reactive({
        dt_kit_set1 <- react_dt_kit_set0()
        dt_kit_set1$AT <- 50
        dt_kit_set1$Model_BSR <- "Best"
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
        dt_kit_set2$select_bsr <- selectColumn(n_mk, 4, choices = c("Best", "Allele", "LUS", "Multi-sequence", "Multiple loci together", "Not consider"), select_default = rep("Best", n_mk), width = "90%")
        dt_kit_set2$select_fsr <- selectColumn(n_mk, 5, choices = c("Best", "Allele", "LUS", "Multi-sequence", "Multiple loci together", "Not consider"), select_default = rep("Multiple loci together", n_mk), width = "90%")
        dt_kit_set2$select_dsr <- selectColumn(n_mk, 6, choices = c("Best", "Allele", "LUS", "Multi-sequence", "Multiple loci together", "Not consider"), select_default = rep("Multiple loci together", n_mk), width = "90%")
        dt_kit_set2$select_m2sr <- selectColumn(n_mk, 7, choices = c("Uniform", "Multiple loci together", "Not consider"), select_default = rep("Not consider", n_mk), width = "90%")
        return(dt_kit_set2)
      })

      # Select kit
      output$kit_par <- renderUI({
        kit_names <- list.files(paste0(path_pack, "/extdata/kit"))
        kit_names <- gsub(".csv", "", kit_names)
        selectInput(session$ns("kit_par"), label = h4("Kit"), choices = kit_names, selected = kit_names[1], multiple = FALSE, selectize = FALSE)
      })

      # Output data.table
      observeEvent(input$kit_par, ignoreInit = TRUE, {
        rv_output$dt_kit_set <- react_dt_kit_set1()
        dt_kit_set2 <- react_dt_kit_set2()

        output$dt_kit_set <- renderDataTable(server = FALSE, {
          datatable(
            dt_kit_set2,
            colnames = c("Marker", "Dye", "AT", "BSR model", "FSR model", "DSR model", "M2SR model"),
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_kit_set2), ordering = FALSE, autowidth = TRUE),
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
      rv_advanced$pu_remove <- TRUE
      rv_advanced$pu_dist <- 1
      rv_advanced$at_sex <- 100
      rv_advanced$min_n_motif_lus <- 2
      rv_advanced$min_n_motif_multi <- 3

      # Make the setting pop-up window
      observeEvent(input$act_advance_par, {
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
                   numericInput(session$ns("at_sex"), label = "AT of the sex-chromosomal markers", value = rv_advanced$at_sex),
                   br(),
                   h4("Stutter model"),
                   numericInput(session$ns("min_n_motif_lus"), label = "Minimum number of repeat motifs for considering the LUS model", value = rv_advanced$min_n_motif_lus),
                   numericInput(session$ns("min_n_motif_multi"), label = "Minimum number of repeat motifs for considering the multi-sequence model", value = rv_advanced$min_n_motif_multi)
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

        observeEvent(input$min_n_motif_lus, {
          min_n_motif_lus <- input$min_n_motif_lus
          if(!is.integer(min_n_motif_lus)){
            showFeedbackDanger(inputId = "min_n_motif_lus", text = "A positive integer is required!")
          }else if(min_n_motif_lus <= 0){
            showFeedbackDanger(inputId = "min_n_motif_lus", text = "A positive integer is required!")
          }else{
            hideFeedback("min_n_motif_lus")
          }
        })

        observeEvent(input$min_n_motif_multi, {
          min_n_motif_multi <- input$min_n_motif_multi
          if(!is.integer(min_n_motif_multi)){
            showFeedbackDanger(inputId = "min_n_motif_multi", text = "A positive integer is required!")
          }else if(min_n_motif_multi <= 0){
            showFeedbackDanger(inputId = "min_n_motif_multi", text = "A positive integer is required!")
          }else{
            hideFeedback("min_n_motif_multi")
          }
        })
      })

      # Save the setting
      observeEvent(input$act_advance_save, {
        req(input$warmup)
        req(input$sampling)
        req(input$chains)
        req(input$pu_dist)
        req(input$at_sex)
        req(input$min_n_motif_lus)
        req(input$min_n_motif_multi)

        isolate(warmup <- input$warmup)
        isolate(sampling <- input$sampling)
        isolate(chains <- input$chains)
        isolate(pu_remove <- input$pu_remove)
        isolate(pu_dist <- input$pu_dist)
        isolate(at_sex <- input$at_sex)
        isolate(min_n_motif_lus <- input$min_n_motif_lus)
        isolate(min_n_motif_multi <- input$min_n_motif_multi)

        if(all(is.integer(warmup), is.integer(sampling), is.integer(chains), is.numeric(pu_dist), is.numeric(at_sex), is.integer(min_n_motif_lus), is.integer(min_n_motif_multi))){
          if(all(warmup > 0, sampling > 0, chains > 0, pu_dist > 0, at_sex > 0, min_n_motif_lus > 0, min_n_motif_multi > 0))

          rv_advanced$warmup <- warmup
          rv_advanced$sampling <- sampling
          rv_advanced$chains <- chains
          rv_advanced$pu_remove <- pu_remove
          rv_advanced$pu_dist <- pu_dist
          rv_advanced$at_sex <- at_sex
          rv_advanced$min_n_motif_lus <- min_n_motif_lus
          rv_advanced$min_n_motif_multi <- min_n_motif_multi

          removeModal()
        }
      })

      #### Analyze ####

      observeEvent(input$act_analyze_est_par, {
        # Fix data.table
        dt_exp <- rv_output$dt_exp
        dt_gt <- rv_output$dt_gt
        dt_seq <- rv_output$dt_seq
        dt_kit <- rv_output$dt_kit
        dt_kit_set <- rv_output$dt_kit_set

        # Fix parameters
        warmup <- rv_advanced$warmup
        sampling <- rv_advanced$sampling
        chains <- rv_advanced$chains
        pu_remove <- rv_advanced$pu_remove
        pu_dist <- rv_advanced$pu_dist
        at_sex <- rv_advanced$at_sex
        min_n_motif_lus <- rv_advanced$min_n_motif_lus
        min_n_motif_multi <- rv_advanced$min_n_motif_multi
        kit_par <- input$kit_par

        # Check input files and setting
        error_message <- check_est_par(dt_exp, dt_gt, dt_seq, dt_kit_set)
        if(error_message != ""){
          showModal(modalDialog(title = "Error", HTML(error_message), easyClose = TRUE, footer = NULL))
        }else{
          # Load the MCMC models
          model_ae <- stan_package_model(name = "model_ae", package = "KongohPlus")
          model_sr_gen <- stan_package_model(name = "model_sr_gen", package = "KongohPlus")
          model_sr_uni <- stan_package_model(name = "model_sr_uni", package = "KongohPlus")
          model_gamma <- stan_package_model(name = "model_gamma", package = "KongohPlus")

          results_est_par <- estimate_parameters(dt_exp, dt_gt, dt_seq, dt_kit, dt_kit_set, kit_par,
                                                 model_ae, model_sr_gen, model_sr_uni, model_gamma,
                                                 warmup, sampling, chains,
                                                 pu_remove, pu_dist, at_sex, min_n_motif_lus, min_n_motif_multi)
          rv_output$results_est_par <- results_est_par

          updateTabsetPanel(session = session_param_top, inputId = "tabset_est_par", selected = "Result")
        }
      })

      return(rv_output)
    }
  )
}
