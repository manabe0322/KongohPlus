#' est_dropin_ui
#'
#' @description The function to create the ui module for estimating drop-in parameters
est_dropin_ui <- function(id){
  ns <- NS(id)

  tabPanel("6. Estimate drop-in parameters",
           titlePanel("6. Estimate drop-in parameters"),
           br(),
           fluidRow(
             column(3,
                    wellPanel(
                      fileInput(ns("ntc_data"), label = h4("Non template control data"), accept = ".csv", buttonLabel = "Load"),
                      br(),
                      h4("Kit"),
                      fluidRow(
                        column(3, br(), shinyFilesButton(ns("load_kit"), label = "Load", title = "Kit", multiple = FALSE)),
                        column(9, uiOutput(ns("name_kit_dropin")))
                      ),
                      br(),
                      numericInput(ns("at_dropin"), label = h4("Analytical threshold"), value = 35),
                      br(),
                      actionButton(ns("act_advance_dropin"), "Advanced setting"),
                      br(),
                      br(),
                      br(),
                      actionButton(ns("act_analyze_est_dropin"), label = "Analyze", class = "btn btn-primary btn-lg"),
                      br(),
                      br(),
                      disabled(actionButton(ns("act_saveas_est_dropin_1"), label = "Save as", class = "btn btn-primary btn-lg"))
                    )
             ),
             column(5,
                    wellPanel(
                      dataTableOutput(ns("dt_dropin"))
                    )
             ),
             column(4)
           )
  )
}

#' est_dropin_server
#'
#' @description The function to create the server module for estimating drop-in parameters
est_dropin_server <- function(id, path_pack){
  moduleServer(
    id,
    function(input, output, session){
      rv_est_dropin <- reactiveValues()
      rv_est_dropin$dt_ntc <- NULL
      rv_est_dropin$fn_ntc <- NULL
      rv_est_dropin$dt_kit <- NULL
      rv_est_dropin$name_kit_dropin <- NULL
      rv_est_dropin$min_shape_di <- 0.01
      rv_est_dropin$max_shape_di <- 20
      rv_est_dropin$min_rate_di <- 0.01
      rv_est_dropin$max_rate_di <- 20
      rv_est_dropin$stop_iter <- 500
      rv_est_dropin$max_iter <- 5000
      rv_est_dropin$num_converged <- 3
      rv_est_dropin$max_mle_runs <- 10
      rv_est_dropin$bool_correct_mle_shape_di <- TRUE
      rv_est_dropin$bool_correct_mle_rate_di <- TRUE
      rv_est_dropin$bool_correct_advanced_1 <- TRUE
      rv_est_dropin$bool_correct_advanced_2 <- TRUE
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
          rv_est_dropin$name_kit_dropin <- str_replace(basename(path_kit), ".{4}$", "")
        }
      })

      output$name_kit_dropin <- renderUI({
        disabled(textInput(session$ns("name_kit_dropin"), label = "", value = rv_est_dropin$name_kit_dropin))
      })

      #### Advanced setting ####

      observeEvent(input$act_advance_dropin, {
        showModal(modalDialog(
          title = "Advanced setting",
          h4("MLE conditions"),
          fluidRow(
            column(6,
                   numericInput(session$ns("min_shape_di"), label = "Minimum shape parameter", value = rv_est_dropin$min_shape_di),
                   numericInput(session$ns("min_rate_di"), label = "Minimum rate parameter", value = rv_est_dropin$min_rate_di),
                   numericInput(session$ns("stop_iter"), label = "Number of iterations without improvement before stopping one MLE run", value = rv_est_dropin$stop_iter),
                   numericInput(session$ns("num_converged"), label = "The number of converged solutions required to terminate the MLE", rv_est_dropin$num_converged)
            ),
            column(6,
                   numericInput(session$ns("max_shape_di"), label = "Maximum shape parameter", value = rv_est_dropin$max_shape_di),
                   numericInput(session$ns("max_rate_di"), label = "Maximum rate parameter", value = rv_est_dropin$max_rate_di),
                   numericInput(session$ns("max_iter"), label = "Number of iterations", value = rv_est_dropin$max_iter),
                   numericInput(session$ns("max_mle_runs"), label = "The maximum number of MLE attempts with different initial values", value = rv_est_dropin$max_mle_runs)
            )
          ),
          footer = tagList(
            actionButton(session$ns("act_advance_dropin_save"), "Save"),
            modalButton("Cancel")
          ),
          size = "l"
        ))
      })

      observeEvent(input$min_shape_di | input$max_shape_di, {
        min_shape_di <- input$min_shape_di
        max_shape_di <- input$max_shape_di

        hideFeedback("min_shape_di")
        hideFeedback("max_shape_di")
        rv_est_dropin$bool_correct_mle_shape_di <- TRUE

        if(!is.numeric(min_shape_di) || min_shape_di <= 0){
          showFeedbackDanger(inputId = "min_shape_di", text = "Value must be positive.")
          rv_est_dropin$bool_correct_mle_shape_di <- FALSE
        }

        if(!is.numeric(max_shape_di) || max_shape_di <= 0){
          showFeedbackDanger(inputId = "max_shape_di", text = "Value must be positive.")
          rv_est_dropin$bool_correct_mle_shape_di <- FALSE
        }

        if(is.numeric(min_shape_di) && is.numeric(max_shape_di) && max_shape_di < min_shape_di){
          showFeedbackDanger(inputId = "min_shape_di", text = "Minimum value cannot be greater than maximum value.")
          showFeedbackDanger(inputId = "max_shape_di", text = "Minimum value cannot be greater than maximum value.")
          rv_est_dropin$bool_correct_mle_shape_di <- FALSE
        }
      })

      observeEvent(input$min_rate_di | input$max_rate_di, {
        min_rate_di <- input$min_rate_di
        max_rate_di <- input$max_rate_di

        hideFeedback("min_rate_di")
        hideFeedback("max_rate_di")
        rv_est_dropin$bool_correct_mle_rate_di <- TRUE

        if(!is.numeric(min_rate_di) || min_rate_di <= 0){
          showFeedbackDanger(inputId = "min_rate_di", text = "Value must be positive.")
          rv_est_dropin$bool_correct_mle_rate_di <- FALSE
        }

        if(!is.numeric(max_rate_di) || max_rate_di <= 0){
          showFeedbackDanger(inputId = "max_rate_di", text = "Value must be positive.")
          rv_est_dropin$bool_correct_mle_rate_di <- FALSE
        }

        if(is.numeric(min_rate_di) && is.numeric(max_rate_di) && max_rate_di < min_rate_di){
          showFeedbackDanger(inputId = "min_rate_di", text = "Minimum value cannot be greater than maximum value.")
          showFeedbackDanger(inputId = "max_rate_di", text = "Minimum value cannot be greater than maximum value.")
          rv_est_dropin$bool_correct_mle_rate_di <- FALSE
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
          showFeedbackDanger(inputId = "stop_iter", text = "Number of iterations without improvement must not exceed maximum number of iterations.")
          showFeedbackDanger(inputId = "max_iter", text = "Number of iterations without improvement must not exceed maximum number of iterations.")
          bool_correct_advanced_1 <- FALSE
        }

        rv_est_dropin$bool_correct_advanced_1 <- bool_correct_advanced_1
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
          showFeedbackDanger(inputId = "num_converged", text = "The number of converged solutions must not exceed the maximum number of MLE attempts.")
          showFeedbackDanger(inputId = "max_mle_runs", text = "The number of converged solutions must not exceed the maximum number of MLE attempts.")
          bool_correct_advanced_2 <- FALSE
        }

        rv_est_dropin$bool_correct_advanced_2 <- bool_correct_advanced_2
      })

      observe({
        req(rv_est_dropin)

        if(all(rv_est_dropin$bool_correct_mle_shape_di,
               rv_est_dropin$bool_correct_mle_rate_di,
               rv_est_dropin$bool_correct_advanced_1,
               rv_est_dropin$bool_correct_advanced_2)){
          enable("act_advance_dropin_save")
        }else{
          disable("act_advance_dropin_save")
        }
      })

      observeEvent(input$act_advance_dropin_save, {
        rv_est_dropin$min_shape_di <- input$min_shape_di
        rv_est_dropin$max_shape_di <- input$max_shape_di
        rv_est_dropin$min_rate_di <- input$min_rate_di
        rv_est_dropin$max_rate_di <- input$max_rate_di
        rv_est_dropin$stop_iter <- input$stop_iter
        rv_est_dropin$max_iter <- input$max_iter
        rv_est_dropin$num_converged <- input$num_converged
        rv_est_dropin$max_mle_runs <- input$max_mle_runs
        removeModal()
      })

      # Start estimating dropin parameters
      observeEvent(input$act_analyze_est_dropin, {
        # Fix data
        dt_ntc <- rv_est_dropin$dt_ntc
        fn_ntc <- rv_est_dropin$fn_ntc
        name_kit_dropin <- rv_est_dropin$name_kit_dropin
        at_dropin <- input$at_dropin
        min_shape_di <- rv_est_dropin$min_shape_di
        max_shape_di <- rv_est_dropin$max_shape_di
        min_rate_di <- rv_est_dropin$min_rate_di
        max_rate_di <- rv_est_dropin$max_rate_di
        stop_iter <- rv_est_dropin$stop_iter
        max_iter <- rv_est_dropin$max_iter
        num_converged <- rv_est_dropin$num_converged
        max_mle_runs <- rv_est_dropin$max_mle_runs
        dt_kit <- rv_est_dropin$dt_kit

        # Check input files and setting
        error_message <- check_est_dropin(dt_ntc, dt_kit)
        if(error_message != ""){
          showModal(modalDialog(title = "Error", HTML(error_message), easyClose = TRUE, footer = NULL))
        }else{
          data_est_dropin <- estimate_dropin_params(dt_ntc, dt_kit, at_dropin, min_shape_di, max_shape_di, min_rate_di, max_rate_di, stop_iter, max_iter, num_converged, max_mle_runs)
          data_est_dropin$fn_ntc <- fn_ntc
          data_est_dropin$name_kit_dropin <- name_kit_dropin
          data_est_dropin$at_dropin <- at_dropin
          data_est_dropin$min_shape_di <- min_shape_di
          data_est_dropin$max_shape_di <- max_shape_di
          data_est_dropin$min_rate_di <- min_rate_di
          data_est_dropin$max_rate_di <- max_rate_di
          data_est_dropin$stop_iter <- stop_iter
          data_est_dropin$max_iter <- max_iter
          data_est_dropin$num_converged <- num_converged
          data_est_dropin$max_mle_runs <- max_mle_runs
          rv_est_dropin$data_est_dropin <- data_est_dropin
        }
      })

      observe({
        req(rv_est_dropin)
        data_est_dropin <- rv_est_dropin$data_est_dropin
        if(!is.null(data_est_dropin)){
          fn_ntc <- data_est_dropin$fn_ntc
          name_kit_dropin <- data_est_dropin$name_kit_dropin
          at_dropin <- data_est_dropin$at_dropin
          min_shape_di <- data_est_dropin$min_shape_di
          max_shape_di <- data_est_dropin$max_shape_di
          min_rate_di <- data_est_dropin$min_rate_di
          max_rate_di <- data_est_dropin$max_rate_di
          stop_iter <- data_est_dropin$stop_iter
          max_iter <- data_est_dropin$max_iter
          num_converged <- data_est_dropin$num_converged
          max_mle_runs <- data_est_dropin$max_mle_runs
          n_sample_di <- data_est_dropin$n_sample
          n_mk_di <- data_est_dropin$n_mk
          n_dropin <- data_est_dropin$n_dropin
          freq_dropin <- data_est_dropin$freq_dropin
          shape_di <- data_est_dropin$shape_di
          rate_di <- data_est_dropin$rate_di
          init_list_di <- data_est_dropin$init_list_di
          opt_list_di <- data_est_dropin$opt_list_di
          enable("act_saveas_est_dropin_1")

          dt_dropin <- data.table(Parameter = c("Number of samples", "Number of markers", "Number of drop-in peaks", "Frequency", "Shape parameter", "Rate parameter"),
                                  Value = c(n_sample_di, n_mk_di, n_dropin, freq_dropin, shape_di, rate_di))

          output$dt_dropin <- renderDataTable({
            datatable(
              dt_dropin,
              selection = "none",
              options = list(dom = "t", iDisplayLength = nrow(dt_dropin), ordering = FALSE, autowidth = TRUE),
              rownames = FALSE
            )
          })

          observeEvent(input$act_saveas_est_dropin_1, {
            showModal(modalDialog(
              title = "Save drop-in data",
              textInput(session$ns("fn_est_dropin"), label = "Enter the file name"),
              footer = tagList(
                disabled(actionButton(session$ns("act_saveas_est_dropin_2"), "Save")),
                modalButton("Close")
              )
            ))
          })

          observeEvent(input$fn_est_dropin, {
            fn_est_dropin <- input$fn_est_dropin
            dropin_file_names <- list.files(paste0(path_pack, "/extdata/dropin"))
            if(is.element(paste0(fn_est_dropin, ".RData"), dropin_file_names)){
              showFeedbackDanger(inputId = "fn_est_dropin", text = "A file with the same name already exists.")
              disable("act_saveas_est_dropin_2")
            }else if(nchar(fn_est_dropin) == 0){
              hideFeedback("fn_est_dropin")
              disable("act_saveas_est_dropin_2")
            }else{
              hideFeedback("fn_est_dropin")
              enable("act_saveas_est_dropin_2")
            }
          })

          observeEvent(input$act_saveas_est_dropin_2, {
            fn_est_dropin <- input$fn_est_dropin

            save(list = c("fn_ntc", "name_kit_dropin", "at_dropin", "min_shape_di", "max_shape_di", "min_rate_di", "max_rate_di", "stop_iter", "max_iter", "num_converged", "max_mle_runs",
                          "n_sample_di", "n_mk_di", "n_dropin", "freq_dropin", "shape_di", "rate_di", "init_list_di", "opt_list_di"),
                 file = paste0(path_pack, "/extdata/dropin/", fn_est_dropin, ".RData"))

            removeModal()
          })
        }else{
          disable("act_saveas_est_dropin_1")
        }
      })
    }
  )
}
