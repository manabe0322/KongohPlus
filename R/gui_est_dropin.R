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
                        column(9, uiOutput(ns("name_kit")))
                      ),
                      br(),
                      numericInput(ns("at_dropin"), label = h4("Analytical threshold"), value = 25),
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
      rv_est_dropin$name_kit <- NULL
      rv_est_dropin$n_conv <- 10
      rv_est_dropin$n_max_init <- 100
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
          rv_est_dropin$name_kit <- str_replace(basename(path_kit), ".{4}$", "")
        }
      })

      output$name_kit <- renderUI({
        disabled(textInput(session$ns("name_kit"), label = "", value = rv_est_dropin$name_kit))
      })

      #### Advanced setting ####

      observeEvent(input$act_advance_dropin, {
        showModal(modalDialog(
          title = "Advanced setting",
          h4("MLE conditions"),
          numericInput(session$ns("n_conv"), label = "Number of successful convergences", value = rv_est_dropin$n_conv),
          numericInput(session$ns("n_max_init"), label = "Maximum number of initializations", value = rv_est_dropin$n_max_init),
          footer = tagList(
            actionButton(session$ns("act_advance_dropin_save"), "Save"),
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
          enable("act_advance_dropin_save")
        }else{
          disable("act_advance_dropin_save")
        }
      })

      observeEvent(input$act_advance_dropin_save, {
        rv_est_dropin$n_conv <- input$n_conv
        rv_est_dropin$n_max_init <- input$n_max_init
        removeModal()
      })

      # Start estimating dropin parameters
      observeEvent(input$act_analyze_est_dropin, {
        # Fix data
        dt_ntc <- rv_est_dropin$dt_ntc
        fn_ntc <- rv_est_dropin$fn_ntc
        name_kit <- rv_est_dropin$name_kit
        at_dropin <- input$at_dropin
        n_conv <- rv_est_dropin$n_conv
        n_max_init <- rv_est_dropin$n_max_init
        dt_kit <- rv_est_dropin$dt_kit

        # Check input files and setting
        error_message <- check_est_dropin(dt_ntc, dt_kit)
        if(error_message != ""){
          showModal(modalDialog(title = "Error", HTML(error_message), easyClose = TRUE, footer = NULL))
        }else{
          # Load the stan model
          model_dropin <- stan_package_model(name = "model_dropin", package = "KongohPlus")

          data_est_dropin <- estimate_dropin_params(dt_ntc, dt_kit, at_dropin, model_dropin, n_conv, n_max_init)
          data_est_dropin$fn_ntc <- fn_ntc
          data_est_dropin$name_kit <- name_kit
          data_est_dropin$at_dropin <- at_dropin
          data_est_dropin$n_conv <- n_conv
          data_est_dropin$n_max_init <- n_max_init
          rv_est_dropin$data_est_dropin <- data_est_dropin
        }
      })

      observe({
        req(rv_est_dropin)
        data_est_dropin <- rv_est_dropin$data_est_dropin
        if(!is.null(data_est_dropin)){
          fn_ntc <- data_est_dropin$fn_ntc
          name_kit <- data_est_dropin$name_kit
          at_dropin <- data_est_dropin$at_dropin
          n_conv <- data_est_dropin$n_conv
          n_max_init <- data_est_dropin$n_max_init
          enable("act_saveas_est_dropin_1")

          dt_dropin <- data.table(Parameter = c("Number of samples", "Number of markers", "Number of drop-in peaks", "Frequency of drop-in", "Shape parameter", "Rate parameter"),
                                  Value = c(data_est_dropin$n_sample, data_est_dropin$n_mk, data_est_dropin$n_dropin, data_est_dropin$freq_dropin, data_est_dropin$shape_di, data_est_dropin$rate_di))
          traceplot_shape_di <- data_est_dropin$traceplot_shape_di
          traceplot_rate_di <- data_est_dropin$traceplot_rate_di

          output$dt_dropin <- renderDataTable({
            datatable(
              dt_dropin,
              selection = "none",
              options = list(dom = "t", iDisplayLength = nrow(dt_dropin), ordering = FALSE, autowidth = TRUE),
              rownames = FALSE
            )
          })

          rv_dropin_file_names <- reactiveVal(NULL)

          observeEvent(input$act_saveas_est_dropin_1, {
            dropin_file_names <- list.files(paste0(path_pack, "/extdata/dropin"))
            rv_dropin_file_names(gsub(".{6}$", "", dropin_file_names))

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
            if(nchar(fn_est_dropin) == 0){
              showFeedbackDanger(inputId = "fn_est_dropin", text = "Required")
              disable("act_saveas_est_dropin_2")
            }else if(is.element(fn_est_dropin, rv_dropin_file_names())){
              showFeedbackDanger(inputId = "fn_est_dropin", text = "This name is already registered.")
              disable("act_saveas_est_dropin_2")
            }else{
              hideFeedback("fn_est_dropin")
              enable("act_saveas_est_dropin_2")
            }
          })

          observeEvent(input$act_saveas_est_dropin_2, {
            fn_est_dropin <- input$fn_est_dropin

            save(list = c("fn_ntc", "name_kit", "at_dropin", "n_conv", "n_max_init", "dt_dropin"),
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
