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
      rv_est_dropin$min_shape_di <- 0.01
      rv_est_dropin$max_shape_di <- 20
      rv_est_dropin$min_rate_di <- 0.01
      rv_est_dropin$max_rate_di <- 20
      rv_est_dropin$initial_point_factor <- 5
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
          fluidRow(
            column(6,
                   numericInput(session$ns("min_shape_di"), label = "Minimum shape parameter", value = rv_est_dropin$min_shape_di),
                   numericInput(session$ns("min_rate_di"), label = "Minimum rate parameter", value = rv_est_dropin$min_rate_di),
                   numericInput(session$ns("initial_point_factor"), label = "Initial point factor", value = rv_est_dropin$initial_point_factor)
            ),
            column(6,
                   numericInput(session$ns("max_shape_di"), label = "Maximum shape parameter", value = rv_est_dropin$max_shape_di),
                   numericInput(session$ns("max_rate_di"), label = "Maximum rate parameter", value = rv_est_dropin$max_rate_di)
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

      observeEvent(input$initial_point_factor, {
        initial_point_factor <- input$initial_point_factor

        hideFeedback("initial_point_factor")
        rv_est_dropin$bool_correct_mle_initial_point_factor <- TRUE

        if(!is.integer(initial_point_factor) || initial_point_factor <= 0){
          showFeedbackDanger(inputId = "initial_point_factor", text = "Value must be an integer greater than zero.")
          rv_est_dropin$bool_correct_mle_initial_point_factor <- FALSE
        }
      })

      observe({
        req(rv_est_dropin)

        if(all(rv_est_dropin$bool_correct_mle_shape_di,
               rv_est_dropin$bool_correct_mle_rate_di,
               rv_est_dropin$bool_correct_mle_initial_point_factor)){
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
        rv_est_dropin$initial_point_factor <- input$initial_point_factor
        removeModal()
      })

      # Start estimating dropin parameters
      observeEvent(input$act_analyze_est_dropin, {
        # Fix data
        dt_ntc <- rv_est_dropin$dt_ntc
        fn_ntc <- rv_est_dropin$fn_ntc
        name_kit <- rv_est_dropin$name_kit
        at_dropin <- input$at_dropin
        min_shape_di <- rv_est_dropin$min_shape_di
        max_shape_di <- rv_est_dropin$max_shape_di
        min_rate_di <- rv_est_dropin$min_rate_di
        max_rate_di <- rv_est_dropin$max_rate_di
        initial_point_factor <- rv_est_dropin$initial_point_factor
        dt_kit <- rv_est_dropin$dt_kit

        # Check input files and setting
        error_message <- check_est_dropin(dt_ntc, dt_kit)
        if(error_message != ""){
          showModal(modalDialog(title = "Error", HTML(error_message), easyClose = TRUE, footer = NULL))
        }else{
          # Load the stan model
          model_dropin <- stan_package_model(name = "model_dropin", package = "KongohPlus")

          data_est_dropin <- estimate_dropin_params(dt_ntc, dt_kit, at_dropin, model_dropin, min_shape_di, max_shape_di, min_rate_di, max_rate_di, initial_point_factor)
          data_est_dropin$fn_ntc <- fn_ntc
          data_est_dropin$name_kit <- name_kit
          data_est_dropin$at_dropin <- at_dropin
          data_est_dropin$min_shape_di <- min_shape_di
          data_est_dropin$max_shape_di <- max_shape_di
          data_est_dropin$min_rate_di <- min_rate_di
          data_est_dropin$max_rate_di <- max_rate_di
          data_est_dropin$initial_point_factor <- initial_point_factor
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
          min_shape_di <- data_est_dropin$min_shape_di
          max_shape_di <- data_est_dropin$max_shape_di
          min_rate_di <- data_est_dropin$min_rate_di
          max_rate_di <- data_est_dropin$max_rate_di
          initial_point_factor <- data_est_dropin$initial_point_factor
          enable("act_saveas_est_dropin_1")

          dt_dropin <- data.table(Parameter = c("Number of samples", "Number of markers", "Number of drop-in peaks", "Frequency of drop-in", "Shape parameter", "Rate parameter"),
                                  Value = c(data_est_dropin$n_sample, data_est_dropin$n_mk, data_est_dropin$n_dropin, data_est_dropin$freq_dropin, data_est_dropin$shape_di, data_est_dropin$rate_di))

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

            save(list = c("fn_ntc", "name_kit", "at_dropin", "min_shape_di", "max_shape_di", "min_rate_di", "max_rate_di", "initial_point_factor", "dt_dropin"),
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
