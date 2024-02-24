#' lr_setting_ui
#'
#' @description The function to create the ui module for setting conditions of the likelihood ratio
lr_setting_ui <- function(id){
  ns <- NS(id)

  tabPanel("Likelihood ratio",
           br(),
           sidebarLayout(
             sidebarPanel(width = 3,
                          uiOutput(ns("use_deconvo")),
                          numericInput(ns("noc_lr"), label = h4("Number of contributors"), value = 1),
                          uiOutput(ns("hp")),
                          uiOutput(ns("hd")),
                          br(),
                          uiOutput(ns("anameth_lr")),
                          br(),
                          actionButton(ns("act_advance_lr"), label = "Advanced setting"),
                          br(),
                          actionButton(ns("act_analyze_lr"), label = "Analyze", class = "btn btn-primary btn-lg")
             ),
             mainPanel(width = 9,
                       fluidRow(
                         column(4,
                                dataTableOutput(ns("dt_at_lr")),
                                br()
                         ),
                         column(8,
                                h4("Kit")
                         )
                       )
             )
           )
  )
}

#' lr_setting_server
#'
#' @description The function to create the server module for setting conditions of the likelihood ratio
lr_setting_server <- function(id, path_pack, rv_file, rv_deconvo){
  moduleServer(
    id,
    function(input, output, session){
      results_deconvo <- rv_deconvo$results

      rv_lr_setting <- reactiveValues()
      rv_lr_setting$refs <- NULL
      rv_lr_setting$dt_anameth <- NULL
      rv_lr_setting$warmup <- 1200
      rv_lr_setting$sampling <- 300
      rv_lr_setting$chains <- 4

      dt_ref <- rv_file$dt_ref
      if(!is.null(dt_ref)){
        rv_lr_setting$refs <- unique(dt_ref[, SampleName])
      }

      output$use_deconvo <- renderUI({
        if(is.null(results_deconvo)){
          disabled(checkboxInput(session$ns("use_deconvo"), label = "Use the result of deconvolution", value = FALSE))
        }else{
          checkboxInput(session$ns("use_deconvo"), label = "Use the result of deconvolution", value = TRUE)
        }
      })

      observeEvent(input$use_deconvo, {
        use_deconvo <- input$use_deconvo
        if(use_deconvo){
          disable("noc_lr")
          updateNumericInput(session, "noc_lr", value = results_deconvo$hnc)
          disable("act_advance_lr")
        }else{
          enable("noc_lr")
          enable("act_advance_lr")
        }
      })

      observeEvent(input$noc_lr, {
        noc_lr <- input$noc_lr
        if(!is.integer(noc_lr) || noc_lr < 1 || noc_lr > 4){
          showFeedbackDanger(inputId = "noc_lr", text = "Only 1, 2, 3, and 4 is allowed.")
          disable("act_analyze_lr")
        }else{
          hideFeedback("noc_lr")
          enable("act_analyze_lr")
        }
      })

      output$hp <- renderUI({
        checkboxGroupInput(session$ns("hp"), label = h4("Prosecutor hypothesis"), choices = rv_lr_setting$refs)
      })

      output$hd <- renderUI({
        checkboxGroupInput(session$ns("hd"), label = h4("Defense hypothesis"), choices = rv_lr_setting$refs)
      })

      observeEvent(input$act_advance_lr, {
        showModal(modalDialog(
          title = "Advanced setting",
          fluidRow(
            column(6,
                   h4("MCMC conditions"),
                   numericInput(session$ns("warmup"), label = "Warmup iterations", value = rv_lr_setting$warmup),
                   numericInput(session$ns("sampling"), label = "Sampling iterations", value = rv_lr_setting$sampling),
                   numericInput(session$ns("chains"), label = "Number of chains", value = rv_lr_setting$chains)
            ),
            column(6)
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

      observeEvent(input$act_advance_save, {
        req(input$warmup)
        req(input$sampling)
        req(input$chains)

        isolate(warmup <- input$warmup)
        isolate(sampling <- input$sampling)
        isolate(chains <- input$chains)

        if(all(is.integer(warmup), is.integer(sampling), is.integer(chains))){
          rv_lr_setting$warmup <- warmup
          rv_lr_setting$sampling <- sampling
          rv_lr_setting$chains <- chains
          removeModal()
        }
      })

      anameth <- list.files(paste0(path_pack, "/extdata/anameth"))

      output$anameth_lr <- renderUI({
        anameth_all <- list.files(paste0(path_pack, "/extdata/anameth"))
        anameth_all <- gsub(".{4}$", "", anameth_all)
        selectInput(session$ns("anameth_lr"), label = h4("Analysis method"), choices = anameth_all, selected = anameth_all[1])
      })

      observeEvent(input$anameth_lr, {
        # Load the selected analysis method
        anameth <- input$anameth_lr
        rv_lr_setting$dt_anameth <- dt_anameth <- fread(paste0(path_pack, "/extdata/anameth/", anameth, ".csv"))

        # Extract data per column
        folders <- dt_anameth[, Folder]
        names <- dt_anameth[, Name]

        rv_lr_setting$dt_at <- dt_at <- fread(paste0(path_pack, "/extdata/at/", names[folders == "at"], ".csv"))

        # Display analytical thresholds
        output$dt_at_lr <- renderDataTable(server = FALSE, {
          datatable(
            dt_at,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_at), ordering = FALSE, autowidth = TRUE),
            rownames = FALSE
          )
        })
      })
    }
  )
}
