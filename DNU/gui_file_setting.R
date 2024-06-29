#' file_setting_ui
#'
#' @description The function to create the ui module for files and settings
file_setting_ui <- function(id){
  ns <- NS(id)

  tabPanel("Files and Settings",
           fluidRow(
             column(4,
                    wellPanel(
                      h2("Load"),
                      br(),
                      fileInput(ns("file_csp"), label = "Crime stain profile", accept = ".csv", buttonLabel = "Load"),
                      fileInput(ns("file_ref"), label = "Reference profile", accept = ".csv", buttonLabel = "Load"),
                      fileInput(ns("file_af"), label = "Allele frequency", accept = ".csv", buttonLabel = "Load")
                    ),
                    wellPanel(
                      h2("Setting"),
                      br(),
                      uiOutput(ns("anameth")),
                      br(),
                      actionButton(ns("act_advance"), label = "Advanced setting")
                    )
             ),
             column(4,
                    wellPanel(
                      h2("Deconvolution"),
                      br(),
                      numericInput(ns("noc_deconvo"), label = "Number of contributors", value = 1),
                      br(),
                      actionButton(ns("act_analyze_deconvo"), label = "Analyze", class = "btn btn-primary btn-lg")
                    )
             ),
             column(4,
                    wellPanel(
                      h2("Likelihood ratio"),
                      uiOutput(ns("use_deconvo")),
                      uiOutput(ns("hp")),
                      numericInput(ns("noc_lr_hp"), label = "Number of contributors", value = 1),
                      uiOutput(ns("hd")),
                      numericInput(ns("noc_lr_hd"), label = "Number of contributors", value = 1),
                      br(),
                      actionButton(ns("act_analyze_lr"), label = "Analyze", class = "btn btn-primary btn-lg")
                    )
             )
           )
  )
}

#' file_setting_server
#'
#' @description The function to create the server module for files and settings
#' @param path_pack Package path
file_setting_server <- function(id, path_pack){
  moduleServer(
    id,
    function(input, output, session){
      rv_file_setting <- reactiveValues()
      rv_file_setting$dt_csp <- NULL
      rv_file_setting$dt_ref <- NULL
      rv_file_setting$dt_af <- NULL

      rv_file_setting$anameth <- NULL
      rv_file_setting$dt_anameth <- NULL
      rv_file_setting$dt_kit <- NULL
      rv_file_setting$dt_at <- NULL
      rv_file_setting$dt_param <- NULL
      rv_file_setting$dt_alcor <- NULL
      rv_file_setting$gtc_fltr <- NULL
      rv_file_setting$gtc_prob <- NULL
      rv_file_setting$dt_ae <- NULL
      rv_file_setting$dt_bsr <- NULL
      rv_file_setting$dt_fsr <- NULL
      rv_file_setting$dt_dsr <- NULL
      rv_file_setting$dt_m2sr <- NULL
      rv_file_setting$dt_gamma <- NULL

      rv_file_setting$warmup <- 1200
      rv_file_setting$sampling <- 300
      rv_file_setting$chains <- 4

      # Load the crime stain profile
      observeEvent(input$file_csp, {
        file_input <- input$file_csp
        if(!is.null(file_input)){
          dt_csp <- fread(file_input$datapath)
          dt_csp <- change_data_type_csp(dt_csp)
          rv_file_setting$dt_csp <- dt_csp
        }
      })

      # Load the reference profile
      observeEvent(input$file_ref, {
        file_input <- input$file_ref
        if(!is.null(file_input)){
          dt_ref <- fread(file_input$datapath)
          dt_ref <- change_data_type_ref(dt_ref)
          rv_file_setting$dt_ref <- dt_ref
        }
      })

      # Load the allele frequency
      observeEvent(input$file_af, {
        file_input <- input$file_af
        if(!is.null(file_input)){
          dt_af <- fread(file_input$datapath)
          dt_af <- change_data_type_af(dt_af)
          rv_file_setting$dt_af <- dt_af
        }
      })

      observeEvent(input$act_advance, {
        showModal(modalDialog(
          title = "Advanced setting",
          fluidRow(
            column(6,
                   h4("MCMC conditions"),
                   numericInput(session$ns("warmup"), label = "Warmup iterations", value = rv_file_setting$warmup),
                   numericInput(session$ns("sampling"), label = "Sampling iterations", value = rv_file_setting$sampling),
                   numericInput(session$ns("chains"), label = "Number of chains", value = rv_file_setting$chains)
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
          rv_file_setting$warmup <- warmup
          rv_file_setting$sampling <- sampling
          rv_file_setting$chains <- chains
          removeModal()
        }
      })

      output$anameth <- renderUI({
        anameth_all <- list.files(paste0(path_pack, "/extdata/anameth"))
        anameth_all <- gsub(".{4}$", "", anameth_all)
        selectInput(session$ns("anameth"), label = "Analysis method", choices = anameth_all, selected = anameth_all[1])
      })

      observeEvent(input$anameth, {
        # Load the selected analysis method
        anameth <- input$anameth
        rv_file_setting$anameth <- anameth
        rv_file_setting$dt_anameth <- dt_anameth <- fread(paste0(path_pack, "/extdata/anameth/", anameth, ".csv"))

        # Extract data per column
        folders <- dt_anameth[, Folder]
        names <- dt_anameth[, Name]

        # Load each data.table
        rv_file_setting$dt_kit <- fread(paste0(path_pack, "/extdata/kit/", names[folders == "kit"], ".csv"))
        rv_file_setting$dt_at <- dt_at <- fread(paste0(path_pack, "/extdata/at/", names[folders == "at"], ".csv"))
        rv_file_setting$dt_param <- fread(paste0(path_pack, "/extdata/param/", names[folders == "param"], ".csv"))
        rv_file_setting$dt_alcor <- fread(paste0(path_pack, "/extdata/alcor/", names[folders == "alcor"], ".csv"))
        rv_file_setting$dt_gtc_fltr <- dt_gtc_fltr <- fread(paste0(path_pack, "/extdata/gtc_fltr/", names[folders == "gtc_fltr"], ".csv"))
        rv_file_setting$dt_gtc_prob <- dt_gtc_prob <- fread(paste0(path_pack, "/extdata/gtc_prob/", names[folders == "gtc_prob"], ".csv"))

        rv_file_setting$dt_ae <- dt_param["AE"]
        rv_file_setting$dt_bsr <- dt_param["BSR"]
        rv_file_setting$dt_fsr <- dt_param["FSR"]
        rv_file_setting$dt_dsr <- dt_param["DSR"]
        rv_file_setting$dt_m2sr <- dt_param["M2SR"]
        rv_file_setting$dt_gamma <- dt_param["Gamma"]
      })

      observeEvent(input$noc_deconvo, {
        noc_deconvo <- input$noc_deconvo
        if(!is.integer(noc_deconvo) || noc_deconvo < 1 || noc_deconvo > 4){
          showFeedbackDanger(inputId = "noc_deconvo", text = "Only 1, 2, 3, and 4 is allowed.")
          disable("act_analyze_deconvo")
        }else{
          hideFeedback("noc_deconvo")
          enable("act_analyze_deconvo")
        }
      })

      observeEvent(input$act_analyze_deconvo, {

      })




      return(rv_file_setting)
    }
  )
}
