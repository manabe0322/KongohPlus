#' deconvo_setting_ui
#'
#' @description The function to create the ui module for setting conditions of deconvolution
deconvo_setting_ui <- function(id){
  ns <- NS(id)

  tabPanel("Setting",
           useWaiter(),
           br(),
           sidebarLayout(
             sidebarPanel(width = 3,
                          numericInput(ns("noc_deconvo"), label = h4("Number of contributors"), value = 1),
                          br(),
                          uiOutput(ns("anameth_deconvo")),
                          br(),
                          actionButton(ns("act_advance_deconvo"), label = "Advanced setting"),
                          br(),
                          br(),
                          actionButton(ns("act_analyze_deconvo"), label = "Analyze", class = "btn btn-primary btn-lg")
             ),
             mainPanel(width = 9,
                       tags$style(HTML("
                         .handsontable .htCore td {
                           vertical-align: middle !important;
                           text-align: center !important;
                         }"
                       )),
                       rHandsontableOutput(ns("dt_info_mk"))
             )
           )
  )
}

#' deconvo_setting_server
#'
#' @description The function to create the server module for setting conditions of deconvolution
deconvo_setting_server <- function(id, path_pack, rv_file){
  moduleServer(
    id,
    function(input, output, session){
      rv_deconvo_setting <- reactiveValues()
      rv_deconvo_setting$dt_csp <- NULL
      rv_deconvo_setting$dt_af <- NULL
      rv_deconvo_setting$fn_csp <- NULL
      rv_deconvo_setting$fn_af <- NULL
      rv_deconvo_setting$dt_anameth <- NULL
      rv_deconvo_setting$warmup <- 1200
      rv_deconvo_setting$sampling <- 300
      rv_deconvo_setting$chains <- 4

      rv_deconvo <- reactiveValues()
      rv_deconvo$data_deconvo <- NULL

      observe({
        req(rv_file)
        rv_deconvo_setting$dt_csp <- rv_file$dt_csp
        rv_deconvo_setting$dt_af <- rv_file$dt_af
        rv_deconvo_setting$fn_csp <- rv_file$fn_csp
        rv_deconvo_setting$fn_af <- rv_file$fn_af
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

      observeEvent(input$act_advance_deconvo, {
        showModal(modalDialog(
          title = "Advanced setting",
          fluidRow(
            column(6,
                   h4("MCMC conditions"),
                   numericInput(session$ns("warmup"), label = "Warmup iterations", value = rv_deconvo_setting$warmup),
                   numericInput(session$ns("sampling"), label = "Sampling iterations", value = rv_deconvo_setting$sampling),
                   numericInput(session$ns("chains"), label = "Number of chains", value = rv_deconvo_setting$chains)
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

        warmup <- input$warmup
        sampling <- input$sampling
        chains <- input$chains

        if(all(is.integer(warmup), is.integer(sampling), is.integer(chains))){
          rv_deconvo_setting$warmup <- warmup
          rv_deconvo_setting$sampling <- sampling
          rv_deconvo_setting$chains <- chains
          removeModal()
        }
      })

      output$anameth_deconvo <- renderUI({
        anameth_all <- list.files(paste0(path_pack, "/extdata/anameth"))
        anameth_all <- gsub(".{4}$", "", anameth_all)
        selectInput(session$ns("anameth_deconvo"), label = h4("Analysis method"), choices = anameth_all, selected = anameth_all[1], selectize = FALSE)
      })

      observeEvent(input$anameth_deconvo, {
        # Load the selected analysis method
        rv_deconvo_setting$dt_anameth <- dt_anameth <- fread(paste0(path_pack, "/extdata/anameth/", input$anameth_deconvo, ".csv"))

        # Create a data.table for information on markers
        dt_at <- data.table(Marker = dt_anameth$Subcategory[dt_anameth$Category == "AT"],
                            AT = as.numeric(dt_anameth$Value[dt_anameth$Category == "AT"]))
        dt_info_mk <- copy(dt_at)
        dt_info_mk <- cbind(
          Ignore = rep(FALSE, nrow(dt_info_mk)),
          dt_info_mk
        )

        output$dt_info_mk <- renderRHandsontable({
          rhandsontable(dt_info_mk, rowHeaders = FALSE) %>%
            hot_table(fillHandle = FALSE, rowHeights = 30) %>%
            hot_col("Ignore", type = "checkbox") %>%
            hot_col("Marker", type = "text", readOnly = TRUE) %>%
            hot_col("AT", type = "numeric", format = "0", readOnly = TRUE)
        })
      })

      observeEvent(input$act_analyze_deconvo, {
        waiter_show(html = spin_3k(), color = "white")

        # Fix setting
        hnc_deconvo <- input$noc_deconvo
        dt_csp <- rv_deconvo_setting$dt_csp
        dt_af <- rv_deconvo_setting$dt_af
        dt_anameth <- rv_deconvo_setting$dt_anameth
        warmup <- rv_deconvo_setting$warmup
        sampling <- rv_deconvo_setting$sampling
        chains <- rv_deconvo_setting$chains
        anameth <- input$anameth_deconvo

        # Extract information on the kit
        dt_kit <- fread(paste0(path_pack, "/extdata/kit/", dt_anameth$Value[dt_anameth$Category == "Kit"], ".csv"))
        dt_kit <- change_data_type_kit(dt_kit)
        dt_kit_auto <- extract_kit_auto(dt_kit)
        markers_auto <- unique(dt_kit_auto[, Marker])

        load(paste0(path_pack, "/extdata/param/", dt_anameth$Value[dt_anameth$Category == "Parameter"], ".Rdata"))
        dt_al_prob <- data.table(Subcategory = dt_anameth$Subcategory[dt_anameth$Category == "Allele_probability"],
                                 Value = dt_anameth$Value[dt_anameth$Category == "Allele_probability"])
        dt_gt_fltr <- data.table(Subcategory = dt_anameth$Subcategory[dt_anameth$Category == "Genotype_filter"],
                                 Value = dt_anameth$Value[dt_anameth$Category == "Genotype_filter"])
        dt_dropin <- data.table(Subcategory = dt_anameth$Subcategory[dt_anameth$Category == "Drop-in"],
                                Value = dt_anameth$Value[dt_anameth$Category == "Drop-in"])

        dt_info_mk <- hot_to_r(input$dt_info_mk)

        # Check input files
        error_message <- check_main_file(dt_csp, NULL, dt_af, markers_auto, method = "deconvo")
        if(error_message != ""){
          waiter_hide()
          showModal(modalDialog(title = "Error", HTML(error_message), easyClose = TRUE, footer = NULL))
        }else{
          # Check registration of the stutter ratio
          message_sr <- check_csp_sr(dt_csp, dt_info_mk, dt_sr)

          if(message_sr == "error"){
            waiter_hide()
            showModal(modalDialog(title = "Error", "All markers are ignored!", easyClose = TRUE, footer = NULL))
          }else{
            data_deconvo <- perform_deconvo(dt_csp, dt_af, dt_kit, dt_info_mk, dt_par, dt_sr, dt_al_prob, dt_gt_fltr, dt_dropin,
                                            hnc_deconvo, warmup, sampling, chains)
            if(!is.null(data_deconvo$error_message)){
              waiter_hide()
              showModal(modalDialog(title = "Error", HTML(data_deconvo$error_message), easyClose = TRUE, footer = NULL))
            }else{
              data_deconvo$message_sr <- message_sr
              data_deconvo$fn_csp <- rv_deconvo_setting$fn_csp
              data_deconvo$fn_af <- rv_deconvo_setting$fn_af
              data_deconvo$anameth <- anameth
              rv_deconvo$data_deconvo <- data_deconvo

              waiter_hide()
            }
          }
        }
      })

      return(rv_deconvo)
    }
  )
}
