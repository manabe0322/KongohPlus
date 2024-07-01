#' deconvo_setting_ui
#'
#' @description The function to create the ui module for setting conditions of deconvolution
deconvo_setting_ui <- function(id){
  ns <- NS(id)

  tabPanel("Setting",
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
                       tags$head(
                         #tags$style(HTML("table {table-layout: fixed; td > div.form-group > div.checkbox {margin:0px; bottom:10px;} td > div.form-group {position:absolute;}")),
                         tags$style(".datatables .display {margin-left: 0;}")
                       ),
                       dataTableOutput(ns("dt_info_mk_display"))
             )
           )
  )
}

#' deconvo_setting_server
#'
#' @description The function to create the server module for setting conditions of deconvolution
deconvo_setting_server <- function(id, path_pack, rv_file, session_deconvo_top){
  moduleServer(
    id,
    function(input, output, session){
      shinyInput <- function(FUN, len, id, ...) {
        inputs <- character(len)
        for (i in seq_len(len)) {
          inputs[i] <- as.character(FUN(paste0(id, i), label = NULL, ...))
        }
        inputs
      }

      js <- function(dtid, ns){
        c(
          "$('[id^=checkb]').on('click', function(){",
          "  var id = this.getAttribute('id');",
          "  var i = parseInt(/checkb(\\d+)/.exec(id)[1]);",
          "  var value = $(this).prop('checked');",
          "  var info = [{row: i, col: 1, value: value}];",
          sprintf(
            "Shiny.setInputValue('%s', info);",
            ns(sprintf("%s_cell_edit:DT.cellInfo", dtid))
          ),
          "})"
        )
      }

      rv_deconvo_setting <- reactiveValues()
      rv_deconvo_setting$dt_csp <- NULL
      rv_deconvo_setting$dt_ref <- NULL
      rv_deconvo_setting$dt_af <- NULL
      rv_deconvo_setting$anameth <- NULL
      #rv_deconvo_setting$dt_anameth <- NULL
      rv_deconvo_setting$dt_kit <- NULL
      rv_deconvo_setting$dt_at <- NULL

      rv_deconvo_setting$dt_info_mk <- NULL

      rv_deconvo_setting$dt_param <- NULL
      rv_deconvo_setting$dt_alcor <- NULL
      rv_deconvo_setting$gtc_fltr <- NULL
      rv_deconvo_setting$gtc_prob <- NULL
      #rv_deconvo_setting$dt_ae <- NULL
      #rv_deconvo_setting$dt_bsr <- NULL
      #rv_deconvo_setting$dt_fsr <- NULL
      #rv_deconvo_setting$dt_dsr <- NULL
      #rv_deconvo_setting$dt_m2sr <- NULL
      #rv_deconvo_setting$dt_gamma <- NULL
      rv_deconvo_setting$warmup <- 1200
      rv_deconvo_setting$sampling <- 300
      rv_deconvo_setting$chains <- 4

      rv_deconvo <- reactiveValues()
      rv_deconvo$data_deconvo <- NULL
      rv_deconvo$dt_csp <- NULL
      rv_deconvo$dt_ref <- NULL
      rv_deconvo$dt_af <- NULL
      rv_deconvo$anameth <- NULL
      rv_deconvo$dt_kit <- NULL
      rv_deconvo$dt_at <- NULL
      rv_deconvo$dt_param <- NULL
      rv_deconvo$dt_alcor <- NULL
      rv_deconvo$dt_gtc_fltr <- NULL
      rv_deconvo$dt_gtc_prob <- NULL
      rv_deconvo$warmup <- NULL
      rv_deconvo$sampling <- NULL
      rv_deconvo$chains <- NULL

      observe({
        req(rv_file)
        rv_deconvo_setting$dt_csp <- rv_file$dt_csp
        rv_deconvo_setting$dt_ref <- rv_file$dt_ref
        rv_deconvo_setting$dt_af <- rv_file$dt_af
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

        isolate(warmup <- input$warmup)
        isolate(sampling <- input$sampling)
        isolate(chains <- input$chains)

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
        selectInput(session$ns("anameth_deconvo"), label = h4("Analysis method"), choices = anameth_all, selected = anameth_all[1])
      })

      observeEvent(input$anameth_deconvo, {
        # Load the selected analysis method
        anameth <- input$anameth_deconvo
        rv_deconvo_setting$anameth <- anameth
        dt_anameth <- fread(paste0(path_pack, "/extdata/anameth/", anameth, ".csv"))

        # Extract data per column
        folders <- dt_anameth[, Folder]
        names <- dt_anameth[, Name]

        # Load each data.table
        rv_deconvo_setting$dt_kit <- fread(paste0(path_pack, "/extdata/kit/", names[folders == "kit"], ".csv"))
        rv_deconvo_setting$dt_at <- dt_at <- fread(paste0(path_pack, "/extdata/at/", names[folders == "at"], ".csv"))
        rv_deconvo_setting$dt_param <- fread(paste0(path_pack, "/extdata/param/", names[folders == "param"], ".csv"))
        rv_deconvo_setting$dt_alcor <- fread(paste0(path_pack, "/extdata/alcor/", names[folders == "alcor"], ".csv"))
        rv_deconvo_setting$dt_gtc_fltr <- dt_gtc_fltr <- fread(paste0(path_pack, "/extdata/gtc_fltr/", names[folders == "gtc_fltr"], ".csv"))
        rv_deconvo_setting$dt_gtc_prob <- dt_gtc_prob <- fread(paste0(path_pack, "/extdata/gtc_prob/", names[folders == "gtc_prob"], ".csv"))

        dt_info_mk <- copy(dt_at)
        dt_info_mk <- cbind(
          Ignore = rep(FALSE, nrow(dt_info_mk)),
          dt_info_mk
        )
        rv_deconvo_setting$dt_info_mk <- dt_info_mk

        dt_info_mk_display <- copy(dt_at)
        dt_info_mk_display <- cbind(
          Ignore = shinyInput(checkboxInput, nrow(dt_info_mk_display), "checkb"),
          dt_info_mk_display
        )

        output$dt_info_mk_display <- renderDataTable(server = FALSE, {
          datatable(
            dt_info_mk_display,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_info_mk_display), ordering = FALSE, autoWidth = TRUE, scrollX = TRUE, scrollY = "650px",
                           columnDefs = list(list(targets = 0:2, width = "100px"))),
            rownames = FALSE,
            escape = FALSE,
            callback = JS(js("dt_info_mk_display", session$ns))
          )# %>% formatStyle(0, target = "row", lineHeight = "75%")
        })
      })

      observeEvent(input[["dt_info_mk_display_cell_edit"]], {
        info <- input[["dt_info_mk_display_cell_edit"]]
        dt_info_mk <- rv_deconvo_setting$dt_info_mk
        dt_info_mk[[info[1, "col"]]][info[1, "row"]] <- info[1, "value"]
        rv_deconvo_setting$dt_info_mk <- dt_info_mk
      })

      observeEvent(input$act_analyze_deconvo, {
        # Fix setting
        hnc_deconvo <- input$noc_deconvo
        dt_csp <- rv_deconvo_setting$dt_csp
        dt_ref <- rv_deconvo_setting$dt_ref
        dt_af <- rv_deconvo_setting$dt_af
        anameth <- rv_deconvo_setting$anameth
        dt_kit <- rv_deconvo_setting$dt_kit
        dt_at <- rv_deconvo_setting$dt_at
        dt_param <- rv_deconvo_setting$dt_param
        dt_alcor <- rv_deconvo_setting$dt_alcor
        dt_gtc_fltr <- rv_deconvo_setting$dt_gtc_fltr
        dt_gtc_prob <- rv_deconvo_setting$dt_gtc_prob
        warmup <- rv_deconvo_setting$warmup
        sampling <- rv_deconvo_setting$sampling
        chains <- rv_deconvo_setting$chains

        # Check input files and setting
        error_message <- check_main_file(dt_csp, dt_ref, dt_af, dt_kit, dt_at, dt_alcor, method = "deconvo")
        if(error_message != ""){
          showModal(modalDialog(title = "Error", HTML(error_message), easyClose = TRUE, footer = NULL))
        }else{
          data_deconvo <- perform_deconvo(path_pack, dt_csp, dt_af, dt_kit, dt_at, dt_param, dt_alcor, dt_gtc_fltr, dt_gtc_prob, hnc_deconvo,
                                          warmup, sampling, chains)
          rv_deconvo$data_deconvo <- data_deconvo
          rv_deconvo$dt_csp <- dt_csp
          rv_deconvo$dt_ref <- dt_ref
          rv_deconvo$dt_af <- dt_af
          rv_deconvo$anameth <- anameth
          rv_deconvo$dt_kit <- dt_kit
          rv_deconvo$dt_at <- dt_at
          rv_deconvo$dt_param <- dt_param
          rv_deconvo$dt_alcor <- dt_alcor
          rv_deconvo$dt_gtc_fltr <- dt_gtc_fltr
          rv_deconvo$dt_gtc_prob <- dt_gtc_prob
          rv_deconvo$warmup <- warmup
          rv_deconvo$sampling <- sampling
          rv_deconvo$chains <- chains

          updateTabsetPanel(session = session_deconvo_top, inputId = "tabset_deconvo", selected = "Result")
        }
      })

      return(rv_deconvo)
    }
  )
}
