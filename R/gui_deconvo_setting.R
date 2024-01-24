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
                          uiOutput(ns("analysis_method_deconvo")),
                          br(),
                          actionButton(ns("act_analyze_deconvo"), label = "Analyze", class = "btn btn-primary btn-lg")
             ),
             mainPanel(width = 9,
                       fluidRow(
                         column(4,
                                dataTableOutput(ns("dt_at")),
                                br()
                         ),
                         column(8,
                                h4("Kit"),
                                textOutput(ns("kit")),
                                br(),
                                h4("Allele probability"),
                                textOutput(ns("alprob")),
                                br(),
                                h4("Theta"),
                                textOutput(ns("theta")),
                                br(),
                                h4("Parameters"),
                                actionButton(ns("act_ae"), label = "AE"),
                                actionButton(ns("act_bsr"), label = "BSR"),
                                actionButton(ns("act_fsr"), label = "FSR"),
                                actionButton(ns("act_dsr"), label = "DSR"),
                                actionButton(ns("act_m2sr"), label = "M2SR"),
                                actionButton(ns("act_gamma"), label = "Gamma"),
                                br(),
                                br(),
                                h4("Conditions to remove unrealistic genotype combinations"),
                                textOutput(ns("bsr_fltr")),
                                textOutput(ns("fsr_fltr")),
                                textOutput(ns("dsr_fltr")),
                                textOutput(ns("m2sr_fltr")),
                                textOutput(ns("phr_fltr")),
                                textOutput(ns("st"))
                         )
                       )
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
      rv_deconvo_setting$dt_analysis_method <- NULL
      rv_deconvo_setting$dt_parameter <- NULL
      rv_deconvo_setting$dt_alcor <- NULL
      rv_deconvo_setting$dt_kit <- NULL
      rv_deconvo_setting$dt_at <- NULL
      rv_deconvo_setting$dt_ae <- NULL
      rv_deconvo_setting$dt_bsr <- NULL
      rv_deconvo_setting$dt_fsr <- NULL
      rv_deconvo_setting$dt_dsr <- NULL
      rv_deconvo_setting$dt_m2sr <- NULL
      rv_deconvo_setting$dt_gamma <- NULL
      rv_deconvo_setting$bsr_fltr <- NULL
      rv_deconvo_setting$fsr_fltr <- NULL
      rv_deconvo_setting$dsr_fltr <- NULL
      rv_deconvo_setting$m2sr_fltr <- NULL
      rv_deconvo_setting$phr_fltr <- NULL
      rv_deconvo_setting$st <- NULL

      rv_deconvo_setting$dt_csp <- NULL
      rv_deconvo_setting$dt_ref <- NULL
      rv_deconvo_setting$dt_af <- NULL

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

      output$analysis_method_deconvo <- renderUI({
        analysis_methods <- list.files(paste0(path_pack, "/extdata/analysis_method"))
        analysis_methods <- gsub(".{4}$", "", analysis_methods)
        selectInput(session$ns("analysis_method_deconvo"), label = h4("Analysis method"), choices = analysis_methods, selected = analysis_methods[1])
      })

      observeEvent(input$analysis_method_deconvo, {
        # Load the selected analysis method
        analysis_method <- input$analysis_method_deconvo
        rv_deconvo_setting$dt_analysis_method <- dt_analysis_method <- fread(paste0(path_pack, "/extdata/analysis_method/", analysis_method, ".csv"))

        # Extract data per column
        names <- dt_analysis_method[, Name]
        values <- dt_analysis_method[, Value]

        # Load the parameters, the corrected alleles, and the kit information
        rv_deconvo_setting$dt_parameter <- fread(paste0(path_pack, "/extdata/parameter/", values[names == "Parameter"], ".csv"))
        rv_deconvo_setting$dt_alcor <- fread(paste0(path_pack, "/extdata/corrected_allele/", values[names == "Corrected allele"], ".csv"))
        rv_deconvo_setting$dt_kit <- fread(paste0(path_pack, "/extdata/kit/", values[names == "Kit"], ".csv"))

        # Extract analytical thresholds
        index_at <- grep("AT_", names)
        rv_deconvo_setting$dt_at <- dt_at <- data.table(Marker = gsub("AT_", "", names[index_at]),
                                                        AT = as.numeric(values[index_at]))

        output$dt_at <- renderDataTable(server = FALSE, {
          datatable(
            dt_at,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_at), ordering = FALSE, autowidth = TRUE),
            rownames = FALSE
          )
        })

        rv_deconvo_setting$bsr_fltr <- bsr_fltr <- as.numeric(values[names == "Filter of back stutter ratio"])
        rv_deconvo_setting$fsr_fltr <- fsr_fltr <- as.numeric(values[names == "Filter of forward stutter ratio"])
        rv_deconvo_setting$dsr_fltr <- dsr_fltr <- as.numeric(values[names == "Filter of double-back stutter ratio"])
        rv_deconvo_setting$m2sr_fltr <- m2sr_fltr <- as.numeric(values[names == "Filter of minus 2-nt stutter ratio"])
        rv_deconvo_setting$phr_fltr <- phr_fltr <- as.numeric(values[names == "Filter of peak height ratio"])
        rv_deconvo_setting$st <- st <- as.numeric(values[names == "Stochastic threshold"])

        output$kit <- renderText({values[names == "Kit"]})
        output$alprob <- renderText({values[names == "Allele probability"]})
        output$theta <- renderText({values[names == "Theta"]})
        output$bsr_fltr <- renderText({paste0("Filter of back stutter ratio: ", bsr_fltr)})
        output$fsr_fltr <- renderText({paste0("Filter of forward stutter ratio: ", fsr_fltr)})
        output$dsr_fltr <- renderText({paste0("Filter of double-back stutter ratio: ", dsr_fltr)})
        output$m2sr_fltr <- renderText({paste0("Filter of minus 2-nt stutter ratio: ", m2sr_fltr)})
        output$phr_fltr <- renderText({paste0("Filter of peak height ratio: ", phr_fltr)})
        output$st <- renderText({paste0("Stochastic threshold: ", st)})
      })

      observeEvent(input$act_ae, {
        dt_parameter <- rv_deconvo_setting$dt_parameter
        setkey(dt_parameter, Parameter)
        rv_deconvo_setting$dt_ae <- dt_parameter["AE"]

        showModal(modalDialog(
          title = "Locus-specific amplification efficiency",
          dataTableOutput(session$ns("dt_ae")),
          easyClose = TRUE,
          footer = NULL
        ))
      })

      output$dt_ae <- renderDataTable(server = FALSE, {
        dt_ae <- rv_deconvo_setting$dt_ae
        datatable(
          dt_ae,
          selection = "none",
          options = list(dom = "t", iDisplayLength = nrow(dt_ae), ordering = FALSE, autoWidth = TRUE,
                         columnDefs = list(list(targets = c(1, 2, 4, 5, 7:14), visible = FALSE))),
          rownames = FALSE
        ) %>%
          formatRound(columns = c("Myu"), digits = 3) %>%
          formatRound(columns = c("Sigma"), digits = 0)
      })

      observeEvent(input$act_bsr, {
        print(input$act_bsr)

        # Extract parameters for BSR
        dt_parameter <- rv_deconvo_setting$dt_parameter
        setkey(dt_parameter, Parameter)
        rv_deconvo_setting$dt_bsr <- dt_bsr <- dt_parameter["BSR"]
        markers <- dt_bsr[, Marker]

        showModal(modalDialog(
          title = "Back stutter ratio",
          selectInput(session$ns("mk_bsr"), label = h4("Marker"), choices = markers, selected = markers[1]),
          h4("Parameter"),
          dataTableOutput(session$ns("dt_param_bsr")),
          br(),
          h4("Corrected alleles and average back stutter ratios"),
          dataTableOutput(session$ns("dt_alcor_bsr")),
          easyClose = TRUE,
          footer = NULL
        ))
      })

      observeEvent(input$mk_bsr, {
        # Extract parameters
        dt_bsr <- rv_deconvo_setting$dt_bsr
        mk_bsr <- input$mk_bsr
        setkey(dt_bsr, Marker)
        dt_bsr_select <- dt_bsr[mk_bsr]
        model_bsr_select <- dt_bsr_select[, Model]

        if(model_bsr_select != "Not consider"){
          # Extract corrected alleles
          dt_alcor <- rv_deconvo_setting$dt_alcor
          setkey(dt_alcor, Marker)
          dt_alcor_select <- dt_alcor[mk_bsr]
        }else{
          dt_alcor_select <- NULL
        }

        output$dt_param_bsr <- renderDataTable(server = FALSE, {
          datatable(
            dt_bsr_select,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_bsr_select), ordering = FALSE, autoWidth = TRUE,
                           columnDefs = list(list(targets = c(0, 1, 8:14), visible = FALSE))),
            rownames = FALSE
          ) %>%
            formatRound(columns = c("Myu", "Slope", "Intercept"), digits = 3) %>%
            formatRound(columns = c("Sigma", "MinUS"), digits = 0)
        })

        output$dt_alcor_bsr <- renderDataTable(server = FALSE, {
          col_names <- names(dt_alcor_select)
          col_names[4] <- "Corrected allele"
          col_names[5] <- "Average back stutter ratio"
          datatable(
            dt_alcor_select,
            colnames = col_names,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_alcor_select), ordering = FALSE, autoWidth = TRUE,
                           columnDefs = list(list(targets = c(0, 1, 5:10), visible = FALSE))),
            rownames = FALSE
          ) %>%
            formatRound(columns = c("Corrected_BSR", "Average_BSR"), digits = 3)
        })
      })

      observeEvent(input$act_fsr, {
        # Extract parameters for FSR
        dt_parameter <- rv_deconvo_setting$dt_parameter
        setkey(dt_parameter, Parameter)
        rv_deconvo_setting$dt_fsr <- dt_fsr <- dt_parameter["FSR"]
        markers <- dt_fsr[, Marker]

        showModal(modalDialog(
          title = "Forward stutter ratio",
          selectInput(session$ns("mk_fsr"), label = h4("Marker"), choices = markers, selected = markers[1]),
          h4("Parameter"),
          dataTableOutput(session$ns("dt_param_fsr")),
          br(),
          h4("Corrected alleles and average forward stutter ratios"),
          dataTableOutput(session$ns("dt_alcor_fsr")),
          easyClose = TRUE,
          footer = NULL
        ))
      })

      observeEvent(input$mk_fsr, {
        # Extract parameters
        dt_fsr <- rv_deconvo_setting$dt_fsr
        mk_fsr <- input$mk_fsr
        setkey(dt_fsr, Marker)
        dt_fsr_select <- dt_fsr[mk_fsr]
        model_fsr_select <- dt_fsr_select[, Model]

        if(model_fsr_select != "Not consider"){
          # Extract corrected alleles
          dt_alcor <- rv_deconvo_setting$dt_alcor
          setkey(dt_alcor, Marker)
          dt_alcor_select <- dt_alcor[mk_fsr]
        }else{
          dt_alcor_select <- NULL
        }

        output$dt_param_fsr <- renderDataTable(server = FALSE, {
          datatable(
            dt_fsr_select,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_fsr_select), ordering = FALSE, autoWidth = TRUE,
                           columnDefs = list(list(targets = c(0, 1, 8:14), visible = FALSE))),
            rownames = FALSE
          ) %>%
            formatRound(columns = c("Myu", "Slope", "Intercept"), digits = 3) %>%
            formatRound(columns = c("Sigma", "MinUS"), digits = 0)
        })

        output$dt_alcor_fsr <- renderDataTable(server = FALSE, {
          col_names <- names(dt_alcor_select)
          col_names[7] <- "Corrected allele"
          col_names[8] <- "Average forward stutter ratio"
          datatable(
            dt_alcor_select,
            colnames = col_names,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_alcor_select), ordering = FALSE, autoWidth = TRUE,
                           columnDefs = list(list(targets = c(0, 1, 3:5, 8:10), visible = FALSE))),
            rownames = FALSE
          ) %>%
            formatRound(columns = c("Corrected_FSR", "Average_FSR"), digits = 3)
        })
      })

      observeEvent(input$act_dsr, {
        # Extract parameters for DSR
        dt_parameter <- rv_deconvo_setting$dt_parameter
        setkey(dt_parameter, Parameter)
        rv_deconvo_setting$dt_dsr <- dt_dsr <- dt_parameter["DSR"]
        markers <- dt_dsr[, Marker]

        showModal(modalDialog(
          title = "Double-back stutter ratio",
          selectInput(session$ns("mk_dsr"), label = h4("Marker"), choices = markers, selected = markers[1]),
          h4("Parameter"),
          dataTableOutput(session$ns("dt_param_dsr")),
          br(),
          h4("Corrected alleles and average double-back stutter ratios"),
          dataTableOutput(session$ns("dt_alcor_dsr")),
          easyClose = TRUE,
          footer = NULL
        ))
      })

      observeEvent(input$mk_dsr, {
        # Extract parameters
        dt_dsr <- rv_deconvo_setting$dt_dsr
        mk_dsr <- input$mk_dsr
        setkey(dt_dsr, Marker)
        dt_dsr_select <- dt_dsr[mk_dsr]
        model_dsr_select <- dt_dsr_select[, Model]

        if(model_dsr_select != "Not consider"){
          # Extract corrected alleles
          dt_alcor <- rv_deconvo_setting$dt_alcor
          setkey(dt_alcor, Marker)
          dt_alcor_select <- dt_alcor[mk_dsr]
        }else{
          dt_alcor_select <- NULL
        }

        output$dt_param_dsr <- renderDataTable(server = FALSE, {
          datatable(
            dt_dsr_select,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_dsr_select), ordering = FALSE, autoWidth = TRUE,
                           columnDefs = list(list(targets = c(0, 1, 8:14), visible = FALSE))),
            rownames = FALSE
          ) %>%
            formatRound(columns = c("Myu", "Slope", "Intercept"), digits = 3) %>%
            formatRound(columns = c("Sigma", "MinUS"), digits = 0)
        })

        output$dt_alcor_dsr <- renderDataTable(server = FALSE, {
          col_names <- names(dt_alcor_select)
          col_names[10] <- "Corrected allele"
          col_names[11] <- "Average double-back stutter ratio"
          datatable(
            dt_alcor_select,
            colnames = col_names,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_alcor_select), ordering = FALSE, autoWidth = TRUE,
                           columnDefs = list(list(targets = c(0, 1, 3:8), visible = FALSE))),
            rownames = FALSE
          ) %>%
            formatRound(columns = c("Corrected_DSR", "Average_DSR"), digits = 3)
        })
      })

      observeEvent(input$act_m2sr, {
        # Extract parameters for M2SR
        dt_parameter <- rv_deconvo_setting$dt_parameter
        setkey(dt_parameter, Parameter)
        rv_deconvo_setting$dt_m2sr <- dt_m2sr <- dt_parameter["M2SR"]
        markers <- dt_m2sr[, Marker]

        showModal(modalDialog(
          title = "Minus 2-nt stutter ratio",
          selectInput(session$ns("mk_m2sr"), label = h4("Marker"), choices = markers, selected = markers[1]),
          dataTableOutput(session$ns("dt_param_m2sr")),
          easyClose = TRUE,
          footer = NULL
        ))
      })

      observeEvent(input$mk_m2sr, {
        # Extract parameters
        dt_m2sr <- rv_deconvo_setting$dt_m2sr
        mk_m2sr <- input$mk_m2sr
        setkey(dt_m2sr, Marker)
        dt_m2sr_select <- dt_m2sr[mk_m2sr]
        model_m2sr_select <- dt_m2sr_select[, Model]

        output$dt_param_m2sr <- renderDataTable(server = FALSE, {
          datatable(
            dt_m2sr_select,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_m2sr_select), ordering = FALSE, autoWidth = TRUE,
                           columnDefs = list(list(targets = c(0, 1, 4, 5, 7:14), visible = FALSE))),
            rownames = FALSE
          ) %>%
            formatRound(columns = c("Myu"), digits = 3) %>%
            formatRound(columns = c("Sigma"), digits = 0)
        })
      })

      observeEvent(input$act_gamma, {
        dt_parameter <- rv_deconvo_setting$dt_parameter
        setkey(dt_parameter, Parameter)
        rv_deconvo_setting$dt_gamma <- dt_parameter["Gamma"]

        showModal(modalDialog(
          title = "Gamma distribution",
          dataTableOutput(session$ns("dt_gamma")),
          easyClose = TRUE,
          footer = NULL
        ))
      })

      output$dt_gamma <- renderDataTable(server = FALSE, {
        dt_gamma <- rv_deconvo_setting$dt_gamma
        datatable(
          dt_gamma,
          selection = "none",
          options = list(dom = "t", iDisplayLength = nrow(dt_gamma), ordering = FALSE, autoWidth = TRUE,
                         columnDefs = list(list(targets = c(0:7, 9:14), visible = FALSE))),
          rownames = FALSE
        ) %>%
          formatRound(columns = c("Beta"), digits = 3)
      })

      observeEvent(input$act_analyze_deconvo, {
        # Fix setting
        noc_deconvo <- input$noc_deconvo

        # Fix data.table
        dt_csp <- rv_deconvo_setting$dt_csp
        dt_ref <- rv_deconvo_setting$dt_ref
        dt_af <- rv_deconvo_setting$dt_af
        #dt_analysis_method <- rv_deconvo_setting$dt_analysis_method
        dt_parameter <- rv_deconvo_setting$dt_parameter
        dt_alcor <- rv_deconvo_setting$dt_alcor
        dt_kit <- rv_deconvo_setting$dt_kit
        dt_at <- rv_deconvo_setting$dt_at
        bsr_fltr <- rv_deconvo_setting$bsr_fltr
        fsr_fltr <- rv_deconvo_setting$fsr_fltr
        dsr_fltr <- rv_deconvo_setting$dsr_fltr
        m2sr_fltr <- rv_deconvo_setting$m2sr_fltr
        phr_fltr <- rv_deconvo_setting$phr_fltr
        st <- rv_deconvo_setting$st

        # Check input files and setting
        error_message <- check_main_file(dt_csp, dt_ref, dt_af, dt_kit, dt_at, dt_alcor, method = "deconvo")
        if(error_message != ""){
          showModal(modalDialog(title = "Error", HTML(error_message), easyClose = TRUE, footer = NULL))
        }else{
          perform_deconvo(path_pack, dt_csp, dt_af, dt_parameter, dt_alcor, dt_kit, dt_at, bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr, phr_fltr, st, noc_deconvo)
        }
      })
    }
  )
}
