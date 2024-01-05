#' deconvo_setting_ui
#'
#' @description The function to create the ui module for setting conditions of deconvolution
deconvo_setting_ui <- function(id){
  ns <- NS(id)

  tabPanel("Setting",
           br(),
           sidebarLayout(
             sidebarPanel(width = 3,
                          numericInput(ns("noc_deconvo"), label = h4("Number of contributors"), value = NULL),
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
                                h4("Conditions to remove unrealistic genotype combinations")
                         )
                       )
             )
           )
  )
}

#' deconvo_setting_server
#'
#' @description The function to create the server module for setting conditions of deconvolution
deconvo_setting_server <- function(id, path_pack){
  moduleServer(
    id,
    function(input, output, session){
      rv_deconvo_setting <- reactiveValues()
      rv_deconvo_setting$dt_analysis_method <- NULL
      rv_deconvo_setting$dt_parameter <- NULL
      rv_deconvo_setting$dt_alcor <- NULL
      rv_deconvo_setting$dt_ae <- NULL
      rv_deconvo_setting$dt_bsr <- NULL
      rv_deconvo_setting$dt_fsr <- NULL
      rv_deconvo_setting$dt_dsr <- NULL
      rv_deconvo_setting$dt_m2sr <- NULL
      rv_deconvo_setting$dt_gamma <- NULL

      output$analysis_method_deconvo <- renderUI({
        analysis_methods <- list.files(paste0(path_pack, "/extdata/analysis_method"))
        analysis_methods <- gsub(".{4}$", "", analysis_methods)
        selectInput(session$ns("analysis_method_deconvo"), label = h4("Analysis method"), choices = analysis_methods, selected = analysis_methods[1])
      })

      observeEvent(input$analysis_method_deconvo, {
        # Load the selected analysis method
        isolate(analysis_method <- input$analysis_method_deconvo)
        rv_deconvo_setting$dt_analysis_method <- dt_analysis_method <- fread(paste0(path_pack, "/extdata/analysis_method/", analysis_method, ".csv"))

        # Extract data per column
        names <- dt_analysis_method[, Name]
        values <- dt_analysis_method[, Value]

        # Load the parameters and the corrected alleles
        rv_deconvo_setting$dt_parameter <- fread(paste0(path_pack, "/extdata/parameter/", values[names == "Parameter"], ".csv"))
        rv_deconvo_setting$dt_alcor <- fread(paste0(path_pack, "/extdata/corrected_allele/", values[names == "Corrected allele"], ".csv"))

        # Extract analytical thresholds
        index_at <- grep("AT_", names)
        dt_at <- data.table(Marker = gsub("AT_", "", names[index_at]),
                            AT = values[index_at])

        output$dt_at <- renderDataTable(server = FALSE, {
          datatable(
            dt_at,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_at), ordering = FALSE, autowidth = TRUE),
            rownames = FALSE
          )
        })

        output$kit <- renderText({values[names == "Kit"]})
        output$alprob <- renderText({values[names == "Allele probability"]})
        output$theta <- renderText({values[names == "Theta"]})
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
        isolate(dt_ae <- rv_deconvo_setting$dt_ae)
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
        isolate(dt_bsr <- rv_deconvo_setting$dt_bsr)
        isolate(mk_bsr <- input$mk_bsr)
        setkey(dt_bsr, Marker)
        dt_bsr_select <- dt_bsr[mk_bsr]
        model_bsr_select <- dt_bsr_select[, Model]

        if(model_bsr_select != "Not consider"){
          # Extract corrected alleles
          isolate(dt_alcor <- rv_deconvo_setting$dt_alcor)
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
        isolate(dt_fsr <- rv_deconvo_setting$dt_fsr)
        isolate(mk_fsr <- input$mk_fsr)
        setkey(dt_fsr, Marker)
        dt_fsr_select <- dt_fsr[mk_fsr]
        model_fsr_select <- dt_fsr_select[, Model]

        if(model_fsr_select != "Not consider"){
          # Extract corrected alleles
          isolate(dt_alcor <- rv_deconvo_setting$dt_alcor)
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
        isolate(dt_dsr <- rv_deconvo_setting$dt_dsr)
        isolate(mk_dsr <- input$mk_dsr)
        setkey(dt_dsr, Marker)
        dt_dsr_select <- dt_dsr[mk_dsr]
        model_dsr_select <- dt_dsr_select[, Model]

        if(model_dsr_select != "Not consider"){
          # Extract corrected alleles
          isolate(dt_alcor <- rv_deconvo_setting$dt_alcor)
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
        isolate(dt_m2sr <- rv_deconvo_setting$dt_m2sr)
        isolate(mk_m2sr <- input$mk_m2sr)
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
        isolate(dt_gamma <- rv_deconvo_setting$dt_gamma)
        datatable(
          dt_gamma,
          selection = "none",
          options = list(dom = "t", iDisplayLength = nrow(dt_gamma), ordering = FALSE, autoWidth = TRUE,
                         columnDefs = list(list(targets = c(0:7, 9:14), visible = FALSE))),
          rownames = FALSE
        ) %>%
          formatRound(columns = c("Beta"), digits = 3)
      })
    }
  )
}
