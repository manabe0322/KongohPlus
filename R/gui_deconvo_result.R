#' deconvo_result_ui
#'
#' @description The function to create the ui module for results of deconvolution
deconvo_result_ui <- function(id){
  ns <- NS(id)

  tabPanel("Result",
           br(),
           fluidRow(
             column(3,
                    wellPanel(
                      h4("Crime stain profile"),
                      textOutput(ns("fn_csp")),
                      br(),
                      h4("Allele frequency"),
                      textOutput(ns("fn_af"))
                    ),
                    wellPanel(
                      h4("Warnings"),
                      htmlOutput(ns("warning_message"))
                    )
             ),
             column(6,
                    wellPanel(
                      h4("Parameters and diagnostics"),
                      br(),
                      fluidRow(
                        column(8,
                               dataTableOutput(ns("dt_mpdeg"))
                        ),
                        column(4,
                               br(),
                               actionButton(ns("act_pardiag"), label = "All parameters"),
                               br(),
                               br(),
                               actionButton(ns("act_traceplot_lp"), label = "Trace plot")
                        )
                      )
                    )
             ),
             column(3,
                    wellPanel(
                      h4("Likelihood"),
                      br(),
                      dataTableOutput(ns("dt_log10like")),
                      br()
                    )
             )
           )
  )
}

#' deconvo_result_server
#'
#' @description The function to create the server module for results of deconvolution
deconvo_result_server <- function(id, rv_deconvo){
  moduleServer(
    id,
    function(input, output, session){
      data_deconvo <- rv_deconvo$data_deconvo
      if(!is.null(data_deconvo)){
        dt_csp <- data_deconvo$dt_csp
        dt_af <- data_deconvo$dt_af
        fn_csp <- data_deconvo$fn_csp
        fn_af <- data_deconvo$fn_af
        hnc <- data_deconvo$hnc
        dt_log10like <- data_deconvo$dt_log10like
        dt_mpdeg <- data_deconvo$dt_mpdeg
        dt_diag <- data_deconvo$dt_diag
        traceplot <- data_deconvo$traceplot
        message_sr <- data_deconvo$message_sr

        output$fn_csp <- renderText({paste0(fn_csp)})
        output$fn_af <- renderText({paste0(fn_af)})

        output$warning_message <- renderUI({HTML(message_sr)})

        output$dt_mpdeg <- renderDataTable(server = FALSE, {
          datatable(
            dt_mpdeg,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_mpdeg), ordering = FALSE, autowidth = TRUE),
            rownames = FALSE,
            colnames = c("Contributor", "Mixture proportion", "Degradation")
          ) %>%
            formatSignif(columns = c("MP", "Deg"), digits = 3)
        })

        observeEvent(input$act_pardiag, {
          showModal(modalDialog(
            title = "Parameters and diagnostics",
            dataTableOutput(session$ns("dt_diag")),
            easyClose = TRUE,
            footer = NULL,
            size = "l"
          ))
        })

        output$dt_diag <- renderDataTable(server = FALSE, {
          datatable(
            dt_diag,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_diag), ordering = FALSE, autowidth = TRUE,
                           columnDefs = list(list(width = "200px", targets = c(0)),
                                             list(width = "100px", targets = c(1, 2, 3, 4)))),
            rownames = FALSE,
            colnames = c("Parameter", "Contributor", "Marker", "Estimated value", "Rhat")
          ) %>%
            formatSignif(columns = c("Estimated_value", "Rhat"), digits = 3)
        })

        observeEvent(input$act_traceplot_lp, {
          showModal(modalDialog(
            title = "Trace plot",
            plotOutput(session$ns("traceplot_lp")),
            footer = tagList(
              modalButton("Close")
            ),
            size = "l"
          ))
        })

        output$traceplot_lp <- renderPlot({
          traceplot
        })

        output$dt_log10like <- renderDataTable(server = FALSE, {
          datatable(
            dt_log10like,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_log10like), ordering = FALSE, autowidth = TRUE),
            rownames = FALSE,
            colnames = c("Marker", "Log<sub>10</sub>(Likelihood)"),
            escape = FALSE
          ) %>%
            formatSignif(columns = c("Log10_Likelihood"), digits = 3)
        })
      }
    }
  )
}
