#' deconvo_result_ui
#'
#' @description The function to create the ui module for results of deconvolution
deconvo_result_ui <- function(id){
  ns <- NS(id)

  tabPanel("Result",
           br(),
           fluidRow(
             column(2,
                    wellPanel(
                      h4(textOutput(ns("title_warnings"))),
                      htmlOutput(ns("warning_message"))
                    )
             ),
             column(2,
                    wellPanel(
                      h4(textOutput(ns("title_like"))),
                      br(),
                      dataTableOutput(ns("dt_log10like")),
                      br()
                    )
             ),
             column(4,
                    wellPanel(
                      h4(textOutput(ns("title_pg"))),
                      fluidRow(
                        column(6,
                               uiOutput(ns("select_meth_pg")),
                               uiOutput(ns("select_mk_pg"))
                        ),
                        column(6,
                               uiOutput(ns("select_noc_pg")),
                               uiOutput(ns("select_contrib_pg"))
                        )
                      ),
                      dataTableOutput(ns("dt_pg"))
                    )
             ),
             column(4,
                    wellPanel(
                      h4(textOutput(ns("title_par_diag"))),
                      actionButton(ns("act_traceplot_lp"), label = "Trace plot"),
                      rHandsontableOutput(ns("dt_diag")),
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
        markers_use <- data_deconvo$markers_use
        message_sr <- data_deconvo$message_sr
        hnc <- data_deconvo$hnc
        dt_log10like <- data_deconvo$dt_log10like
        dt_pg_list <- data_deconvo$dt_pg_list
        dt_diag <- data_deconvo$dt_diag
        traceplot <- data_deconvo$traceplot

        rv_tmp <- reactiveValues()
        rv_tmp$update <- FALSE

        output$title_warnings <- renderText({"Warnings"})
        output$title_like <- renderText({"Likelihood"})
        output$title_pg <- renderText({"Probabilistic genotyping"})
        output$title_par_diag <- renderText({"Parameters and diagnostics"})

        output$warning_message <- renderUI({HTML(message_sr)})

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

        ## Probabilistic genotyping ##

        output$select_meth_pg <- renderUI({
          selectInput(session$ns("select_meth_pg"), label = "Method", choices = c("Per combination", "Per contributor"), selected = "Per combination", selectize = FALSE)
        })

        output$select_noc_pg <- renderUI({
          disabled(selectInput(session$ns("select_noc_pg"), label = "Number of contributors", choices = hnc, selected = hnc[1], selectize = FALSE))
        })

        output$select_mk_pg <- renderUI({
          selectInput(session$ns("select_mk_pg"), label = "Marker", choices = markers_use, selected = markers_use[1], selectize = FALSE)
        })

        output$select_contrib_pg <- renderUI({
          cand_contrib <- paste0("U", 1:hnc)
          selectInput(session$ns("select_contrib_pg"), label = "Contributor", choices = cand_contrib, selected = cand_contrib[1], selectize = FALSE)
        })

        observeEvent(input$select_meth_pg, {
          rv_tmp$update <- TRUE
        })

        observeEvent(input$select_mk_pg, {
          rv_tmp$update <- TRUE
        })

        observeEvent(input$select_contrib_pg, {
          rv_tmp$update <- TRUE
        })

        observe({
          if(rv_tmp$update){
            meth_pg <- input$select_meth_pg
            if(meth_pg == "Per combination"){
              disable("select_contrib_pg")

              output$dt_pg <- renderDataTable(server = FALSE, {
                datatable(
                  dt_pg_list[[input$select_mk_pg]],
                  selection = list(mode = "single", target = "row"),
                  options = list(iDisplayLength = 10, ordering = FALSE, autowidth = TRUE, searching = FALSE,
                                 columnDefs = list(list(targets = c(hnc, hnc + 1, hnc + 2), visible = FALSE))),
                  rownames = FALSE
                ) %>%
                  formatSignif(columns = c("Probability"), digits = 3)
              })
            }else{
              enable("select_contrib_pg")

              dt_pg_select <- dt_pg_list[[input$select_mk_pg]]
              aggregated_data <- aggregate_pg_contrib(dt_pg_select, input$select_contrib_pg)

              output$dt_pg <- renderDataTable(server = FALSE, {
                datatable(
                  aggregated_data,
                  selection = list(mode = "single", target = "row"),
                  options = list(iDisplayLength = 10, ordering = FALSE, autowidth = TRUE, searching = FALSE),
                  rownames = FALSE
                ) %>%
                  formatSignif(columns = c("Probability"), digits = 3)
              })
            }

            rv_tmp$update <- FALSE
          }
        })

        ## Parameters and diagnostics

        output$dt_diag <- renderRHandsontable({
          rhandsontable(dt_diag, rowHeaders = FALSE) %>%
            hot_table(fillHandle = FALSE, rowHeights = 30) %>%
            hot_col("Parameter", type = "text", readOnly = TRUE, width = 175) %>%
            hot_col("Contributor", type = "text", readOnly = TRUE, width = 90) %>%
            hot_col("Marker", type = "text", readOnly = TRUE, width = 100) %>%
            hot_col("Estimated value", type = "numeric", readOnly = TRUE, width = 125) %>%
            hot_col("Rhat", type = "numeric", readOnly = TRUE, width = 50)
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
      }
    }
  )
}
