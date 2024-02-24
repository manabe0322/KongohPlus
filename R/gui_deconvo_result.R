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
                      h4("Estimated parameter"),
                      br(),
                      dataTableOutput(ns("dt_mpdeg")),
                      br(),
                      textOutput(ns("tph")),
                      br(),
                      textOutput(ns("lambda")),
                      br(),
                      actionButton(ns("act_ae"), label = "AE")
                    ),
                    wellPanel(
                      h4("Diagnostics"),
                      br(),
                      actionButton(ns("act_rhat"), label = "Rhat"),
                      br(),
                      actionButton(ns("act_traceplot_lp"), label = "Trace plot")
                    ),
                    br(),
                    actionButton(ns("act_proceed_lr"), label = "Likelihood ratio", class = "btn btn-primary btn-lg")
             ),
             column(6,
                    wellPanel(
                      h4("Probabilistic genotyping"),
                      br(),
                      uiOutput(ns("mk_pg")),
                      br(),
                      dataTableOutput(ns("dt_pg")),
                      br()
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
      results <- rv_deconvo$results
      if(!is.null(results)){
        dt_pg_list <- results$dt_pg_list
        dt_log10like <- results$dt_log10like
        dt_mpdeg <- results$dt_mpdeg
        tph <- results$tph
        lambda <- results$lambda
        dt_ae <- results$dt_ae
        dt_rhat <- results$dt_rhat
        traceplot <- results$traceplot

        hnc <- nrow(dt_mpdeg)
        markers_auto <- names(dt_pg_list)

        output$dt_mpdeg <- renderDataTable(server = FALSE, {
          datatable(
            dt_mpdeg,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_mpdeg), ordering = FALSE, autowidth = TRUE),
            rownames = FALSE
          ) %>%
            formatSignif(columns = c("MP", "Deg"), digits = 3)
        })

        output$tph <- renderText({paste0("tph: ", round(tph, 0))})

        output$lambda <- renderText({paste0("lambda: ", paste(round(lambda, 3), collapse = ", "))})

        observeEvent(input$act_ae, {
          showModal(modalDialog(
            title = "AE",
            dataTableOutput(session$ns("dt_ae")),
            easyClose = TRUE,
            footer = NULL
          ))
        })

        output$dt_ae <- renderDataTable(server = FALSE, {
          datatable(
            dt_ae,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_ae), ordering = FALSE, autoWidth = TRUE),
            rownames = FALSE
          ) %>%
            formatRound(columns = c("AE"), digits = 3)
        })

        observeEvent(input$act_rhat, {
          showModal(modalDialog(
            title = "Rhat",
            dataTableOutput(session$ns("dt_rhat")),
            easyClose = TRUE,
            footer = NULL
          ))
        })

        output$dt_rhat <- renderDataTable(server = FALSE, {
          datatable(
            dt_rhat,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_rhat), ordering = FALSE, autoWidth = TRUE),
            rownames = FALSE
          ) %>%
            formatRound(columns = c("Rhat"), digits = 3)
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

        output$mk_pg <- renderUI({
          selectInput(session$ns("mk_pg"), label = h4("Marker"), choices = markers_auto, selected = markers_auto[1])
        })

        observeEvent(input$mk_pg, {
          mk_select <- input$mk_pg

          output$dt_pg <- renderDataTable(server = FALSE, {
            datatable(
              dt_pg_list[[mk_select]],
              selection = list(mode = "single", target = "row"),
              options = list(iDisplayLength = 25, ordering = FALSE, autowidth = TRUE,
                             columnDefs = list(list(targets = c(hnc + 1, hnc + 2), visible = FALSE))),
              rownames = FALSE
            ) %>%
              formatSignif(columns = c("Weight"), digits = 3)
          })
        })

        output$dt_log10like <- renderDataTable(server = FALSE, {
          datatable(
            dt_log10like,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_log10like), ordering = FALSE, autowidth = TRUE),
            rownames = FALSE
          ) %>%
            formatSignif(columns = c("Likelihood"), digits = 3)
        })
      }
    }
  )
}
