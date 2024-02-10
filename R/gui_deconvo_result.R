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
                      br()
                    ),
                    wellPanel(
                      h4("Diagnostics"),
                      br(),
                      textOutput(ns("rhat_lp")),
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
                      dataTableOutput(ns("dt_likelihood")),
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
        dt_likelihood <- results$dt_likelihood
        dt_mpdeg <- results$dt_mpdeg
        tph <- results$tph
        rhat_lp <- results$rhat_lp
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

        output$rhat_lp <- renderText({paste0("Rhat: ", round(rhat_lp, 3))})

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

        output$dt_likelihood <- renderDataTable(server = FALSE, {
          datatable(
            dt_likelihood,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_likelihood), ordering = FALSE, autowidth = TRUE),
            rownames = FALSE
          ) %>%
            formatSignif(columns = c("Likelihood"), digits = 3)
        })
      }
    }
  )
}
