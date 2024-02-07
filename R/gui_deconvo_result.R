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
                      h4("Estimated parameter")
                    )
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
                      h4("Likelihood")
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

        markers_auto <- names(dt_pg_list)

        output$mk_pg <- renderUI({
          selectInput(session$ns("mk_pg"), label = h4("Marker"), choices = markers_auto, selected = markers_auto[1])
        })

        observeEvent(input$mk_pg, {
          mk_select <- input$mk_pg

          output$dt_pg <- renderDataTable(server = FALSE, {
            datatable(
              dt_pg_list[[mk_select]],
              selection = list(mode = "single", target = "row"),
              options = list(iDisplayLength = 20, ordering = FALSE, autowidth = TRUE),
              rownames = FALSE
            )
          })
        })
      }
    }
  )
}
