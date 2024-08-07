#' deconvo_pg_ui
#'
#' @description The function to create the ui module for probabilistic genotyping (deconvolution)
deconvo_pg_ui <- function(id){
  ns <- NS(id)

  tabPanel("Probabilistic genotyping",
           br(),
           fluidRow(
             column(6,
                    wellPanel(
                      h4("Probabilistic genotyping (per combination)"),
                      br(),
                      uiOutput(ns("mk_pg_comb")),
                      br(),
                      br(),
                      br(),
                      br(),
                      br(),
                      br(),
                      dataTableOutput(ns("dt_pg_comb")),
                      br()
                    )
             ),
             column(6,
                    wellPanel(
                      h4("Probabilistic genotyping (per contributor)"),
                      br(),
                      uiOutput(ns("mk_pg_contrib")),
                      br(),
                      uiOutput(ns("contrib_pg_contrib")),
                      br(),
                      dataTableOutput(ns("dt_pg_contrib")),
                      br()
                    )
             )
           )
  )
}

#' deconvo_pg_server
#'
#' @description The function to create the server module for probabilistic genotyping (deconvolution)
deconvo_pg_server <- function(id, rv_deconvo){
  moduleServer(
    id,
    function(input, output, session){
      data_deconvo <- rv_deconvo$data_deconvo
      if(!is.null(data_deconvo)){
        rv_tmp <- reactiveValues(updated = FALSE)

        hnc <- data_deconvo$hnc
        markers_use <- data_deconvo$markers_use
        dt_pg_list <- data_deconvo$dt_pg_list
        contributors <- names(dt_pg_list[[1]])[1:hnc]

        output$mk_pg_comb <- renderUI({
          selectInput(session$ns("mk_pg_comb"), label = "Marker", choices = markers_use, selected = markers_use[1], selectize = FALSE)
        })

        observeEvent(input$mk_pg_comb, {
          mk_select <- input$mk_pg_comb

          output$dt_pg_comb <- renderDataTable(server = FALSE, {
            datatable(
              dt_pg_list[[mk_select]],
              selection = list(mode = "single", target = "row"),
              options = list(iDisplayLength = 10, ordering = FALSE, autowidth = TRUE, searching = FALSE,
                             columnDefs = list(list(targets = c(hnc, hnc + 1, hnc + 2), visible = FALSE))),
              rownames = FALSE
            ) %>%
              formatSignif(columns = c("Probability"), digits = 3)
          })
        })

        output$mk_pg_contrib <- renderUI({
          selectInput(session$ns("mk_pg_contrib"), label = "Marker", choices = markers_use, selected = markers_use[1], selectize = FALSE)
        })

        output$contrib_pg_contrib <- renderUI({
          selectInput(session$ns("contrib_pg_contrib"), label = "Contributor", choices = contributors, selected = contributors[1], selectize = FALSE)
        })

        observeEvent(input$mk_pg_contrib, {
          rv_tmp$updated <- TRUE
        })

        observeEvent(input$contrib_pg_contrib, {
          rv_tmp$updated <- TRUE
        })

        observe({
          if(rv_tmp$updated){
            mk_select <- input$mk_pg_contrib
            dt_pg_select <- dt_pg_list[[mk_select]]
            contrib_select <- input$contrib_pg_contrib
            aggregated_data <- aggregate_pg_contrib(dt_pg_select, contrib_select)

            output$dt_pg_contrib <- renderDataTable(server = FALSE, {
              datatable(
                aggregated_data,
                selection = list(mode = "single", target = "row"),
                options = list(iDisplayLength = 10, ordering = FALSE, autowidth = TRUE, searching = FALSE),
                rownames = FALSE
              ) %>%
                formatSignif(columns = c("Probability"), digits = 3)
            })

            rv_tmp$updated <- FALSE
          }
        })
      }
    }
  )
}

aggregate_pg_contrib <- function(dt_pg, contrib){
  aggregated_data <- dt_pg %>%
    group_by_at(which(names(dt_pg) == contrib)) %>%
    summarise(Probability = sum(Probability))
  aggregated_data <- aggregated_data %>%
    arrange(desc(Probability))
  return(aggregated_data)
}

