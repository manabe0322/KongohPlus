#' deconvo_pg_ui
#'
#' @description The function to create the ui module for probabilistic genotyping
deconvo_pg_ui <- function(id){
  ns <- NS(id)

  tabPanel("Probabilistic genotyping",
           br(),
           fluidRow(
             column(6,
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
             column(6)
           )
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

#' deconvo_pg_server
#'
#' @description The function to create the server module for probabilistic genotyping
deconvo_pg_server <- function(id, rv_deconvo){
  moduleServer(
    id,
    function(input, output, session){
      data_deconvo <- rv_deconvo$data_deconvo
      if(!is.null(data_deconvo)){
        #output$select_meth_pg <- renderUI({
        #  selectInput(session$ns("select_meth_pg"), label = "Method", choices = c("Per combination", "Per contributor"), selected = "Per combination", selectize = FALSE)
        #})

        #output$select_noc_pg <- renderUI({
        #  disabled(selectInput(session$ns("select_noc_pg"), label = "Number of contributors", choices = hnc, selected = hnc[1], selectize = FALSE))
        #})

        #output$select_mk_pg <- renderUI({
        #  selectInput(session$ns("select_mk_pg"), label = "Marker", choices = markers_use, selected = markers_use[1], selectize = FALSE)
        #})

        #output$select_contrib_pg <- renderUI({
        #  cand_contrib <- paste0("U", 1:hnc)
        #  selectInput(session$ns("select_contrib_pg"), label = "Contributor", choices = cand_contrib, selected = cand_contrib[1], selectize = FALSE)
        #})

        #observeEvent(input$select_meth_pg, ignoreInit = TRUE, {
        #  rv_tmp$update <- TRUE
        #})

        #observeEvent(input$select_mk_pg, ignoreInit = TRUE, {
        #  rv_tmp$update <- TRUE
        #})

        #observeEvent(input$select_contrib_pg, ignoreInit = TRUE, {
        #  rv_tmp$update <- TRUE
        #})

        #observe({
        #  if(rv_tmp$update){
        #    meth_pg <- input$select_meth_pg
        #    if(meth_pg == "Per combination"){
        #      disable("select_contrib_pg")
        #      dt_pg_select <- dt_pg_list[[input$select_mk_pg]]

        #      output$dt_pg <- renderDataTable(server = FALSE, {
        #        datatable(
        #          dt_pg_select,
        #          selection = "none",
        #          options = list(iDisplayLength = 10, ordering = FALSE, autowidth = TRUE, searching = FALSE,
        #                         columnDefs = list(list(targets = c(hnc, hnc + 1, hnc + 2), visible = FALSE))),
        #          rownames = FALSE
        #        ) %>%
        #          formatSignif(columns = c("Probability"), digits = 3) %>%
        #          formatStyle(names(dt_pg_select), lineHeight = "80%")
        #      })
        #    }else{
        #      enable("select_contrib_pg")

        #      dt_pg_select <- dt_pg_list[[input$select_mk_pg]]
        #      aggregated_data <- aggregate_pg_contrib(dt_pg_select, input$select_contrib_pg)

        #      output$dt_pg <- renderDataTable(server = FALSE, {
        #        datatable(
        #          aggregated_data,
        #          selection = "none",
        #          options = list(iDisplayLength = 10, ordering = FALSE, autowidth = TRUE, searching = FALSE),
        #          rownames = FALSE
        #        ) %>%
        #          formatSignif(columns = c("Probability"), digits = 3) %>%
        #          formatStyle(names(aggregated_data), lineHeight = "80%")
        #      })
        #    }

        #    rv_tmp$update <- FALSE
        #  }
        #})
      }
    }
  )
}
