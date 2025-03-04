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
                      uiOutput(ns("pg"))
                    )
             ),
             column(6)
           )
  )
}

#' aggregate_pg_contrib
#'
#' @description The function to change the format of the pg table from per combination to per contributor
#' @param dt_pg The data.table of the probabilistic genotyping in one marker
#' @param contrib The name of the selected contributor
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
        markers_use <- data_deconvo$markers_use
        hnc_all <- data_deconvo$hnc_all
        data_hnc_list <- data_deconvo$data_hnc_list

        rv_deconvo_pg <- reactiveValues()
        rv_deconvo_pg$cand_contrib <- paste0("U", 1:hnc_all[1])
        rv_deconvo_pg$dt_pg_per_comb <- NULL
        rv_deconvo_pg$dt_pg_per_contrib <- NULL
        rv_deconvo_pg$error_message <- NULL

        output$select_meth_pg <- renderUI({
          selectInput(session$ns("select_meth_pg"), label = "Method", choices = c("Per combination", "Per contributor"), selected = "Per combination", selectize = FALSE)
        })

        output$select_noc_pg <- renderUI({
          selectInput(session$ns("select_noc_pg"), label = "Number of contributors", choices = hnc_all, selected = hnc_all[1], selectize = FALSE)
        })

        output$select_mk_pg <- renderUI({
          selectInput(session$ns("select_mk_pg"), label = "Marker", choices = markers_use, selected = markers_use[1], selectize = FALSE)
        })

        output$select_contrib_pg <- renderUI({
          selectInput(session$ns("select_contrib_pg"), label = "Contributor", choices = rv_deconvo_pg$cand_contrib, selected = rv_deconvo_pg$cand_contrib[1], selectize = FALSE)
        })

        react_input_pg <- reactive({
          list(input$select_meth_pg, input$select_noc_pg, input$select_mk_pg, input$select_contrib_pg)
        })

        observeEvent(react_input_pg(), ignoreInit = TRUE, {
          hnc <- as.numeric(input$select_noc_pg)
          rv_deconvo_pg$cand_contrib <- paste0("U", 1:hnc)

          index_hnc <- which(hnc_all == hnc)
          data_hnc <- data_hnc_list[[index_hnc]]
          rv_deconvo_pg$error_message <- data_hnc$error_message
          if(is.null(rv_deconvo_pg$error_message)){
            dt_pg_list <- data_hnc$dt_pg_list

            if(input$select_meth_pg == "Per combination"){
              rv_deconvo_pg$dt_pg_per_comb <- dt_pg_list[[input$select_mk_pg]]
              rv_deconvo_pg$dt_pg_per_contrib <- NULL
            }else{
              rv_deconvo_pg$dt_pg_per_comb <- NULL
              rv_deconvo_pg$dt_pg_per_contrib <- aggregate_pg_contrib(dt_pg_list[[input$select_mk_pg]], input$select_contrib_pg)
            }
          }else{
            rv_deconvo_pg$dt_pg_per_comb <- NULL
            rv_deconvo_pg$dt_pg_per_contrib <- NULL
          }
        })

        output$pg <- renderUI({
          dt_pg_per_comb <- rv_deconvo_pg$dt_pg_per_comb
          dt_pg_per_contrib <- rv_deconvo_pg$dt_pg_per_contrib
          error_message <- rv_deconvo_pg$error_message
          if(is.null(error_message)){
            if(is.null(dt_pg_per_contrib)){
              dataTableOutput(session$ns("dt_pg_per_comb"))
            }else{
              dataTableOutput(session$ns("dt_pg_per_contrib"))
            }
          }else{
            HTML(error_message)
          }
        })

        output$dt_pg_per_comb <- renderDataTable(server = FALSE, {
          dt_pg_per_comb <- rv_deconvo_pg$dt_pg_per_comb
          hnc <- as.numeric(input$select_noc_pg)
          if(!is.null(dt_pg_per_comb)){
            datatable(
              dt_pg_per_comb,
              selection = "none",
              options = list(iDisplayLength = 10, ordering = FALSE, autowidth = TRUE, searching = FALSE,
                             columnDefs = list(list(targets = c(hnc, hnc + 1, hnc + 2), visible = FALSE))),
              rownames = FALSE
            ) %>%
              formatSignif(columns = c("Probability"), digits = 3) %>%
              formatStyle(names(dt_pg_per_comb), lineHeight = "80%")
          }
        })

        output$dt_pg_per_contrib <- renderDataTable(server = FALSE, {
          dt_pg_per_contrib <- rv_deconvo_pg$dt_pg_per_contrib
          if(!is.null(dt_pg_per_contrib)){
            datatable(
              dt_pg_per_contrib,
              selection = "none",
              options = list(iDisplayLength = 10, ordering = FALSE, autowidth = TRUE, searching = FALSE),
              rownames = FALSE
            ) %>%
              formatSignif(columns = c("Probability"), digits = 3) %>%
              formatStyle(names(dt_pg_per_contrib), lineHeight = "80%")
          }
        })

        observeEvent(input$dt_pg_per_comb_rows_selected, {
          index_select <- input$dt_pg_per_comb_rows_selected

          if(is.null(index_select)){

          }else{

          }
        })

        observeEvent(input$dt_pg_per_contrib_rows_selected, {
          index_select <- input$dt_pg_per_contrib_rows_selected

          if(is.null(index_select)){

          }else{

          }
        })
      }
    }
  )
}
