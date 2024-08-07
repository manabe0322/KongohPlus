#' lr_pg_ui
#'
#' @description The function to create the ui module for probabilistic genotyping (LR)
lr_pg_ui <- function(id){
  ns <- NS(id)

  tabPanel("Probabilistic genotyping",
           br(),
           fluidRow(
             column(6,
                    wellPanel(
                      h4("Probabilistic genotyping (per combination)"),
                      br(),
                      uiOutput(ns("hyp_pg_comb")),
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
                      uiOutput(ns("hyp_pg_contrib")),
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

#' lr_pg_server
#'
#' @description The function to create the server module for probabilistic genotyping (LR)
lr_pg_server <- function(id, rv_lr){
  moduleServer(
    id,
    function(input, output, session){
      data_lr <- rv_lr$data_lr
      if(!is.null(data_lr)){
        rv_tmp <- reactiveValues()
        rv_tmp$updated_comb <- FALSE
        rv_tmp$updated_contrib <- FALSE
        rv_tmp$contributors <- NULL

        dt_pg_list_hp <- data_lr$dt_pg_list_hp
        data_deconvo_hp <- rv_lr$data_deconvo_hp
        hnc_hp <- data_deconvo_hp$hnc
        rv_tmp$contributors <- contributors_hp <- names(dt_pg_list_hp[[1]])[1:hnc_hp]
        markers_use <- data_deconvo_hp$markers_use

        dt_pg_list_hd <- data_lr$dt_pg_list_hd
        data_deconvo_hd <- rv_lr$data_deconvo_hd
        hnc_hd <- data_deconvo_hd$hnc
        contributors_hd <- names(dt_pg_list_hd[[1]])[1:hnc_hd]

        output$hyp_pg_comb <- renderUI({
          selectInput(session$ns("hyp_pg_comb"), label = "Hypothesis", choices = c("Hp", "Hd"), selected = "Hp", selectize = FALSE)
        })

        output$mk_pg_comb <- renderUI({
          selectInput(session$ns("mk_pg_comb"), label = "Marker", choices = markers_use, selected = markers_use[1], selectize = FALSE)
        })

        observeEvent(input$hyp_pg_comb, {
          rv_tmp$updated_comb <- TRUE
        })

        observeEvent(input$mk_pg_comb, {
          rv_tmp$updated_comb <- TRUE
        })

        observe({
          if(rv_tmp$updated_comb){
            if(input$hyp_pg_comb == "Hp"){
              dt_pg_display <- dt_pg_list_hp
              hnc <- hnc_hp
            }else{
              dt_pg_display <- dt_pg_list_hd
              hnc <- hnc_hd
            }

            output$dt_pg_comb <- renderDataTable(server = FALSE, {
              datatable(
                dt_pg_display[[input$mk_pg_comb]],
                selection = list(mode = "single", target = "row"),
                options = list(iDisplayLength = 10, ordering = FALSE, autowidth = TRUE, searching = FALSE,
                               columnDefs = list(list(targets = c(hnc, hnc + 1, hnc + 2), visible = FALSE))),
                rownames = FALSE
              ) %>%
                formatSignif(columns = c("Probability"), digits = 3)
            })

            rv_tmp$updated_comb <- FALSE
          }
        })

        output$hyp_pg_contrib <- renderUI({
          selectInput(session$ns("hyp_pg_contrib"), label = "Hypothesis", choices = c("Hp", "Hd"), selected = "Hp", selectize = FALSE)
        })

        output$mk_pg_contrib <- renderUI({
          selectInput(session$ns("mk_pg_contrib"), label = "Marker", choices = markers_use, selected = markers_use[1], selectize = FALSE)
        })

        output$contrib_pg_contrib <- renderUI({
          selectInput(session$ns("contrib_pg_contrib"), label = "Contributor", choices = rv_tmp$contributors, selected = rv_tmp$contributors[1], selectize = FALSE)
        })

        observeEvent(input$hyp_pg_contrib, {
          if(input$hyp_pg_contrib == "Hp"){
            rv_tmp$contributors <- contributors_hp
          }else{
            rv_tmp$contributors <- contributors_hd
          }
          rv_tmp$updated_contrib <- TRUE
        })

        observeEvent(input$mk_pg_contrib, {
          rv_tmp$updated_contrib <- TRUE
        })

        observeEvent(input$contrib_pg_contrib, {
          rv_tmp$updated_contrib <- TRUE
        })

        observe({
          if(rv_tmp$updated_contrib){
            if(input$hyp_pg_contrib == "Hp"){
              dt_pg_display <- dt_pg_list_hp
              hnc <- hnc_hp
            }else{
              dt_pg_display <- dt_pg_list_hd
              hnc <- hnc_hd
            }

            dt_pg_select <- dt_pg_display[[input$mk_pg_contrib]]
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

            rv_tmp$updated_contrib <- FALSE
          }
        })
      }
    }
  )
}
