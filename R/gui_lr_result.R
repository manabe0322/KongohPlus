#' lr_result_ui
#'
#' @description The function to create the ui module for results of the likelihood ratio
lr_result_ui <- function(id){
  ns <- NS(id)

  tabPanel("Result",
           br(),
           fluidRow(
             column(3,
                    wellPanel(
                      h4("Parameters for contributors"),
                      br(),
                      dataTableOutput(ns("dt_mpdeg_hp")),
                      br(),
                      dataTableOutput(ns("dt_mpdeg_hd")),
                    ),
                    br(),
                    actionButton(ns("act_pardiag"), label = "Parameters and diagnostics"),
                    br(),
                    br(),
                    actionButton(ns("act_traceplot_lp"), label = "Trace plot"),
                    br(),
                    br(),
                    downloadButton(ns("download_report"), "Report")
             ),
             column(6,
                    wellPanel(
                      h4("Probabilistic genotyping"),
                      br(),
                      uiOutput(ns("hyp_pg")),
                      br(),
                      uiOutput(ns("mk_pg")),
                      br(),
                      dataTableOutput(ns("dt_pg")),
                      br()
                    )
             ),
             column(3,
                    wellPanel(
                      h4("Likelihood ratio"),
                      br(),
                      dataTableOutput(ns("dt_lr")),
                      br()
                    )
             )
           )
  )
}

#' lr_result_server
#'
#' @description The function to create the server module for results of the likelihood ratio
lr_result_server <- function(id, rv_lr){
  moduleServer(
    id,
    function(input, output, session){
      data_lr <- rv_lr$data_lr
      if(!is.null(data_lr)){
        dt_mpdeg_hp <- data_lr$dt_mpdeg_hp
        dt_mpdeg_hd <- data_lr$dt_mpdeg_hd
        dt_pg_list_hp <- data_lr$dt_pg_list_hp
        dt_pg_list_hd <- data_lr$dt_pg_list_hd
        dt_lr <- data_lr$dt_lr

        hnc <- nrow(dt_mpdeg_hd)
        markers_use <- names(dt_pg_list_hd)

        rv_lr_result <- reactiveValues()
        rv_lr_result$dt_pg_display <- dt_pg_list_hp

        output$dt_mpdeg_hp <- renderDataTable(server = FALSE, {
          datatable(
            dt_mpdeg_hp,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_mpdeg_hp), ordering = FALSE, autowidth = TRUE),
            rownames = FALSE
          ) %>%
            formatSignif(columns = c("MP", "Deg"), digits = 3)
        })

        output$dt_mpdeg_hd <- renderDataTable(server = FALSE, {
          datatable(
            dt_mpdeg_hd,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_mpdeg_hd), ordering = FALSE, autowidth = TRUE),
            rownames = FALSE
          ) %>%
            formatSignif(columns = c("MP", "Deg"), digits = 3)
        })

        output$hyp_pg <- renderUI({
          selectInput(session$ns("hyp_pg"), label = "Hypothesis", choices = c("Hp", "Hd"), selected = "Hp", selectize = FALSE)
        })

        output$mk_pg <- renderUI({
          selectInput(session$ns("mk_pg"), label = "Marker", choices = markers_use, selected = markers_use[1], selectize = FALSE)
        })

        observeEvent(input$hyp_pg, {
          hyp_select <- input$hyp_pg
          if(hyp_select == "Hp"){
            rv_lr_result$dt_pg_display <- dt_pg_list_hp
          }else{
            rv_lr_result$dt_pg_display <- dt_pg_list_hd
          }

          mk_select <- input$mk_pg
          isolate(dt_pg_display <- rv_lr_result$dt_pg_display)

          output$dt_pg <- renderDataTable(server = FALSE, {
            datatable(
              dt_pg_display[[mk_select]],
              selection = list(mode = "single", target = "row"),
              options = list(iDisplayLength = 25, ordering = FALSE, autowidth = TRUE,
                             columnDefs = list(list(targets = c(hnc, hnc + 1, hnc + 2), visible = FALSE))),
              rownames = FALSE
            ) %>%
              formatSignif(columns = c("post_prob"), digits = 3)
          })
        })

        observeEvent(input$mk_pg, {
          mk_select <- input$mk_pg
          isolate(dt_pg_display <- rv_lr_result$dt_pg_display)

          output$dt_pg <- renderDataTable(server = FALSE, {
            datatable(
              dt_pg_display[[mk_select]],
              selection = list(mode = "single", target = "row"),
              options = list(iDisplayLength = 25, ordering = FALSE, autowidth = TRUE,
                             columnDefs = list(list(targets = c(hnc, hnc + 1, hnc + 2), visible = FALSE))),
              rownames = FALSE
            ) %>%
              formatSignif(columns = c("post_prob"), digits = 3)
          })
        })

        output$dt_lr <- renderDataTable(server = FALSE, {
          datatable(
            dt_lr,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_lr), ordering = FALSE, autowidth = TRUE,
                           columnDefs = list(list(targets = c(1, 2), visible = FALSE))),
            rownames = FALSE
          ) %>%
            formatSignif(columns = c("LR"), digits = 3)
        })
      }
    }
  )
}
