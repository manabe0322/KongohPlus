#' lr_result_ui
#'
#' @description The function to create the ui module for results of the likelihood ratio
lr_result_ui <- function(id){
  ns <- NS(id)

  tabPanel("Result",
           br(),
           fluidRow(
             column(2,
                    wellPanel(
                      uiOutput(ns("hp")),
                      uiOutput(ns("noc_lr_hp")),
                      uiOutput(ns("hd")),
                      uiOutput(ns("noc_lr_hd"))
                    )
             ),
             column(2,
                    dataTableOutput(ns("dt_lr")),
                    br()
             ),
             column(4,
                    wellPanel(
                      h4("Parameters and diagnostics"),
                      fluidRow(
                        column(6,
                               uiOutput(ns("select_noc_diag"))
                        ),
                        column(6,
                               actionButton(ns("act_traceplot_lp"), label = "Trace plot")
                        )
                      ),
                      rHandsontableOutput(ns("dt_diag")),
                      br()
                    ),
                    wellPanel(
                      h4("Warnings"),
                      htmlOutput(ns("warning_message"))
                    )
             ),
             column(4,
                    wellPanel(
                      h4("Probabilistic genotyping"),
                      fluidRow(
                        column(6,
                        ),
                        column(6,
                        )
                      )
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
        dt_lr <- data_lr$dt_lr

        data_deconvo_hp <- rv_lr$data_deconvo_hp
        hnc_hp <- data_deconvo_hp$hnc
        dt_diag_hp <- data_deconvo_hp$dt_diag

        data_deconvo_hd <- rv_lr$data_deconvo_hd
        hnc_hd <- data_deconvo_hd$hnc
        dt_diag_hd <- data_deconvo_hd$dt_diag

        output$noc_lr_hp <- renderUI({
          disabled(numericInput(session$ns("noc_lr_hp"), label = "Number of contributors", value = hnc_hp))
        })

        output$noc_lr_hd <- renderUI({
          disabled(numericInput(session$ns("noc_lr_hd"), label = "Number of contributors", value = hnc_hd))
        })

        ## Parameters and diagnostics

        output$select_noc_diag <- renderUI({
          cand_noc <- sort(unique(c(hnc_hp, hnc_hd)))
          selectInput(session$ns("select_noc_diag"), label = "Number of contributors", choices = cand_noc, selected = cand_noc[1], multiple = FALSE, selectize = FALSE)
        })

        observeEvent(input$select_noc_diag, {
          noc_diag <- input$select_noc_diag
          if(noc_diag == hnc_hp){
            dt_diag <- dt_diag_hp
          }else{
            dt_diag <- dt_diag_hd
          }

          output$dt_diag <- renderRHandsontable({
            rhandsontable(dt_diag, rowHeaders = FALSE) %>%
              hot_table(fillHandle = FALSE, rowHeights = 30) %>%
              hot_col("Parameter", type = "text", readOnly = TRUE, width = 175) %>%
              hot_col("Contributor", type = "text", readOnly = TRUE, width = 90) %>%
              hot_col("Marker", type = "text", readOnly = TRUE, width = 100) %>%
              hot_col("Estimated value", type = "numeric", readOnly = TRUE, width = 125) %>%
              hot_col("Rhat", type = "numeric", readOnly = TRUE, width = 50)
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
