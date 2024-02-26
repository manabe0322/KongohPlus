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
                      h4("Estimated parameter"),
                      br(),
                      dataTableOutput(ns("dt_mpdeg_hp")),
                      br(),
                      dataTableOutput(ns("dt_mpdeg_hd")),
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
                    )
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
        dt_lr <- data_lr$dt_lr

        output$dt_lr <- renderDataTable(server = FALSE, {
          datatable(
            dt_lr,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_lr), ordering = FALSE, autowidth = TRUE),
            rownames = FALSE
          ) %>%
            formatSignif(columns = c("LR"), digits = 3)
        })
      }
    }
  )
}
