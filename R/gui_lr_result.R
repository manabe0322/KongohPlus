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
             column(5,
                    dataTableOutput(ns("dt_lr")),
                    br()
             ),
             column(5,
                    wellPanel(
                      h4("Parameters and diagnostics"),
                      rHandsontableOutput(ns("dt_par")),
                      br(),
                      actionButton(ns("act_traceplot_lp"), label = "Trace plot")
                    ),
                    wellPanel(
                      h4("Warnings"),
                      htmlOutput(ns("warning_message"))
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

        #output$hp <- renderUI({
        #  checkboxGroupInput(session$ns("hp"), label = h4("Prosecutor hypothesis"), choices = c("A", "B"), selected = NULL)
        #})

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
