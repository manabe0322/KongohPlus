est_par_result_ui <- function(id){
  ns <- NS(id)

  tabPanel("Result",
           br(),
           tabsetPanel(id = ns("tabset_est_par_result"),
                       type = "tabs",
                       tabPanel("AE",
                                sidebarLayout(
                                  sidebarPanel(width = 2,
                                               actionButton(ns("act_trace_ae"), label = "Trace plot"),
                                               br(),
                                               br(),
                                               br(),
                                               actionButton(ns("act_plot_ae"), label = "Data plot")
                                  ),
                                  mainPanel(width = 10,
                                            dataTableOutput(ns("dt_result_ae"))
                                  )
                                )
                       ),
                       tabPanel("BSR",
                                sidebarLayout(
                                  sidebarPanel(width = 2,
                                               actionButton(ns("change_model_bsr"), label = "Change model"),
                                               br(),
                                               br(),
                                               br(),
                                               actionButton(ns("act_trace_bsr"), label = "Trace plot"),
                                               br(),
                                               br(),
                                               br(),
                                               actionButton(ns("act_plot_bsr"), label = "Data plot")
                                  ),
                                  mainPanel(width = 10,
                                            dataTableOutput(ns("dt_result_bsr"))
                                  )
                                )
                       ),
                       tabPanel("FSR",
                                sidebarLayout(
                                  sidebarPanel(width = 2,
                                               actionButton(ns("change_model_fsr"), label = "Change model"),
                                               br(),
                                               br(),
                                               br(),
                                               actionButton(ns("act_trace_fsr"), label = "Trace plot"),
                                               br(),
                                               br(),
                                               br(),
                                               actionButton(ns("act_plot_fsr"), label = "Data plot")
                                  ),
                                  mainPanel(width = 10,
                                            dataTableOutput(ns("dt_result_fsr"))
                                  )
                                )
                       ),
                       tabPanel("DSR",
                                sidebarLayout(
                                  sidebarPanel(width = 2,
                                               actionButton(ns("change_model_dsr"), label = "Change model"),
                                               br(),
                                               br(),
                                               br(),
                                               actionButton(ns("act_trace_dsr"), label = "Trace plot"),
                                               br(),
                                               br(),
                                               br(),
                                               actionButton(ns("act_plot_dsr"), label = "Data plot")
                                  ),
                                  mainPanel(width = 10,
                                            dataTableOutput(ns("dt_result_dsr"))
                                  )
                                )
                       ),
                       tabPanel("M2SR",
                                sidebarLayout(
                                  sidebarPanel(width = 2,
                                               actionButton(ns("change_model_m2sr"), label = "Change model"),
                                               br(),
                                               br(),
                                               br(),
                                               actionButton(ns("act_trace_m2sr"), label = "Trace plot"),
                                               br(),
                                               br(),
                                               br(),
                                               actionButton(ns("act_plot_m2sr"), label = "Data plot")
                                  ),
                                  mainPanel(width = 10,
                                            dataTableOutput(ns("dt_result_m2sr"))
                                  )
                                )
                       )
           )
  )
}

est_par_result_server <- function(id){
  moduleServer(
    id,
    function(input, output, session){
    }
  )
}
