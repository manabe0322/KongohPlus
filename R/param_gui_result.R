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

est_par_result_server <- function(id, react_results_est_par){
  moduleServer(
    id,
    function(input, output, session){
      dt_result_ae <- NULL
      if(!is.null(react_results_est_par())){
        results_est_par <- react_results_est_par()
        markers_par <- results_est_par$markers_par
        ae_median_myu_all_mk <- results_est_par$ae_median_myu_all_mk
        ae_median_sigma_all_mk <- results_est_par$ae_median_sigma_all_mk
        ae_rhat_myu_all_mk <- results_est_par$ae_rhat_myu_all_mk
        ae_rhat_sigma_all_mk <- results_est_par$ae_rhat_sigma_all_mk
        ae_traceplot_all_mk <- results_est_par$ae_traceplot_all_mk
        dt_result_ae <- data.table(Marker = markers_par,
                                   Myu = ae_median_myu_all_mk,
                                   Sigma = ae_median_sigma_all_mk,
                                   Rhat_Myu = ae_rhat_myu_all_mk,
                                   Rhat_Sigma = ae_rhat_sigma_all_mk)
      }

      output$dt_result_ae <- renderDataTable(server = FALSE, {
        datatable(dt_result_ae)
      })

#      ae_median_myu_all_mk <- results_est_par$ae_median_myu_all_mk
#      ae_median_sigma_all_mk <- results_est_par$ae_median_sigma_all_mk
#      ae_rhat_myu_all_mk <- results_est_par$ae_rhat_myu_all_mk
#      ae_rhat_sigma_all_mk <- results_est_par$ae_rhat_sigma_all_mk
#      ae_traceplot_all_mk <- results_est_par$ae_traceplot_all_mk
#      bsr_allele_median_slope_all_mk <- results_est_par$bsr_allele_median_slope_all_mk
#      bsr_allele_median_intercept_all_mk <- results_est_par$bsr_allele_median_intercept_all_mk
#      bsr_allele_median_sigma_all_mk <- results_est_par$bsr_allele_median_sigma_all_mk
#      bsr_allele_rhat_slope_all_mk <- results_est_par$bsr_allele_rhat_slope_all_mk
#      bsr_allele_rhat_intercept_all_mk <- results_est_par$bsr_allele_rhat_intercept_all_mk
#      bsr_allele_rhat_sigma_all_mk <- results_est_par$bsr_allele_rhat_sigma_all_mk
#      bsr_allele_trace <- results_est_par$bsr_allele_trace
#      bsr_lus_median_slope_all_mk <- results_est_par$bsr_lus_median_slope_all_mk
#      bsr_lus_median_intercept_all_mk <- results_est_par$bsr_lus_median_intercept_all_mk
#      bsr_lus_median_sigma_all_mk <- results_est_par$bsr_lus_median_sigma_all_mk
#      bsr_lus_rhat_slope_all_mk <- results_est_par$bsr_lus_rhat_slope_all_mk
#      bsr_lus_rhat_intercept_all_mk <- results_est_par$bsr_lus_rhat_intercept_all_mk
#      bsr_lus_rhat_sigma_all_mk <- results_est_par$bsr_lus_rhat_sigma_all_mk
#      bsr_lus_trace <- results_est_par$bsr_lus_trace
#      bsr_multi_median_slope_all_mk <- results_est_par$bsr_multi_median_slope_all_mk
#      bsr_multi_median_intercept_all_mk <- results_est_par$bsr_multi_median_intercept_all_mk
#      bsr_multi_median_sigma_all_mk <- results_est_par$bsr_multi_median_sigma_all_mk
#      bsr_multi_min_us_all_mk <- results_est_par$bsr_multi_min_us_all_mk
#      bsr_multi_rhat_slope_all_mk <- results_est_par$bsr_multi_rhat_slope_all_mk
#      bsr_multi_rhat_intercept_all_mk <- results_est_par$bsr_multi_rhat_intercept_all_mk
#      bsr_multi_rhat_sigma_all_mk <- results_est_par$bsr_multi_rhat_sigma_all_mk
#      bsr_multi_trace <- results_est_par$bsr_multi_trace
#      bsr_best_model_all_mk <- results_est_par$bsr_best_model_all_mk
#      bsr_best_median_myu_all_mk <- results_est_par$bsr_best_median_myu_all_mk
#      bsr_best_median_slope_all_mk <- results_est_par$bsr_best_median_slope_all_mk
#      bsr_best_median_intercept_all_mk <- results_est_par$bsr_best_median_intercept_all_mk
#      bsr_best_median_sigma_all_mk <- results_est_par$bsr_best_median_sigma_all_mk
#      bsr_best_min_us_all_mk <- results_est_par$bsr_best_min_us_all_mk
#      bsr_best_rhat_myu_all_mk <- results_est_par$bsr_best_rhat_myu_all_mk
#      bsr_best_rhat_slope_all_mk <- results_est_par$bsr_best_rhat_slope_all_mk
#      bsr_best_rhat_intercept_all_mk <- results_est_par$bsr_best_rhat_intercept_all_mk
#      bsr_best_rhat_sigma_all_mk <- results_est_par$bsr_best_rhat_sigma_all_mk
#      bsr_best_trace <- results_est_par$bsr_best_trace
#      fsr_allele_median_slope_all_mk <- results_est_par$fsr_allele_median_slope_all_mk
#      fsr_allele_median_intercept_all_mk <- results_est_par$fsr_allele_median_intercept_all_mk
#      fsr_allele_median_sigma_all_mk <- results_est_par$fsr_allele_median_sigma_all_mk
#      fsr_allele_rhat_slope_all_mk <- results_est_par$fsr_allele_rhat_slope_all_mk
#      fsr_allele_rhat_intercept_all_mk <- results_est_par$fsr_allele_rhat_intercept_all_mk
#      fsr_allele_rhat_sigma_all_mk <- results_est_par$fsr_allele_rhat_sigma_all_mk
#      fsr_allele_trace <- results_est_par$fsr_allele_trace
#      fsr_lus_median_slope_all_mk <- results_est_par$fsr_lus_median_slope_all_mk
#      fsr_lus_median_intercept_all_mk <- results_est_par$fsr_lus_median_intercept_all_mk
#      fsr_lus_median_sigma_all_mk <- results_est_par$fsr_lus_median_sigma_all_mk
#      fsr_lus_rhat_slope_all_mk <- results_est_par$fsr_lus_rhat_slope_all_mk
#      fsr_lus_rhat_intercept_all_mk <- results_est_par$fsr_lus_rhat_intercept_all_mk
#      fsr_lus_rhat_sigma_all_mk <- results_est_par$fsr_lus_rhat_sigma_all_mk
#      fsr_lus_trace <- results_est_par$fsr_lus_trace
#      fsr_multi_median_slope_all_mk <- results_est_par$fsr_multi_median_slope_all_mk
#      fsr_multi_median_intercept_all_mk <- results_est_par$fsr_multi_median_intercept_all_mk
#      fsr_multi_median_sigma_all_mk <- results_est_par$fsr_multi_median_sigma_all_mk
#      fsr_multi_min_us_all_mk <- results_est_par$fsr_multi_min_us_all_mk
#      fsr_multi_rhat_slope_all_mk <- results_est_par$fsr_multi_rhat_slope_all_mk
#      fsr_multi_rhat_intercept_all_mk <- results_est_par$fsr_multi_rhat_intercept_all_mk
#      fsr_multi_rhat_sigma_all_mk <- results_est_par$fsr_multi_rhat_sigma_all_mk
#      fsr_multi_trace <- results_est_par$fsr_multi_trace
#      fsr_best_model_all_mk <- results_est_par$fsr_best_model_all_mk
#      fsr_best_median_myu_all_mk <- results_est_par$fsr_best_median_myu_all_mk
#      fsr_best_median_slope_all_mk <- results_est_par$fsr_best_median_slope_all_mk
#      fsr_best_median_intercept_all_mk <- results_est_par$fsr_best_median_intercept_all_mk
#      fsr_best_median_sigma_all_mk <- results_est_par$fsr_best_median_sigma_all_mk
#      fsr_best_min_us_all_mk <- results_est_par$fsr_best_min_us_all_mk
#      fsr_best_rhat_myu_all_mk <- results_est_par$fsr_best_rhat_myu_all_mk
#      fsr_best_rhat_slope_all_mk <- fsr_best_rhat_slope_all_mk
#      fsr_best_rhat_intercept_all_mk <- fsr_best_rhat_intercept_all_mk
#      fsr_best_rhat_sigma_all_mk <- fsr_best_rhat_sigma_all_mk
#      fsr_best_trace <- fsr_best_trace
#      dsr_allele_median_slope_all_mk <- dsr_allele_median_slope_all_mk
#      dsr_allele_median_intercept_all_mk <- dsr_allele_median_intercept_all_mk
#      dsr_allele_median_sigma_all_mk <- dsr_allele_median_sigma_all_mk
#      dsr_allele_rhat_slope_all_mk <- dsr_allele_rhat_slope_all_mk
#      dsr_allele_rhat_intercept_all_mk <- dsr_allele_rhat_intercept_all_mk
#      dsr_allele_rhat_sigma_all_mk <- dsr_allele_rhat_sigma_all_mk
#      dsr_allele_trace <- dsr_allele_trace
#      dsr_lus_median_slope_all_mk <- dsr_lus_median_slope_all_mk
#      dsr_lus_median_intercept_all_mk <- dsr_lus_median_intercept_all_mk
#      dsr_lus_median_sigma_all_mk <- dsr_lus_median_sigma_all_mk
#      dsr_lus_rhat_slope_all_mk <- dsr_lus_rhat_slope_all_mk
#      dsr_lus_rhat_intercept_all_mk <- dsr_lus_rhat_intercept_all_mk
#      dsr_lus_rhat_sigma_all_mk <- dsr_lus_rhat_sigma_all_mk
#      dsr_lus_trace <- dsr_lus_trace
#      dsr_multi_median_slope_all_mk <- dsr_multi_median_slope_all_mk
#      dsr_multi_median_intercept_all_mk <- dsr_multi_median_intercept_all_mk
#      dsr_multi_median_sigma_all_mk <- dsr_multi_median_sigma_all_mk
#      dsr_multi_min_us_all_mk <- dsr_multi_min_us_all_mk
#      dsr_multi_rhat_slope_all_mk <- dsr_multi_rhat_slope_all_mk
#      dsr_multi_rhat_intercept_all_mk <- dsr_multi_rhat_intercept_all_mk
#      dsr_multi_rhat_sigma_all_mk <- dsr_multi_rhat_sigma_all_mk
#      dsr_multi_trace <- dsr_multi_trace
#      dsr_best_model_all_mk <- dsr_best_model_all_mk
#      dsr_best_median_myu_all_mk <- dsr_best_median_myu_all_mk
#      dsr_best_median_slope_all_mk <- dsr_best_median_slope_all_mk
#      dsr_best_median_intercept_all_mk <- dsr_best_median_intercept_all_mk
#      dsr_best_median_sigma_all_mk <- dsr_best_median_sigma_all_mk
#      dsr_best_min_us_all_mk <- dsr_best_min_us_all_mk
#      dsr_best_rhat_myu_all_mk <- dsr_best_rhat_myu_all_mk
#      dsr_best_rhat_slope_all_mk <- dsr_best_rhat_slope_all_mk
#      dsr_best_rhat_intercept_all_mk <- dsr_best_rhat_intercept_all_mk
#      dsr_best_rhat_sigma_all_mk <- dsr_best_rhat_sigma_all_mk
#      dsr_best_trace <- dsr_best_trace
#      m2sr_best_model_all_mk <- m2sr_best_model_all_mk
#      m2sr_best_median_myu_all_mk <- m2sr_best_median_myu_all_mk
#      m2sr_best_median_sigma_all_mk <- m2sr_best_median_sigma_all_mk
#      m2sr_best_rhat_myu_all_mk <- m2sr_best_rhat_myu_all_mk
#      m2sr_best_rhat_sigma_all_mk <- m2sr_best_rhat_sigma_all_mk
#      m2sr_best_trace <- m2sr_best_trace
    }
  )
}
