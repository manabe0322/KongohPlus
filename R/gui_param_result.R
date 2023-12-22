est_par_result_ui <- function(id){
  ns <- NS(id)

  tabPanel("Result",
           br(),
           tabsetPanel(id = ns("tabset_est_par_result"),
                       type = "tabs",
                       tabPanel("AE",
                                br(),
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
                                br(),
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
                                br(),
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
                                br(),
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
                                br(),
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

est_par_result_server <- function(id, rv_output){
  moduleServer(
    id,
    function(input, output, session){
      results_est_par <- rv_output$results_est_par

      rv_results <- reactiveValues()
      rv_results$dt_result_ae <- NULL
      rv_results$dt_result_bsr <- NULL
      rv_results$dt_result_fsr <- NULL
      rv_results$dt_result_dsr <- NULL
      rv_results$dt_result_m2sr <- NULL
      if(is.null(results_est_par)){
        markers_par <- NULL

        ae_data <- NULL
#        ae_median_myu_all_mk <- NULL
#        ae_median_sigma_all_mk <- NULL
#        ae_rhat_myu_all_mk <- NULL
#        ae_rhat_sigma_all_mk <- NULL
        ae_trace_all_mk <- NULL

        bsr_data_all_mk <- NULL
#        bsr_allele_median_slope_all_mk <- NULL
#        bsr_allele_median_intercept_all_mk <- NULL
#        bsr_allele_median_sigma_all_mk <- NULL
#        bsr_allele_rhat_slope_all_mk <- NULL
#        bsr_allele_rhat_intercept_all_mk <- NULL
#        bsr_allele_rhat_sigma_all_mk <- NULL
#        bsr_allele_trace <- NULL
#        bsr_lus_median_slope_all_mk <- NULL
#        bsr_lus_median_intercept_all_mk <- NULL
#        bsr_lus_median_sigma_all_mk <- NULL
#        bsr_lus_rhat_slope_all_mk <- NULL
#        bsr_lus_rhat_intercept_all_mk <- NULL
#        bsr_lus_rhat_sigma_all_mk <- NULL
#        bsr_lus_trace <- NULL
#        bsr_multi_median_slope_all_mk <- NULL
#        bsr_multi_median_intercept_all_mk <- NULL
#        bsr_multi_median_sigma_all_mk <- NULL
#        bsr_multi_min_us_all_mk <- NULL
#        bsr_multi_rhat_slope_all_mk <- NULL
#        bsr_multi_rhat_intercept_all_mk <- NULL
#        bsr_multi_rhat_sigma_all_mk <- NULL
#        bsr_multi_trace <- NULL
#        bsr_best_model_all_mk <- NULL
#        bsr_best_median_myu_all_mk <- NULL
#        bsr_best_median_slope_all_mk <- NULL
#        bsr_best_median_intercept_all_mk <- NULL
#        bsr_best_median_sigma_all_mk <- NULL
#        bsr_best_min_us_all_mk <- NULL
#        bsr_best_rhat_myu_all_mk <- NULL
#        bsr_best_rhat_slope_all_mk <- NULL
#        bsr_best_rhat_intercept_all_mk <- NULL
#        bsr_best_rhat_sigma_all_mk <- NULL
        bsr_best_trace_all_mk <- NULL

        fsr_data_all_mk <- NULL
#        fsr_allele_median_slope_all_mk <- NULL
#        fsr_allele_median_intercept_all_mk <- NULL
#        fsr_allele_median_sigma_all_mk <- NULL
#        fsr_allele_rhat_slope_all_mk <- NULL
#        fsr_allele_rhat_intercept_all_mk <- NULL
#        fsr_allele_rhat_sigma_all_mk <- NULL
#        fsr_allele_trace <- NULL
#        fsr_lus_median_slope_all_mk <- NULL
#        fsr_lus_median_intercept_all_mk <- NULL
#        fsr_lus_median_sigma_all_mk <- NULL
#        fsr_lus_rhat_slope_all_mk <- NULL
#        fsr_lus_rhat_intercept_all_mk <- NULL
#        fsr_lus_rhat_sigma_all_mk <- NULL
#        fsr_lus_trace <- NULL
#        fsr_multi_median_slope_all_mk <- NULL
#        fsr_multi_median_intercept_all_mk <- NULL
#        fsr_multi_median_sigma_all_mk <- NULL
#        fsr_multi_min_us_all_mk <- NULL
#        fsr_multi_rhat_slope_all_mk <- NULL
#        fsr_multi_rhat_intercept_all_mk <- NULL
#        fsr_multi_rhat_sigma_all_mk <- NULL
#        fsr_multi_trace <- NULL
#        fsr_best_model_all_mk <- NULL
#        fsr_best_median_myu_all_mk <- NULL
#        fsr_best_median_slope_all_mk <- NULL
#        fsr_best_median_intercept_all_mk <- NULL
#        fsr_best_median_sigma_all_mk <- NULL
#        fsr_best_min_us_all_mk <- NULL
#        fsr_best_rhat_myu_all_mk <- NULL
#        fsr_best_rhat_slope_all_mk <- NULL
#        fsr_best_rhat_intercept_all_mk <- NULL
#        fsr_best_rhat_sigma_all_mk <- NULL
        fsr_best_trace_all_mk <- NULL

        dsr_data_all_mk <- NULL
#        dsr_allele_median_slope_all_mk <- NULL
#        dsr_allele_median_intercept_all_mk <- NULL
#        dsr_allele_median_sigma_all_mk <- NULL
#        dsr_allele_rhat_slope_all_mk <- NULL
#        dsr_allele_rhat_intercept_all_mk <- NULL
#        dsr_allele_rhat_sigma_all_mk <- NULL
#        dsr_allele_trace <- NULL
#        dsr_lus_median_slope_all_mk <- NULL
#        dsr_lus_median_intercept_all_mk <- NULL
#        dsr_lus_median_sigma_all_mk <- NULL
#        dsr_lus_rhat_slope_all_mk <- NULL
#        dsr_lus_rhat_intercept_all_mk <- NULL
#        dsr_lus_rhat_sigma_all_mk <- NULL
#        dsr_lus_trace <- NULL
#        dsr_multi_median_slope_all_mk <- NULL
#        dsr_multi_median_intercept_all_mk <- NULL
#        dsr_multi_median_sigma_all_mk <- NULL
#        dsr_multi_min_us_all_mk <- NULL
#        dsr_multi_rhat_slope_all_mk <- NULL
#        dsr_multi_rhat_intercept_all_mk <- NULL
#        dsr_multi_rhat_sigma_all_mk <- NULL
#        dsr_multi_trace <- NULL
#        dsr_best_model_all_mk <- NULL
#        dsr_best_median_myu_all_mk <- NULL
#        dsr_best_median_slope_all_mk <- NULL
#        dsr_best_median_intercept_all_mk <- NULL
#        dsr_best_median_sigma_all_mk <- NULL
#        dsr_best_min_us_all_mk <- NULL
#        dsr_best_rhat_myu_all_mk <- NULL
#        dsr_best_rhat_slope_all_mk <- NULL
#        dsr_best_rhat_intercept_all_mk <- NULL
#        dsr_best_rhat_sigma_all_mk <- NULL
        dsr_best_trace_all_mk <- NULL

        m2sr_data_all_mk <- NULL
#        m2sr_best_model_all_mk <- NULL
#        m2sr_best_median_myu_all_mk <- NULL
#        m2sr_best_median_sigma_all_mk <- NULL
#        m2sr_best_rhat_myu_all_mk <- NULL
#        m2sr_best_rhat_sigma_all_mk <- NULL
        m2sr_best_trace_all_mk <- NULL
      }else{
        markers_par <- results_est_par$markers_par

        ae_data <- results_est_par$ae_data
        ae_median_myu_all_mk <- results_est_par$ae_median_myu_all_mk
        ae_median_sigma_all_mk <- results_est_par$ae_median_sigma_all_mk
        ae_rhat_myu_all_mk <- results_est_par$ae_rhat_myu_all_mk
        ae_rhat_sigma_all_mk <- results_est_par$ae_rhat_sigma_all_mk
        ae_trace_all_mk <- results_est_par$ae_trace_all_mk
        rv_results$dt_result_ae <- data.table(Marker = markers_par,
                                              Myu = signif(ae_median_myu_all_mk, 3),
                                              Sigma = signif(ae_median_sigma_all_mk, 3),
                                              Rhat_Myu = signif(ae_rhat_myu_all_mk, 3),
                                              Rhat_Sigma = signif(ae_rhat_sigma_all_mk, 3)
                                   )

        bsr_data_all_mk <- results_est_par$bsr_data_all_mk
#        bsr_allele_median_slope_all_mk <- results_est_par$bsr_allele_median_slope_all_mk
#        bsr_allele_median_intercept_all_mk <- results_est_par$bsr_allele_median_intercept_all_mk
#        bsr_allele_median_sigma_all_mk <- results_est_par$bsr_allele_median_sigma_all_mk
#        bsr_allele_rhat_slope_all_mk <- results_est_par$bsr_allele_rhat_slope_all_mk
#        bsr_allele_rhat_intercept_all_mk <- results_est_par$bsr_allele_rhat_intercept_all_mk
#        bsr_allele_rhat_sigma_all_mk <- results_est_par$bsr_allele_rhat_sigma_all_mk
#        bsr_allele_trace <- results_est_par$bsr_allele_trace
#        bsr_lus_median_slope_all_mk <- results_est_par$bsr_lus_median_slope_all_mk
#        bsr_lus_median_intercept_all_mk <- results_est_par$bsr_lus_median_intercept_all_mk
#        bsr_lus_median_sigma_all_mk <- results_est_par$bsr_lus_median_sigma_all_mk
#        bsr_lus_rhat_slope_all_mk <- results_est_par$bsr_lus_rhat_slope_all_mk
#        bsr_lus_rhat_intercept_all_mk <- results_est_par$bsr_lus_rhat_intercept_all_mk
#        bsr_lus_rhat_sigma_all_mk <- results_est_par$bsr_lus_rhat_sigma_all_mk
#        bsr_lus_trace <- results_est_par$bsr_lus_trace
#        bsr_multi_median_slope_all_mk <- results_est_par$bsr_multi_median_slope_all_mk
#        bsr_multi_median_intercept_all_mk <- results_est_par$bsr_multi_median_intercept_all_mk
#        bsr_multi_median_sigma_all_mk <- results_est_par$bsr_multi_median_sigma_all_mk
#        bsr_multi_min_us_all_mk <- results_est_par$bsr_multi_min_us_all_mk
#        bsr_multi_rhat_slope_all_mk <- results_est_par$bsr_multi_rhat_slope_all_mk
#        bsr_multi_rhat_intercept_all_mk <- results_est_par$bsr_multi_rhat_intercept_all_mk
#        bsr_multi_rhat_sigma_all_mk <- results_est_par$bsr_multi_rhat_sigma_all_mk
#        bsr_multi_trace <- results_est_par$bsr_multi_trace
        bsr_best_model_all_mk <- results_est_par$bsr_best_model_all_mk
        bsr_best_median_myu_all_mk <- results_est_par$bsr_best_median_myu_all_mk
        bsr_best_median_slope_all_mk <- results_est_par$bsr_best_median_slope_all_mk
        bsr_best_median_intercept_all_mk <- results_est_par$bsr_best_median_intercept_all_mk
        bsr_best_median_sigma_all_mk <- results_est_par$bsr_best_median_sigma_all_mk
        bsr_best_min_us_all_mk <- results_est_par$bsr_best_min_us_all_mk
        bsr_best_rhat_myu_all_mk <- results_est_par$bsr_best_rhat_myu_all_mk
        bsr_best_rhat_slope_all_mk <- results_est_par$bsr_best_rhat_slope_all_mk
        bsr_best_rhat_intercept_all_mk <- results_est_par$bsr_best_rhat_intercept_all_mk
        bsr_best_rhat_sigma_all_mk <- results_est_par$bsr_best_rhat_sigma_all_mk
        bsr_best_trace_all_mk <- results_est_par$bsr_best_trace_all_mk
        rv_results$dt_result_bsr <- data.table(Marker = markers_par,
                                               Model = bsr_best_model_all_mk,
                                               Myu = signif(bsr_best_median_myu_all_mk, 3),
                                               Slope = signif(bsr_best_median_slope_all_mk, 3),
                                               Intercept = signif(bsr_best_median_intercept_all_mk, 3),
                                               Sigma = signif(bsr_best_median_sigma_all_mk, 3),
                                               Rhat_Myu = signif(bsr_best_rhat_myu_all_mk, 3),
                                               Rhat_Slope = signif(bsr_best_rhat_slope_all_mk, 3),
                                               Rhat_Intercept = signif(bsr_best_rhat_intercept_all_mk, 3),
                                               Rhat_Sigma = signif(bsr_best_rhat_sigma_all_mk, 3)
                                    )

        fsr_data_all_mk <- results_est_par$fsr_data_all_mk
#        fsr_allele_median_slope_all_mk <- results_est_par$fsr_allele_median_slope_all_mk
#        fsr_allele_median_intercept_all_mk <- results_est_par$fsr_allele_median_intercept_all_mk
#        fsr_allele_median_sigma_all_mk <- results_est_par$fsr_allele_median_sigma_all_mk
#        fsr_allele_rhat_slope_all_mk <- results_est_par$fsr_allele_rhat_slope_all_mk
#        fsr_allele_rhat_intercept_all_mk <- results_est_par$fsr_allele_rhat_intercept_all_mk
#        fsr_allele_rhat_sigma_all_mk <- results_est_par$fsr_allele_rhat_sigma_all_mk
#        fsr_allele_trace <- results_est_par$fsr_allele_trace
#        fsr_lus_median_slope_all_mk <- results_est_par$fsr_lus_median_slope_all_mk
#        fsr_lus_median_intercept_all_mk <- results_est_par$fsr_lus_median_intercept_all_mk
#        fsr_lus_median_sigma_all_mk <- results_est_par$fsr_lus_median_sigma_all_mk
#        fsr_lus_rhat_slope_all_mk <- results_est_par$fsr_lus_rhat_slope_all_mk
#        fsr_lus_rhat_intercept_all_mk <- results_est_par$fsr_lus_rhat_intercept_all_mk
#        fsr_lus_rhat_sigma_all_mk <- results_est_par$fsr_lus_rhat_sigma_all_mk
#        fsr_lus_trace <- results_est_par$fsr_lus_trace
#        fsr_multi_median_slope_all_mk <- results_est_par$fsr_multi_median_slope_all_mk
#        fsr_multi_median_intercept_all_mk <- results_est_par$fsr_multi_median_intercept_all_mk
#        fsr_multi_median_sigma_all_mk <- results_est_par$fsr_multi_median_sigma_all_mk
#        fsr_multi_min_us_all_mk <- results_est_par$fsr_multi_min_us_all_mk
#        fsr_multi_rhat_slope_all_mk <- results_est_par$fsr_multi_rhat_slope_all_mk
#        fsr_multi_rhat_intercept_all_mk <- results_est_par$fsr_multi_rhat_intercept_all_mk
#        fsr_multi_rhat_sigma_all_mk <- results_est_par$fsr_multi_rhat_sigma_all_mk
#        fsr_multi_trace <- results_est_par$fsr_multi_trace
        fsr_best_model_all_mk <- results_est_par$fsr_best_model_all_mk
        fsr_best_median_myu_all_mk <- results_est_par$fsr_best_median_myu_all_mk
        fsr_best_median_slope_all_mk <- results_est_par$fsr_best_median_slope_all_mk
        fsr_best_median_intercept_all_mk <- results_est_par$fsr_best_median_intercept_all_mk
        fsr_best_median_sigma_all_mk <- results_est_par$fsr_best_median_sigma_all_mk
        fsr_best_min_us_all_mk <- results_est_par$fsr_best_min_us_all_mk
        fsr_best_rhat_myu_all_mk <- results_est_par$fsr_best_rhat_myu_all_mk
        fsr_best_rhat_slope_all_mk <- results_est_par$fsr_best_rhat_slope_all_mk
        fsr_best_rhat_intercept_all_mk <- results_est_par$fsr_best_rhat_intercept_all_mk
        fsr_best_rhat_sigma_all_mk <- results_est_par$fsr_best_rhat_sigma_all_mk
        fsr_best_trace_all_mk <- results_est_par$fsr_best_trace_all_mk
        rv_results$dt_result_fsr <- data.table(Marker = markers_par,
                                               Model = fsr_best_model_all_mk,
                                               Myu = signif(fsr_best_median_myu_all_mk, 3),
                                               Slope = signif(fsr_best_median_slope_all_mk, 3),
                                               Intercept = signif(fsr_best_median_intercept_all_mk, 3),
                                               Sigma = signif(fsr_best_median_sigma_all_mk, 3),
                                               Rhat_Myu = signif(fsr_best_rhat_myu_all_mk, 3),
                                               Rhat_Slope = signif(fsr_best_rhat_slope_all_mk, 3),
                                               Rhat_Intercept = signif(fsr_best_rhat_intercept_all_mk, 3),
                                               Rhat_Sigma = signif(fsr_best_rhat_sigma_all_mk, 3)
                                    )

        dsr_data_all_mk <- results_est_par$dsr_data_all_mk
#        dsr_allele_median_slope_all_mk <- results_est_par$dsr_allele_median_slope_all_mk
#        dsr_allele_median_intercept_all_mk <- results_est_par$dsr_allele_median_intercept_all_mk
#        dsr_allele_median_sigma_all_mk <- results_est_par$dsr_allele_median_sigma_all_mk
#        dsr_allele_rhat_slope_all_mk <- results_est_par$dsr_allele_rhat_slope_all_mk
#        dsr_allele_rhat_intercept_all_mk <- results_est_par$dsr_allele_rhat_intercept_all_mk
#        dsr_allele_rhat_sigma_all_mk <- results_est_par$dsr_allele_rhat_sigma_all_mk
#        dsr_allele_trace <- results_est_par$dsr_allele_trace
#        dsr_lus_median_slope_all_mk <- results_est_par$dsr_lus_median_slope_all_mk
#        dsr_lus_median_intercept_all_mk <- results_est_par$dsr_lus_median_intercept_all_mk
#        dsr_lus_median_sigma_all_mk <- results_est_par$dsr_lus_median_sigma_all_mk
#        dsr_lus_rhat_slope_all_mk <- results_est_par$dsr_lus_rhat_slope_all_mk
#        dsr_lus_rhat_intercept_all_mk <- results_est_par$dsr_lus_rhat_intercept_all_mk
#        dsr_lus_rhat_sigma_all_mk <- results_est_par$dsr_lus_rhat_sigma_all_mk
#        dsr_lus_trace <- results_est_par$dsr_lus_trace
#        dsr_multi_median_slope_all_mk <- results_est_par$dsr_multi_median_slope_all_mk
#        dsr_multi_median_intercept_all_mk <- results_est_par$dsr_multi_median_intercept_all_mk
#        dsr_multi_median_sigma_all_mk <- results_est_par$dsr_multi_median_sigma_all_mk
#        dsr_multi_min_us_all_mk <- results_est_par$dsr_multi_min_us_all_mk
#        dsr_multi_rhat_slope_all_mk <- results_est_par$dsr_multi_rhat_slope_all_mk
#        dsr_multi_rhat_intercept_all_mk <- results_est_par$dsr_multi_rhat_intercept_all_mk
#        dsr_multi_rhat_sigma_all_mk <- results_est_par$dsr_multi_rhat_sigma_all_mk
#        dsr_multi_trace <- results_est_par$dsr_multi_trace
        dsr_best_model_all_mk <- results_est_par$dsr_best_model_all_mk
        dsr_best_median_myu_all_mk <- results_est_par$dsr_best_median_myu_all_mk
        dsr_best_median_slope_all_mk <- results_est_par$dsr_best_median_slope_all_mk
        dsr_best_median_intercept_all_mk <- results_est_par$dsr_best_median_intercept_all_mk
        dsr_best_median_sigma_all_mk <- results_est_par$dsr_best_median_sigma_all_mk
        dsr_best_min_us_all_mk <- results_est_par$dsr_best_min_us_all_mk
        dsr_best_rhat_myu_all_mk <- results_est_par$dsr_best_rhat_myu_all_mk
        dsr_best_rhat_slope_all_mk <- results_est_par$dsr_best_rhat_slope_all_mk
        dsr_best_rhat_intercept_all_mk <- results_est_par$dsr_best_rhat_intercept_all_mk
        dsr_best_rhat_sigma_all_mk <- results_est_par$dsr_best_rhat_sigma_all_mk
        dsr_best_trace_all_mk <- results_est_par$dsr_best_trace_all_mk
        rv_results$dt_result_dsr <- data.table(Marker = markers_par,Model = dsr_best_model_all_mk,
                                               Myu = signif(dsr_best_median_myu_all_mk, 3),
                                               Slope = signif(dsr_best_median_slope_all_mk, 3),
                                               Intercept = signif(dsr_best_median_intercept_all_mk, 3),
                                               Sigma = signif(dsr_best_median_sigma_all_mk, 3),
                                               Rhat_Myu = signif(dsr_best_rhat_myu_all_mk, 3),
                                               Rhat_Slope = signif(dsr_best_rhat_slope_all_mk, 3),
                                               Rhat_Intercept = signif(dsr_best_rhat_intercept_all_mk, 3),
                                               Rhat_Sigma = signif(dsr_best_rhat_sigma_all_mk, 3)
                                    )

        m2sr_data_all_mk <- results_est_par$m2sr_data_all_mk
        m2sr_best_model_all_mk <- results_est_par$m2sr_best_model_all_mk
        m2sr_best_median_myu_all_mk <- results_est_par$m2sr_best_median_myu_all_mk
        m2sr_best_median_sigma_all_mk <- results_est_par$m2sr_best_median_sigma_all_mk
        m2sr_best_rhat_myu_all_mk <- results_est_par$m2sr_best_rhat_myu_all_mk
        m2sr_best_rhat_sigma_all_mk <- results_est_par$m2sr_best_rhat_sigma_all_mk
        m2sr_best_trace_all_mk <- results_est_par$m2sr_best_trace_all_mk
        rv_results$dt_result_m2sr <- data.table(Marker = markers_par,
                                                Model = m2sr_best_model_all_mk,
                                                Myu = signif(m2sr_best_median_myu_all_mk, 3),
                                                Sigma = signif(m2sr_best_median_sigma_all_mk, 3),
                                                Rhat_Myu = signif(m2sr_best_rhat_myu_all_mk, 3),
                                                Rhat_Sigma = signif(m2sr_best_rhat_sigma_all_mk, 3)
                                     )
      }

      rv_graph <- reactiveValues()
      rv_graph$trace_ae <- NULL
      rv_graph$trace_bsr <- NULL
      rv_graph$xval_bsr <- NULL
      rv_graph$yval_bsr <- NULL
      rv_graph$trace_fsr <- NULL
      rv_graph$xval_fsr <- NULL
      rv_graph$yval_fsr <- NULL
      rv_graph$trace_dsr <- NULL
      rv_graph$xval_dsr <- NULL
      rv_graph$yval_dsr <- NULL
      rv_graph$trace_m2sr <- NULL
      rv_graph$xval_m2sr <- NULL
      rv_graph$yval_m2sr <- NULL

      observeEvent(input$act_trace_ae, {
        if(!is.null(ae_trace_all_mk)){
          row_index_select <- input$dt_result_ae_rows_selected
          if(isTruthy(row_index_select)){
            rv_graph$trace_ae <- ae_trace_all_mk[[row_index_select]]

            showModal(modalDialog(
              title = "Trace plot",
              plotOutput(session$ns("trace_ae")),
              footer = tagList(
                modalButton("Close")
              ),
              size = "l"
            ))

          }else{
            showModal(modalDialog(title = "Error", "Select a row!", easyClose = TRUE, footer = NULL))
          }
        }
      })

      output$trace_ae <- renderPlot({
        rv_graph$trace_ae
      })

      observeEvent(input$act_plot_ae, {
        if(!is.null(ae_data)){
          showModal(modalDialog(
            title = "Data plot",
            plotOutput(session$ns("dataplot_ae")),
            footer = tagList(
              modalButton("Close")
            ),
            size = "l"
          ))
        }
      })

      output$dataplot_ae <- renderPlot({
        boxplot(ae_data)
      })

      observeEvent(input$act_trace_bsr, {
        if(!is.null(bsr_best_trace_all_mk)){
          row_index_select <- input$dt_result_bsr_rows_selected
          if(isTruthy(row_index_select)){
            isolate(dt_result_bsr <- rv_results$dt_result_bsr)
            if(dt_result_bsr$Model[row_index_select] != "Not consider"){
              rv_graph$trace_bsr <- bsr_best_trace_all_mk[[row_index_select]]

              showModal(modalDialog(
                title = "Trace plot",
                plotOutput(session$ns("trace_bsr")),
                footer = tagList(
                  modalButton("Close")
                ),
                size = "l"
              ))
            }else{
              showModal(modalDialog(title = "Information", "The back stutter ratio is not considered in the selected marker!", easyClose = TRUE, footer = NULL))
            }
          }else{
            showModal(modalDialog(title = "Error", "Select a row!", easyClose = TRUE, footer = NULL))
          }
        }
      })

      output$trace_bsr <- renderPlot({
        rv_graph$trace_bsr
      })

      observeEvent(input$act_plot_bsr, {
        if(!is.null(bsr_data_all_mk)){
          row_index_select <- input$dt_result_bsr_rows_selected
          if(isTruthy(row_index_select)){
            isolate(dt_result_bsr <- rv_results$dt_result_bsr)
            model_bsr <- dt_result_bsr$Model[row_index_select]
            if(model_bsr != "Not consider"){
              bsr_data_selected <- bsr_data_all_mk[[row_index_select]]
              if(model_bsr == "LUS"){
                rv_graph$xval_bsr <- bsr_data_selected$bsr_lus
              }else if(model_bsr == "Multi-sequence"){
                rv_graph$xval_bsr <- bsr_data_selected$bsr_ca
              }else{
                rv_graph$xval_bsr <- bsr_data_selected$bsr_al
              }
              rv_graph$yval_bsr <- bsr_data_selected$bsr

              showModal(modalDialog(
                title = "Data plot",
                plotOutput(session$ns("dataplot_bsr")),
                footer = tagList(
                  modalButton("Close")
                ),
                size = "l"
              ))
            }else{
              showModal(modalDialog(title = "Information", "The back stutter ratio is not considered in the selected marker!", easyClose = TRUE, footer = NULL))
            }
          }else{
            showModal(modalDialog(title = "Error", "Select a row!", easyClose = TRUE, footer = NULL))
          }
        }
      })

      output$dataplot_bsr <- renderPlot({
        plot(rv_graph$xval_bsr, rv_graph$yval_bsr)
      })

      observeEvent(input$act_trace_fsr, {
        if(!is.null(fsr_best_trace_all_mk)){
          row_index_select <- input$dt_result_fsr_rows_selected
          if(isTruthy(row_index_select)){
            isolate(dt_result_fsr <- rv_results$dt_result_fsr)
            if(dt_result_fsr$Model[row_index_select] != "Not consider"){
              rv_graph$trace_fsr <- fsr_best_trace_all_mk[[row_index_select]]

              showModal(modalDialog(
                title = "Trace plot",
                plotOutput(session$ns("trace_fsr")),
                footer = tagList(
                  modalButton("Close")
                ),
                size = "l"
              ))
            }else{
              showModal(modalDialog(title = "Information", "The forward stutter ratio is not considered in the selected marker!", easyClose = TRUE, footer = NULL))
            }
          }else{
            showModal(modalDialog(title = "Error", "Select a row!", easyClose = TRUE, footer = NULL))
          }
        }
      })

      output$trace_fsr <- renderPlot({
        rv_graph$trace_fsr
      })

      observeEvent(input$act_trace_dsr, {
        if(!is.null(dsr_best_trace_all_mk)){
          row_index_select <- input$dt_result_dsr_rows_selected
          if(isTruthy(row_index_select)){
            isolate(dt_result_dsr <- rv_results$dt_result_dsr)
            if(dt_result_dsr$Model[row_index_select] != "Not consider"){
              rv_graph$trace_dsr <- dsr_best_trace_all_mk[[row_index_select]]

              showModal(modalDialog(
                title = "Trace plot",
                plotOutput(session$ns("trace_dsr")),
                footer = tagList(
                  modalButton("Close")
                ),
                size = "l"
              ))
            }else{
              showModal(modalDialog(title = "Information", "The double-back stutter ratio is not considered in the selected marker!", easyClose = TRUE, footer = NULL))
            }
          }else{
            showModal(modalDialog(title = "Error", "Select a row!", easyClose = TRUE, footer = NULL))
          }
        }
      })

      output$trace_dsr <- renderPlot({
        rv_graph$trace_dsr
      })

      observeEvent(input$act_trace_m2sr, {
        if(!is.null(m2sr_best_trace_all_mk)){
          row_index_select <- input$dt_result_m2sr_rows_selected
          if(isTruthy(row_index_select)){
            isolate(dt_result_m2sr <- rv_results$dt_result_m2sr)
            if(dt_result_m2sr$Model[row_index_select] != "Not consider"){
              rv_graph$trace_m2sr <- m2sr_best_trace_all_mk[[row_index_select]]

              showModal(modalDialog(
                title = "Trace plot",
                plotOutput(session$ns("trace_m2sr")),
                footer = tagList(
                  modalButton("Close")
                ),
                size = "l"
              ))
            }else{
              showModal(modalDialog(title = "Information", "The minus 2-nt stutter ratio is not considered in the selected marker!", easyClose = TRUE, footer = NULL))
            }
          }else{
            showModal(modalDialog(title = "Error", "Select a row!", easyClose = TRUE, footer = NULL))
          }
        }
      })

      output$trace_m2sr <- renderPlot({
        rv_graph$trace_m2sr
      })

      output$dt_result_ae <- renderDataTable(server = FALSE, {
        datatable(
          rv_results$dt_result_ae,
          selection = list(mode = "single", target = "row"),
          options = list(iDisplayLength = 50, ordering = FALSE, autowidth = TRUE),
          rownames = FALSE
        )
      })

      output$dt_result_bsr <- renderDataTable(server = FALSE, {
        datatable(
          rv_results$dt_result_bsr,
          selection = list(mode = "single", target = "row"),
          options = list(iDisplayLength = 50, ordering = FALSE, autowidth = TRUE),
          rownames = FALSE
        )
      })

      output$dt_result_fsr <- renderDataTable(server = FALSE, {
        datatable(
          rv_results$dt_result_fsr,
          selection = list(mode = "single", target = "row"),
          options = list(iDisplayLength = 50, ordering = FALSE, autowidth = TRUE),
          rownames = FALSE
        )
      })

      output$dt_result_dsr <- renderDataTable(server = FALSE, {
        datatable(
          rv_results$dt_result_dsr,
          selection = list(mode = "single", target = "row"),
          options = list(iDisplayLength = 50, ordering = FALSE, autowidth = TRUE),
          rownames = FALSE
        )
      })

      output$dt_result_m2sr <- renderDataTable(server = FALSE, {
        datatable(
          rv_results$dt_result_m2sr,
          selection = list(mode = "single", target = "row"),
          options = list(iDisplayLength = 50, ordering = FALSE, autowidth = TRUE),
          rownames = FALSE
        )
      })
    }
  )
}
