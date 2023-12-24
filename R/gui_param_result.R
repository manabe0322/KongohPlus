est_par_result_ui <- function(id){
  ns <- NS(id)

  tabPanel("Result",
           br(),
           sidebarLayout(
             sidebarPanel(width = 2,
                          actionButton(ns("act_change_model"), label = "Change model"),
                          br(),
                          br(),
                          br(),
                          actionButton(ns("act_traceplot"), label = "Trace plot"),
                          br(),
                          br(),
                          br(),
                          actionButton(ns("act_dataplot"), label = "Data plot")
             ),
             mainPanel(width = 10,
                       tabsetPanel(id = ns("tabset_est_par_result"),
                                   type = "tabs",
                                   tabPanel("AE",
                                            br(),
                                            dataTableOutput(ns("dt_result_ae"))
                                   ),
                                   tabPanel("BSR",
                                            br(),
                                            dataTableOutput(ns("dt_result_bsr"))
                                   ),
                                   tabPanel("FSR",
                                            br(),
                                            dataTableOutput(ns("dt_result_fsr"))
                                   ),
                                   tabPanel("DSR",
                                            br(),
                                            dataTableOutput(ns("dt_result_dsr"))
                                   ),
                                   tabPanel("M2SR",
                                            br(),
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
      rv_results$traceplot_bsr <- NULL
      rv_results$traceplot_fsr <- NULL
      rv_results$traceplot_dsr <- NULL
      if(is.null(results_est_par)){
        markers_par <- NULL

        ae_data <- NULL
        ae_median_myu_all_mk <- NULL
        ae_median_sigma_all_mk <- NULL
        ae_rhat_myu_all_mk <- NULL
        ae_rhat_sigma_all_mk <- NULL
        ae_rhat_lp_all_mk <- NULL
        ae_trace_all_mk <- NULL

        bsr_data_all_mk <- NULL
        aics_bsr_all_mk <- NULL
        bsr_allele_median_slope_all_mk <- NULL
        bsr_allele_median_intercept_all_mk <- NULL
        bsr_allele_median_sigma_all_mk <- NULL
        bsr_allele_rhat_slope_all_mk <- NULL
        bsr_allele_rhat_intercept_all_mk <- NULL
        bsr_allele_rhat_sigma_all_mk <- NULL
        bsr_allele_rhat_lp_all_mk <- NULL
        bsr_allele_trace_all_mk <- NULL
        bsr_lus_median_slope_all_mk <- NULL
        bsr_lus_median_intercept_all_mk <- NULL
        bsr_lus_median_sigma_all_mk <- NULL
        bsr_lus_rhat_slope_all_mk <- NULL
        bsr_lus_rhat_intercept_all_mk <- NULL
        bsr_lus_rhat_sigma_all_mk <- NULL
        bsr_lus_rhat_lp_all_mk <- NULL
        bsr_lus_trace_all_mk <- NULL
        bsr_multi_median_slope_all_mk <- NULL
        bsr_multi_median_intercept_all_mk <- NULL
        bsr_multi_median_sigma_all_mk <- NULL
        bsr_multi_min_us_all_mk <- NULL
        bsr_multi_rhat_slope_all_mk <- NULL
        bsr_multi_rhat_intercept_all_mk <- NULL
        bsr_multi_rhat_sigma_all_mk <- NULL
        bsr_multi_rhat_lp_all_mk <- NULL
        bsr_multi_trace_all_mk <- NULL
        bsr_best_model_all_mk <- NULL
        bsr_best_median_myu_all_mk <- NULL
        bsr_best_median_slope_all_mk <- NULL
        bsr_best_median_intercept_all_mk <- NULL
        bsr_best_median_sigma_all_mk <- NULL
        bsr_best_min_us_all_mk <- NULL
        bsr_best_rhat_myu_all_mk <- NULL
        bsr_best_rhat_slope_all_mk <- NULL
        bsr_best_rhat_intercept_all_mk <- NULL
        bsr_best_rhat_sigma_all_mk <- NULL
        bsr_best_rhat_lp_all_mk <- NULL
        bsr_best_trace_all_mk <- NULL

        fsr_data_all_mk <- NULL
        aics_fsr_all_mk <- NULL
        fsr_allele_median_slope_all_mk <- NULL
        fsr_allele_median_intercept_all_mk <- NULL
        fsr_allele_median_sigma_all_mk <- NULL
        fsr_allele_rhat_slope_all_mk <- NULL
        fsr_allele_rhat_intercept_all_mk <- NULL
        fsr_allele_rhat_sigma_all_mk <- NULL
        fsr_allele_rhat_lp_all_mk <- NULL
        fsr_allele_trace_all_mk <- NULL
        fsr_lus_median_slope_all_mk <- NULL
        fsr_lus_median_intercept_all_mk <- NULL
        fsr_lus_median_sigma_all_mk <- NULL
        fsr_lus_rhat_slope_all_mk <- NULL
        fsr_lus_rhat_intercept_all_mk <- NULL
        fsr_lus_rhat_sigma_all_mk <- NULL
        fsr_lus_rhat_lp_all_mk <- NULL
        fsr_lus_trace_all_mk <- NULL
        fsr_multi_median_slope_all_mk <- NULL
        fsr_multi_median_intercept_all_mk <- NULL
        fsr_multi_median_sigma_all_mk <- NULL
        fsr_multi_min_us_all_mk <- NULL
        fsr_multi_rhat_slope_all_mk <- NULL
        fsr_multi_rhat_intercept_all_mk <- NULL
        fsr_multi_rhat_sigma_all_mk <- NULL
        fsr_multi_rhat_lp_all_mk <- NULL
        fsr_multi_trace_all_mk <- NULL
        fsr_best_model_all_mk <- NULL
        fsr_best_median_myu_all_mk <- NULL
        fsr_best_median_slope_all_mk <- NULL
        fsr_best_median_intercept_all_mk <- NULL
        fsr_best_median_sigma_all_mk <- NULL
        fsr_best_min_us_all_mk <- NULL
        fsr_best_rhat_myu_all_mk <- NULL
        fsr_best_rhat_slope_all_mk <- NULL
        fsr_best_rhat_intercept_all_mk <- NULL
        fsr_best_rhat_sigma_all_mk <- NULL
        fsr_best_rhat_lp_all_mk <- NULL
        fsr_best_trace_all_mk <- NULL

        dsr_data_all_mk <- NULL
        aics_dsr_all_mk <- NULL
        dsr_allele_median_slope_all_mk <- NULL
        dsr_allele_median_intercept_all_mk <- NULL
        dsr_allele_median_sigma_all_mk <- NULL
        dsr_allele_rhat_slope_all_mk <- NULL
        dsr_allele_rhat_intercept_all_mk <- NULL
        dsr_allele_rhat_sigma_all_mk <- NULL
        dsr_allele_rhat_lp_all_mk <- NULL
        dsr_allele_trace_all_mk <- NULL
        dsr_lus_median_slope_all_mk <- NULL
        dsr_lus_median_intercept_all_mk <- NULL
        dsr_lus_median_sigma_all_mk <- NULL
        dsr_lus_rhat_slope_all_mk <- NULL
        dsr_lus_rhat_intercept_all_mk <- NULL
        dsr_lus_rhat_sigma_all_mk <- NULL
        dsr_lus_rhat_lp_all_mk <- NULL
        dsr_lus_trace_all_mk <- NULL
        dsr_multi_median_slope_all_mk <- NULL
        dsr_multi_median_intercept_all_mk <- NULL
        dsr_multi_median_sigma_all_mk <- NULL
        dsr_multi_min_us_all_mk <- NULL
        dsr_multi_rhat_slope_all_mk <- NULL
        dsr_multi_rhat_intercept_all_mk <- NULL
        dsr_multi_rhat_sigma_all_mk <- NULL
        dsr_multi_rhat_lp_all_mk <- NULL
        dsr_multi_trace_all_mk <- NULL
        dsr_best_model_all_mk <- NULL
        dsr_best_median_myu_all_mk <- NULL
        dsr_best_median_slope_all_mk <- NULL
        dsr_best_median_intercept_all_mk <- NULL
        dsr_best_median_sigma_all_mk <- NULL
        dsr_best_min_us_all_mk <- NULL
        dsr_best_rhat_myu_all_mk <- NULL
        dsr_best_rhat_slope_all_mk <- NULL
        dsr_best_rhat_intercept_all_mk <- NULL
        dsr_best_rhat_sigma_all_mk <- NULL
        dsr_best_rhat_lp_all_mk <- NULL
        dsr_best_trace_all_mk <- NULL

        m2sr_data_all_mk <- NULL
        m2sr_best_model_all_mk <- NULL
        m2sr_best_median_myu_all_mk <- NULL
        m2sr_best_median_sigma_all_mk <- NULL
        m2sr_best_rhat_myu_all_mk <- NULL
        m2sr_best_rhat_sigma_all_mk <- NULL
        m2sr_best_rhat_lp_all_mk <- NULL
        m2sr_best_trace_all_mk <- NULL
      }else{
        markers_par <- results_est_par$markers_par

        ae_data <- results_est_par$ae_data
        ae_median_myu_all_mk <- results_est_par$ae_median_myu_all_mk
        ae_median_sigma_all_mk <- results_est_par$ae_median_sigma_all_mk
        ae_rhat_myu_all_mk <- results_est_par$ae_rhat_myu_all_mk
        ae_rhat_sigma_all_mk <- results_est_par$ae_rhat_sigma_all_mk
        ae_rhat_lp_all_mk <- results_est_par$ae_rhat_lp_all_mk
        ae_trace_all_mk <- results_est_par$ae_trace_all_mk
        rv_results$dt_result_ae <- data.table(Marker = markers_par,
                                              Myu = signif(ae_median_myu_all_mk, 3),
                                              Sigma = signif(ae_median_sigma_all_mk, 3),
                                              Rhat_Myu = signif(ae_rhat_myu_all_mk, 3),
                                              Rhat_Sigma = signif(ae_rhat_sigma_all_mk, 3),
                                              Rhat_Likelihood = signif(ae_rhat_lp_all_mk, 3)
                                              )

        bsr_data_all_mk <- results_est_par$bsr_data_all_mk
        aics_bsr_all_mk <- results_est_par$aics_bsr_all_mk
        bsr_allele_median_slope_all_mk <- results_est_par$bsr_allele_median_slope_all_mk
        bsr_allele_median_intercept_all_mk <- results_est_par$bsr_allele_median_intercept_all_mk
        bsr_allele_median_sigma_all_mk <- results_est_par$bsr_allele_median_sigma_all_mk
        bsr_allele_rhat_slope_all_mk <- results_est_par$bsr_allele_rhat_slope_all_mk
        bsr_allele_rhat_intercept_all_mk <- results_est_par$bsr_allele_rhat_intercept_all_mk
        bsr_allele_rhat_sigma_all_mk <- results_est_par$bsr_allele_rhat_sigma_all_mk
        bsr_allele_rhat_lp_all_mk <- results_est_par$bsr_allele_rhat_lp_all_mk
        bsr_allele_trace_all_mk <- results_est_par$bsr_allele_trace_all_mk
        bsr_lus_median_slope_all_mk <- results_est_par$bsr_lus_median_slope_all_mk
        bsr_lus_median_intercept_all_mk <- results_est_par$bsr_lus_median_intercept_all_mk
        bsr_lus_median_sigma_all_mk <- results_est_par$bsr_lus_median_sigma_all_mk
        bsr_lus_rhat_slope_all_mk <- results_est_par$bsr_lus_rhat_slope_all_mk
        bsr_lus_rhat_intercept_all_mk <- results_est_par$bsr_lus_rhat_intercept_all_mk
        bsr_lus_rhat_sigma_all_mk <- results_est_par$bsr_lus_rhat_sigma_all_mk
        bsr_lus_rhat_lp_all_mk <- results_est_par$bsr_lus_rhat_lp_all_mk
        bsr_lus_trace_all_mk <- results_est_par$bsr_lus_trace_all_mk
        bsr_multi_median_slope_all_mk <- results_est_par$bsr_multi_median_slope_all_mk
        bsr_multi_median_intercept_all_mk <- results_est_par$bsr_multi_median_intercept_all_mk
        bsr_multi_median_sigma_all_mk <- results_est_par$bsr_multi_median_sigma_all_mk
        bsr_multi_min_us_all_mk <- results_est_par$bsr_multi_min_us_all_mk
        bsr_multi_rhat_slope_all_mk <- results_est_par$bsr_multi_rhat_slope_all_mk
        bsr_multi_rhat_intercept_all_mk <- results_est_par$bsr_multi_rhat_intercept_all_mk
        bsr_multi_rhat_sigma_all_mk <- results_est_par$bsr_multi_rhat_sigma_all_mk
        bsr_multi_rhat_lp_all_mk <- results_est_par$bsr_multi_rhat_lp_all_mk
        bsr_multi_trace_all_mk <- results_est_par$bsr_multi_trace_all_mk
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
        bsr_best_rhat_lp_all_mk <- results_est_par$bsr_best_rhat_lp_all_mk
        bsr_best_trace_all_mk <- results_est_par$bsr_best_trace_all_mk
        rv_results$dt_result_bsr <- data.table(Marker = markers_par,
                                               Model = bsr_best_model_all_mk,
                                               Myu = signif(bsr_best_median_myu_all_mk, 3),
                                               Slope = signif(bsr_best_median_slope_all_mk, 3),
                                               Intercept = signif(bsr_best_median_intercept_all_mk, 3),
                                               Sigma = signif(bsr_best_median_sigma_all_mk, 3),
                                               MinUS = bsr_best_min_us_all_mk,
                                               Rhat_Myu = signif(bsr_best_rhat_myu_all_mk, 3),
                                               Rhat_Slope = signif(bsr_best_rhat_slope_all_mk, 3),
                                               Rhat_Intercept = signif(bsr_best_rhat_intercept_all_mk, 3),
                                               Rhat_Sigma = signif(bsr_best_rhat_sigma_all_mk, 3),
                                               Rhat_Likelihood = signif(bsr_best_rhat_lp_all_mk, 3)
                                               )
        rv_results$traceplot_bsr <- bsr_best_trace_all_mk

        fsr_data_all_mk <- results_est_par$fsr_data_all_mk
        aics_fsr_all_mk <- results_est_par$aics_fsr_all_mk
        fsr_allele_median_slope_all_mk <- results_est_par$fsr_allele_median_slope_all_mk
        fsr_allele_median_intercept_all_mk <- results_est_par$fsr_allele_median_intercept_all_mk
        fsr_allele_median_sigma_all_mk <- results_est_par$fsr_allele_median_sigma_all_mk
        fsr_allele_rhat_slope_all_mk <- results_est_par$fsr_allele_rhat_slope_all_mk
        fsr_allele_rhat_intercept_all_mk <- results_est_par$fsr_allele_rhat_intercept_all_mk
        fsr_allele_rhat_sigma_all_mk <- results_est_par$fsr_allele_rhat_sigma_all_mk
        fsr_allele_rhat_lp_all_mk <- results_est_par$fsr_allele_rhat_lp_all_mk
        fsr_allele_trace_all_mk <- results_est_par$fsr_allele_trace_all_mk
        fsr_lus_median_slope_all_mk <- results_est_par$fsr_lus_median_slope_all_mk
        fsr_lus_median_intercept_all_mk <- results_est_par$fsr_lus_median_intercept_all_mk
        fsr_lus_median_sigma_all_mk <- results_est_par$fsr_lus_median_sigma_all_mk
        fsr_lus_rhat_slope_all_mk <- results_est_par$fsr_lus_rhat_slope_all_mk
        fsr_lus_rhat_intercept_all_mk <- results_est_par$fsr_lus_rhat_intercept_all_mk
        fsr_lus_rhat_sigma_all_mk <- results_est_par$fsr_lus_rhat_sigma_all_mk
        fsr_lus_rhat_lp_all_mk <- results_est_par$fsr_lus_rhat_lp_all_mk
        fsr_lus_trace_all_mk <- results_est_par$fsr_lus_trace_all_mk
        fsr_multi_median_slope_all_mk <- results_est_par$fsr_multi_median_slope_all_mk
        fsr_multi_median_intercept_all_mk <- results_est_par$fsr_multi_median_intercept_all_mk
        fsr_multi_median_sigma_all_mk <- results_est_par$fsr_multi_median_sigma_all_mk
        fsr_multi_min_us_all_mk <- results_est_par$fsr_multi_min_us_all_mk
        fsr_multi_rhat_slope_all_mk <- results_est_par$fsr_multi_rhat_slope_all_mk
        fsr_multi_rhat_intercept_all_mk <- results_est_par$fsr_multi_rhat_intercept_all_mk
        fsr_multi_rhat_sigma_all_mk <- results_est_par$fsr_multi_rhat_sigma_all_mk
        fsr_multi_rhat_lp_all_mk <- results_est_par$fsr_multi_rhat_lp_all_mk
        fsr_multi_trace_all_mk <- results_est_par$fsr_multi_trace_all_mk
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
        fsr_best_rhat_lp_all_mk <- results_est_par$fsr_best_rhat_lp_all_mk
        fsr_best_trace_all_mk <- results_est_par$fsr_best_trace_all_mk
        rv_results$dt_result_fsr <- data.table(Marker = markers_par,
                                               Model = fsr_best_model_all_mk,
                                               Myu = signif(fsr_best_median_myu_all_mk, 3),
                                               Slope = signif(fsr_best_median_slope_all_mk, 3),
                                               Intercept = signif(fsr_best_median_intercept_all_mk, 3),
                                               Sigma = signif(fsr_best_median_sigma_all_mk, 3),
                                               MinUS = fsr_best_min_us_all_mk,
                                               Rhat_Myu = signif(fsr_best_rhat_myu_all_mk, 3),
                                               Rhat_Slope = signif(fsr_best_rhat_slope_all_mk, 3),
                                               Rhat_Intercept = signif(fsr_best_rhat_intercept_all_mk, 3),
                                               Rhat_Sigma = signif(fsr_best_rhat_sigma_all_mk, 3),
                                               Rhat_Likelihood = signif(fsr_best_rhat_lp_all_mk, 3)
                                               )
        rv_results$traceplot_fsr <- fsr_best_trace_all_mk

        dsr_data_all_mk <- results_est_par$dsr_data_all_mk
        aics_dsr_all_mk <- results_est_par$aics_dsr_all_mk
        dsr_allele_median_slope_all_mk <- results_est_par$dsr_allele_median_slope_all_mk
        dsr_allele_median_intercept_all_mk <- results_est_par$dsr_allele_median_intercept_all_mk
        dsr_allele_median_sigma_all_mk <- results_est_par$dsr_allele_median_sigma_all_mk
        dsr_allele_rhat_slope_all_mk <- results_est_par$dsr_allele_rhat_slope_all_mk
        dsr_allele_rhat_intercept_all_mk <- results_est_par$dsr_allele_rhat_intercept_all_mk
        dsr_allele_rhat_sigma_all_mk <- results_est_par$dsr_allele_rhat_sigma_all_mk
        dsr_allele_rhat_lp_all_mk <- results_est_par$dsr_allele_rhat_lp_all_mk
        dsr_allele_trace_all_mk <- results_est_par$dsr_allele_trace_all_mk
        dsr_lus_median_slope_all_mk <- results_est_par$dsr_lus_median_slope_all_mk
        dsr_lus_median_intercept_all_mk <- results_est_par$dsr_lus_median_intercept_all_mk
        dsr_lus_median_sigma_all_mk <- results_est_par$dsr_lus_median_sigma_all_mk
        dsr_lus_rhat_slope_all_mk <- results_est_par$dsr_lus_rhat_slope_all_mk
        dsr_lus_rhat_intercept_all_mk <- results_est_par$dsr_lus_rhat_intercept_all_mk
        dsr_lus_rhat_sigma_all_mk <- results_est_par$dsr_lus_rhat_sigma_all_mk
        dsr_lus_rhat_lp_all_mk <- results_est_par$dsr_lus_rhat_lp_all_mk
        dsr_lus_trace_all_mk <- results_est_par$dsr_lus_trace_all_mk
        dsr_multi_median_slope_all_mk <- results_est_par$dsr_multi_median_slope_all_mk
        dsr_multi_median_intercept_all_mk <- results_est_par$dsr_multi_median_intercept_all_mk
        dsr_multi_median_sigma_all_mk <- results_est_par$dsr_multi_median_sigma_all_mk
        dsr_multi_min_us_all_mk <- results_est_par$dsr_multi_min_us_all_mk
        dsr_multi_rhat_slope_all_mk <- results_est_par$dsr_multi_rhat_slope_all_mk
        dsr_multi_rhat_intercept_all_mk <- results_est_par$dsr_multi_rhat_intercept_all_mk
        dsr_multi_rhat_sigma_all_mk <- results_est_par$dsr_multi_rhat_sigma_all_mk
        dsr_multi_rhat_lp_all_mk <- results_est_par$dsr_multi_rhat_lp_all_mk
        dsr_multi_trace_all_mk <- results_est_par$dsr_multi_trace_all_mk
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
        dsr_best_rhat_lp_all_mk <- results_est_par$dsr_best_rhat_lp_all_mk
        dsr_best_trace_all_mk <- results_est_par$dsr_best_trace_all_mk
        rv_results$dt_result_dsr <- data.table(Marker = markers_par,Model = dsr_best_model_all_mk,
                                               Myu = signif(dsr_best_median_myu_all_mk, 3),
                                               Slope = signif(dsr_best_median_slope_all_mk, 3),
                                               Intercept = signif(dsr_best_median_intercept_all_mk, 3),
                                               Sigma = signif(dsr_best_median_sigma_all_mk, 3),
                                               MinUS = dsr_best_min_us_all_mk,
                                               Rhat_Myu = signif(dsr_best_rhat_myu_all_mk, 3),
                                               Rhat_Slope = signif(dsr_best_rhat_slope_all_mk, 3),
                                               Rhat_Intercept = signif(dsr_best_rhat_intercept_all_mk, 3),
                                               Rhat_Sigma = signif(dsr_best_rhat_sigma_all_mk, 3),
                                               Rhat_Likelihood = signif(dsr_best_rhat_lp_all_mk, 3)
                                               )
        rv_results$traceplot_dsr <- dsr_best_trace_all_mk

        m2sr_data_all_mk <- results_est_par$m2sr_data_all_mk
        m2sr_best_model_all_mk <- results_est_par$m2sr_best_model_all_mk
        m2sr_best_median_myu_all_mk <- results_est_par$m2sr_best_median_myu_all_mk
        m2sr_best_median_sigma_all_mk <- results_est_par$m2sr_best_median_sigma_all_mk
        m2sr_best_rhat_myu_all_mk <- results_est_par$m2sr_best_rhat_myu_all_mk
        m2sr_best_rhat_sigma_all_mk <- results_est_par$m2sr_best_rhat_sigma_all_mk
        m2sr_best_rhat_lp_all_mk <- results_est_par$m2sr_best_rhat_lp_all_mk
        m2sr_best_trace_all_mk <- results_est_par$m2sr_best_trace_all_mk
        rv_results$dt_result_m2sr <- data.table(Marker = markers_par,
                                                Model = m2sr_best_model_all_mk,
                                                Myu = signif(m2sr_best_median_myu_all_mk, 3),
                                                Sigma = signif(m2sr_best_median_sigma_all_mk, 3),
                                                Rhat_Myu = signif(m2sr_best_rhat_myu_all_mk, 3),
                                                Rhat_Sigma = signif(m2sr_best_rhat_sigma_all_mk, 3),
                                                Rhat_Likelihood = signif(m2sr_best_rhat_lp_all_mk, 3)
                                                )
      }

      #### Change model ####

      # Define reactive values for displaying AIC values and parameters of each candidate model
      rv_dt_change_model <- reactiveValues()
      rv_dt_change_model$bsr <- NULL
      rv_dt_change_model$fsr <- NULL
      rv_dt_change_model$dsr <- NULL

      # Change the status of the "Change model" button
      observeEvent(input$tabset_est_par_result, {
        if(is.element(input$tabset_est_par_result, c("AE", "M2SR"))){
          disable("act_change_model")
        }else{
          enable("act_change_model")
        }
      })

      # Show the pop-up window to change the model
      observeEvent(input$act_change_model, {
        tabset_est_par_result <- input$tabset_est_par_result

        # BSR
        if(tabset_est_par_result == "BSR"){
          if(!is.null(aics_bsr_all_mk)){
            row_index_select <- input$dt_result_bsr_rows_selected
            if(isTruthy(row_index_select)){
              aics <- aics_bsr_all_mk[[row_index_select]]
              if(!is.null(aics)){
                models_cand <- names(aics)

                bsr_slopes <- bsr_allele_median_slope_all_mk[row_index_select]
                bsr_intercepts <- bsr_allele_median_intercept_all_mk[row_index_select]
                bsr_sigmas <- bsr_allele_median_sigma_all_mk[row_index_select]
                bsr_min_us <- NA
                bsr_rhat_slopes <- bsr_allele_rhat_slope_all_mk[row_index_select]
                bsr_rhat_intercepts <- bsr_allele_rhat_intercept_all_mk[row_index_select]
                bsr_rhat_sigmas <- bsr_allele_rhat_sigma_all_mk[row_index_select]
                bsr_rhat_lps <- bsr_allele_rhat_lp_all_mk[row_index_select]

                if(is.element("LUS", models_cand)){
                  bsr_slopes <- c(bsr_slopes, bsr_lus_median_slope_all_mk[row_index_select])
                  bsr_intercepts <- c(bsr_intercepts, bsr_lus_median_intercept_all_mk[row_index_select])
                  bsr_sigmas <- c(bsr_sigmas, bsr_lus_median_sigma_all_mk[row_index_select])
                  bsr_min_us <- c(bsr_min_us, NA)
                  bsr_rhat_slopes <- c(bsr_rhat_slopes, bsr_lus_rhat_slope_all_mk[row_index_select])
                  bsr_rhat_intercepts <- c(bsr_rhat_intercepts, bsr_lus_rhat_intercept_all_mk[row_index_select])
                  bsr_rhat_sigmas <- c(bsr_rhat_sigmas, bsr_lus_rhat_sigma_all_mk[row_index_select])
                  bsr_rhat_lps <- c(bsr_rhat_lps, bsr_lus_rhat_lp_all_mk[row_index_select])
                }

                if(is.element("Multi-sequence", models_cand)){
                  bsr_slopes <- c(bsr_slopes, bsr_multi_median_slope_all_mk[row_index_select])
                  bsr_intercepts <- c(bsr_intercepts, bsr_multi_median_intercept_all_mk[row_index_select])
                  bsr_sigmas <- c(bsr_sigmas, bsr_multi_median_sigma_all_mk[row_index_select])
                  bsr_min_us <- c(bsr_min_us, bsr_multi_min_us_all_mk[row_index_select])
                  bsr_rhat_slopes <- c(bsr_rhat_slopes, bsr_multi_rhat_slope_all_mk[row_index_select])
                  bsr_rhat_intercepts <- c(bsr_rhat_intercepts, bsr_multi_rhat_intercept_all_mk[row_index_select])
                  bsr_rhat_sigmas <- c(bsr_rhat_sigmas, bsr_multi_rhat_sigma_all_mk[row_index_select])
                  bsr_rhat_lps <- c(bsr_rhat_lps, bsr_multi_rhat_lp_all_mk[row_index_select])
                }

                rv_dt_change_model$bsr <- data.table(Model = models_cand,
                                                     AIC = signif(aics, 3),
                                                     Slope = signif(bsr_slopes, 3),
                                                     Intercept = signif(bsr_intercepts, 3),
                                                     Sigma = signif(bsr_sigmas, 3),
                                                     MinUS = bsr_min_us,
                                                     Rhat_Slope = signif(bsr_rhat_slopes, 3),
                                                     Rhat_Intercept = signif(bsr_rhat_intercepts, 3),
                                                     Rhat_Sigma = signif(bsr_rhat_sigmas, 3),
                                                     Rhat_Likelihood = signif(bsr_rhat_lps, 3)
                                                     )

                showModal(modalDialog(
                  title = "Change model",
                  dataTableOutput(session$ns("dt_change_model_bsr")),
                  br(),
                  h5(div(paste0("Current model: ", bsr_best_model_all_mk[row_index_select]), style = "color:#555555;font-weight:bold;")),
                  br(),
                  selectInput(session$ns("select_adopt_model_bsr"), label = h4("Select a model"), choices = models_cand, selected = bsr_best_model_all_mk[row_index_select], multiple = FALSE),
                  footer = tagList(
                    actionButton(session$ns("act_adopt_model_bsr"), "Adopt"),
                    modalButton("Close")
                  ),
                  size = "l"
                ))
              }else if(bsr_best_model_all_mk[row_index_select] == "Multiple loci together"){
                showModal(modalDialog(title = "Information", "The model of the selected marker is 'multiple loci together' which cannot be changed to other models.", easyClose = TRUE, footer = NULL))
              }else{
                showModal(modalDialog(title = "Information", "The back stutter ratio is not considered in the selected marker.", easyClose = TRUE, footer = NULL))
              }
            }else{
              showModal(modalDialog(title = "Error", "Select a row!", easyClose = TRUE, footer = NULL))
            }
          }
        }

        # FSR
        else if(tabset_est_par_result == "FSR"){
          if(!is.null(aics_fsr_all_mk)){
            row_index_select <- input$dt_result_fsr_rows_selected
            if(isTruthy(row_index_select)){
              aics <- aics_fsr_all_mk[[row_index_select]]
              if(!is.null(aics)){
                models_cand <- names(aics)

                fsr_slopes <- fsr_allele_median_slope_all_mk[row_index_select]
                fsr_intercepts <- fsr_allele_median_intercept_all_mk[row_index_select]
                fsr_sigmas <- fsr_allele_median_sigma_all_mk[row_index_select]
                fsr_min_us <- NA
                fsr_rhat_slopes <- fsr_allele_rhat_slope_all_mk[row_index_select]
                fsr_rhat_intercepts <- fsr_allele_rhat_intercept_all_mk[row_index_select]
                fsr_rhat_sigmas <- fsr_allele_rhat_sigma_all_mk[row_index_select]
                fsr_rhat_lps <- fsr_allele_rhat_lp_all_mk[row_index_select]

                if(is.element("LUS", models_cand)){
                  fsr_slopes <- c(fsr_slopes, fsr_lus_median_slope_all_mk[row_index_select])
                  fsr_intercepts <- c(fsr_intercepts, fsr_lus_median_intercept_all_mk[row_index_select])
                  fsr_sigmas <- c(fsr_sigmas, fsr_lus_median_sigma_all_mk[row_index_select])
                  fsr_min_us <- c(fsr_min_us, NA)
                  fsr_rhat_slopes <- c(fsr_rhat_slopes, fsr_lus_rhat_slope_all_mk[row_index_select])
                  fsr_rhat_intercepts <- c(fsr_rhat_intercepts, fsr_lus_rhat_intercept_all_mk[row_index_select])
                  fsr_rhat_sigmas <- c(fsr_rhat_sigmas, fsr_lus_rhat_sigma_all_mk[row_index_select])
                  fsr_rhat_lps <- c(fsr_rhat_lps, fsr_lus_rhat_lp_all_mk[row_index_select])
                }

                if(is.element("Multi-sequence", models_cand)){
                  fsr_slopes <- c(fsr_slopes, fsr_multi_median_slope_all_mk[row_index_select])
                  fsr_intercepts <- c(fsr_intercepts, fsr_multi_median_intercept_all_mk[row_index_select])
                  fsr_sigmas <- c(fsr_sigmas, fsr_multi_median_sigma_all_mk[row_index_select])
                  fsr_min_us <- c(fsr_min_us, fsr_multi_min_us_all_mk[row_index_select])
                  fsr_rhat_slopes <- c(fsr_rhat_slopes, fsr_multi_rhat_slope_all_mk[row_index_select])
                  fsr_rhat_intercepts <- c(fsr_rhat_intercepts, fsr_multi_rhat_intercept_all_mk[row_index_select])
                  fsr_rhat_sigmas <- c(fsr_rhat_sigmas, fsr_multi_rhat_sigma_all_mk[row_index_select])
                  fsr_rhat_lps <- c(fsr_rhat_lps, fsr_multi_rhat_lp_all_mk[row_index_select])
                }

                rv_dt_change_model$fsr <- data.table(Model = models_cand,
                                                     AIC = signif(aics, 3),
                                                     Slope = signif(fsr_slopes, 3),
                                                     Intercept = signif(fsr_intercepts, 3),
                                                     Sigma = signif(fsr_sigmas, 3),
                                                     MinUS = fsr_min_us,
                                                     Rhat_Slope = signif(fsr_rhat_slopes, 3),
                                                     Rhat_Intercept = signif(fsr_rhat_intercepts, 3),
                                                     Rhat_Sigma = signif(fsr_rhat_sigmas, 3),
                                                     Rhat_Likelihood = signif(fsr_rhat_lps, 3)
                                                     )

                showModal(modalDialog(
                  title = "Change model",
                  dataTableOutput(session$ns("dt_change_model_fsr")),
                  br(),
                  h5(div(paste0("Current model: ", fsr_best_model_all_mk[row_index_select]), style = "color:#555555;font-weight:bold;")),
                  br(),
                  selectInput(session$ns("select_adopt_model_fsr"), label = h4("Select a model"), choices = models_cand, selected = fsr_best_model_all_mk[row_index_select], multiple = FALSE),
                  footer = tagList(
                    actionButton(session$ns("act_adopt_model_fsr"), "Adopt"),
                    modalButton("Close")
                  ),
                  size = "l"
                ))
              }else if(fsr_best_model_all_mk[row_index_select] == "Multiple loci together"){
                showModal(modalDialog(title = "Information", "The model of the selected marker is 'multiple loci together' which cannot be changed to other models.", easyClose = TRUE, footer = NULL))
              }else{
                showModal(modalDialog(title = "Information", "The forward stutter ratio is not considered in the selected marker.", easyClose = TRUE, footer = NULL))
              }
            }else{
              showModal(modalDialog(title = "Error", "Select a row!", easyClose = TRUE, footer = NULL))
            }
          }
        }

        # DSR
        else if(tabset_est_par_result == "DSR"){
          if(!is.null(aics_dsr_all_mk)){
            row_index_select <- input$dt_result_dsr_rows_selected
            if(isTruthy(row_index_select)){
              aics <- aics_dsr_all_mk[[row_index_select]]
              if(!is.null(aics)){
                models_cand <- names(aics)

                dsr_slopes <- dsr_allele_median_slope_all_mk[row_index_select]
                dsr_intercepts <- dsr_allele_median_intercept_all_mk[row_index_select]
                dsr_sigmas <- dsr_allele_median_sigma_all_mk[row_index_select]
                dsr_min_us <- NA
                dsr_rhat_slopes <- dsr_allele_rhat_slope_all_mk[row_index_select]
                dsr_rhat_intercepts <- dsr_allele_rhat_intercept_all_mk[row_index_select]
                dsr_rhat_sigmas <- dsr_allele_rhat_sigma_all_mk[row_index_select]
                dsr_rhat_lps <- dsr_allele_rhat_lp_all_mk[row_index_select]

                if(is.element("LUS", models_cand)){
                  dsr_slopes <- c(dsr_slopes, dsr_lus_median_slope_all_mk[row_index_select])
                  dsr_intercepts <- c(dsr_intercepts, dsr_lus_median_intercept_all_mk[row_index_select])
                  dsr_sigmas <- c(dsr_sigmas, dsr_lus_median_sigma_all_mk[row_index_select])
                  dsr_min_us <- c(dsr_min_us, NA)
                  dsr_rhat_slopes <- c(dsr_rhat_slopes, dsr_lus_rhat_slope_all_mk[row_index_select])
                  dsr_rhat_intercepts <- c(dsr_rhat_intercepts, dsr_lus_rhat_intercept_all_mk[row_index_select])
                  dsr_rhat_sigmas <- c(dsr_rhat_sigmas, dsr_lus_rhat_sigma_all_mk[row_index_select])
                  dsr_rhat_lps <- c(dsr_rhat_lps, dsr_lus_rhat_lp_all_mk[row_index_select])
                }

                if(is.element("Multi-sequence", models_cand)){
                  dsr_slopes <- c(dsr_slopes, dsr_multi_median_slope_all_mk[row_index_select])
                  dsr_intercepts <- c(dsr_intercepts, dsr_multi_median_intercept_all_mk[row_index_select])
                  dsr_sigmas <- c(dsr_sigmas, dsr_multi_median_sigma_all_mk[row_index_select])
                  dsr_min_us <- c(dsr_min_us, dsr_multi_min_us_all_mk[row_index_select])
                  dsr_rhat_slopes <- c(dsr_rhat_slopes, dsr_multi_rhat_slope_all_mk[row_index_select])
                  dsr_rhat_intercepts <- c(dsr_rhat_intercepts, dsr_multi_rhat_intercept_all_mk[row_index_select])
                  dsr_rhat_sigmas <- c(dsr_rhat_sigmas, dsr_multi_rhat_sigma_all_mk[row_index_select])
                  dsr_rhat_lps <- c(dsr_rhat_lps, dsr_multi_rhat_lp_all_mk[row_index_select])
                }

                rv_dt_change_model$dsr <- data.table(Model = models_cand,
                                                     AIC = signif(aics, 3),
                                                     Slope = signif(dsr_slopes, 3),
                                                     Intercept = signif(dsr_intercepts, 3),
                                                     Sigma = signif(dsr_sigmas, 3),
                                                     MinUS = dsr_min_us,
                                                     Rhat_Slope = signif(dsr_rhat_slopes, 3),
                                                     Rhat_Intercept = signif(dsr_rhat_intercepts, 3),
                                                     Rhat_Sigma = signif(dsr_rhat_sigmas, 3),
                                                     Rhat_Likelihood = signif(dsr_rhat_lps, 3)
                                                     )

                showModal(modalDialog(
                  title = "Change model",
                  dataTableOutput(session$ns("dt_change_model_dsr")),
                  br(),
                  h5(div(paste0("Current model: ", dsr_best_model_all_mk[row_index_select]), style = "color:#555555;font-weight:bold;")),
                  br(),
                  selectInput(session$ns("select_adopt_model_dsr"), label = h4("Select a model"), choices = models_cand, selected = dsr_best_model_all_mk[row_index_select], multiple = FALSE),
                  footer = tagList(
                    actionButton(session$ns("act_adopt_model_dsr"), "Adopt"),
                    modalButton("Close")
                  ),
                  size = "l"
                ))
              }else if(dsr_best_model_all_mk[row_index_select] == "Multiple loci together"){
                showModal(modalDialog(title = "Information", "The model of the selected marker is 'multiple loci together' which cannot be changed to other models.", easyClose = TRUE, footer = NULL))
              }else{
                showModal(modalDialog(title = "Information", "The double-back stutter ratio is not considered in the selected marker.", easyClose = TRUE, footer = NULL))
              }
            }else{
              showModal(modalDialog(title = "Error", "Select a row!", easyClose = TRUE, footer = NULL))
            }
          }
        }
      })

      # Display information on AIC values and parameters of each back stutter model
      output$dt_change_model_bsr <- renderDataTable(server = FALSE, {
        dt_change_model_bsr <- rv_dt_change_model$bsr
        datatable(
          dt_change_model_bsr,
          selection = "none",
          options = list(dom = "t", iDisplayLength = nrow(dt_change_model_bsr), ordering = FALSE, autowidth = TRUE),
          rownames = FALSE
        )
      })

      # Display information on AIC values and parameters of each forward stutter model
      output$dt_change_model_fsr <- renderDataTable(server = FALSE, {
        dt_change_model_fsr <- rv_dt_change_model$fsr
        datatable(
          dt_change_model_fsr,
          selection = "none",
          options = list(dom = "t", iDisplayLength = nrow(dt_change_model_fsr), ordering = FALSE, autowidth = TRUE),
          rownames = FALSE
        )
      })

      # Display information on AIC values and parameters of each double-back stutter model
      output$dt_change_model_dsr <- renderDataTable(server = FALSE, {
        dt_change_model_dsr <- rv_dt_change_model$dsr
        datatable(
          dt_change_model_dsr,
          selection = "none",
          options = list(dom = "t", iDisplayLength = nrow(dt_change_model_dsr), ordering = FALSE, autowidth = TRUE),
          rownames = FALSE
        )
      })

      # Adopt model for BSR (Update the result table and the trace plot of the best model)
      observeEvent(input$act_adopt_model_bsr, {
        row_index_select <- input$dt_result_bsr_rows_selected
        dt_result_bsr <- rv_results$dt_result_bsr
        adopt_model <- input$select_adopt_model_bsr
        dt_result_bsr$Model[row_index_select] <- adopt_model
        if(adopt_model == "Allele"){
          dt_result_bsr$Slope[row_index_select] <- bsr_allele_median_slope_all_mk[row_index_select]
          dt_result_bsr$Intercept[row_index_select] <- bsr_allele_median_intercept_all_mk[row_index_select]
          dt_result_bsr$Sigma[row_index_select] <- bsr_allele_median_sigma_all_mk[row_index_select]
          dt_result_bsr$MinUS[row_index_select] <- NA
          dt_result_bsr$Rhat_Slope[row_index_select] <- bsr_allele_rhat_slope_all_mk[row_index_select]
          dt_result_bsr$Rhat_Intercept[row_index_select] <- bsr_allele_rhat_intercept_all_mk[row_index_select]
          dt_result_bsr$Rhat_Sigma[row_index_select] <- bsr_allele_rhat_sigma_all_mk[row_index_select]
          dt_result_bsr$Rhat_Likelihood[row_index_select] <- bsr_allele_rhat_lp_all_mk[row_index_select]
          rv_results$dt_result_bsr <- dt_result_bsr
          rv_results$traceplot_bsr[[row_index_select]] <- bsr_allele_trace_all_mk[[row_index_select]]
        }else if(adopt_model == "LUS"){
          dt_result_bsr$Slope[row_index_select] <- bsr_lus_median_slope_all_mk[row_index_select]
          dt_result_bsr$Intercept[row_index_select] <- bsr_lus_median_intercept_all_mk[row_index_select]
          dt_result_bsr$Sigma[row_index_select] <- bsr_lus_median_sigma_all_mk[row_index_select]
          dt_result_bsr$MinUS[row_index_select] <- NA
          dt_result_bsr$Rhat_Slope[row_index_select] <- bsr_lus_rhat_slope_all_mk[row_index_select]
          dt_result_bsr$Rhat_Intercept[row_index_select] <- bsr_lus_rhat_intercept_all_mk[row_index_select]
          dt_result_bsr$Rhat_Sigma[row_index_select] <- bsr_lus_rhat_sigma_all_mk[row_index_select]
          dt_result_bsr$Rhat_Likelihood[row_index_select] <- bsr_lus_rhat_lp_all_mk[row_index_select]
          rv_results$dt_result_bsr <- dt_result_bsr
          rv_results$traceplot_bsr[[row_index_select]] <- bsr_lus_trace_all_mk[[row_index_select]]
        }else if(adopt_model == "Multi-sequence"){
          dt_result_bsr$Slope[row_index_select] <- bsr_multi_median_slope_all_mk[row_index_select]
          dt_result_bsr$Intercept[row_index_select] <- bsr_multi_median_intercept_all_mk[row_index_select]
          dt_result_bsr$Sigma[row_index_select] <- bsr_multi_median_sigma_all_mk[row_index_select]
          dt_result_bsr$MinUS[row_index_select] <- bsr_multi_min_us_all_mk[row_index_select]
          dt_result_bsr$Rhat_Slope[row_index_select] <- bsr_multi_rhat_slope_all_mk[row_index_select]
          dt_result_bsr$Rhat_Intercept[row_index_select] <- bsr_multi_rhat_intercept_all_mk[row_index_select]
          dt_result_bsr$Rhat_Sigma[row_index_select] <- bsr_multi_rhat_sigma_all_mk[row_index_select]
          dt_result_bsr$Rhat_Likelihood[row_index_select] <- bsr_multi_rhat_lp_all_mk[row_index_select]
          rv_results$dt_result_bsr <- dt_result_bsr
          rv_results$traceplot_bsr[[row_index_select]] <- bsr_multi_trace_all_mk[[row_index_select]]
        }
      })

      # Adopt model for FSR (Update the result table and the trace plot of the best model)
      observeEvent(input$act_adopt_model_fsr, {
        row_index_select <- input$dt_result_fsr_rows_selected
        dt_result_fsr <- rv_results$dt_result_fsr
        adopt_model <- input$select_adopt_model_fsr
        dt_result_fsr$Model[row_index_select] <- adopt_model
        if(adopt_model == "Allele"){
          dt_result_fsr$Slope[row_index_select] <- fsr_allele_median_slope_all_mk[row_index_select]
          dt_result_fsr$Intercept[row_index_select] <- fsr_allele_median_intercept_all_mk[row_index_select]
          dt_result_fsr$Sigma[row_index_select] <- fsr_allele_median_sigma_all_mk[row_index_select]
          dt_result_fsr$MinUS[row_index_select] <- NA
          dt_result_fsr$Rhat_Slope[row_index_select] <- fsr_allele_rhat_slope_all_mk[row_index_select]
          dt_result_fsr$Rhat_Intercept[row_index_select] <- fsr_allele_rhat_intercept_all_mk[row_index_select]
          dt_result_fsr$Rhat_Sigma[row_index_select] <- fsr_allele_rhat_sigma_all_mk[row_index_select]
          dt_result_fsr$Rhat_Likelihood[row_index_select] <- fsr_allele_rhat_lp_all_mk[row_index_select]
          rv_results$dt_result_fsr <- dt_result_fsr
          rv_results$traceplot_fsr[[row_index_select]] <- fsr_allele_trace_all_mk[[row_index_select]]
        }else if(adopt_model == "LUS"){
          dt_result_fsr$Slope[row_index_select] <- fsr_lus_median_slope_all_mk[row_index_select]
          dt_result_fsr$Intercept[row_index_select] <- fsr_lus_median_intercept_all_mk[row_index_select]
          dt_result_fsr$Sigma[row_index_select] <- fsr_lus_median_sigma_all_mk[row_index_select]
          dt_result_fsr$MinUS[row_index_select] <- NA
          dt_result_fsr$Rhat_Slope[row_index_select] <- fsr_lus_rhat_slope_all_mk[row_index_select]
          dt_result_fsr$Rhat_Intercept[row_index_select] <- fsr_lus_rhat_intercept_all_mk[row_index_select]
          dt_result_fsr$Rhat_Sigma[row_index_select] <- fsr_lus_rhat_sigma_all_mk[row_index_select]
          dt_result_fsr$Rhat_Likelihood[row_index_select] <- fsr_lus_rhat_lp_all_mk[row_index_select]
          rv_results$dt_result_fsr <- dt_result_fsr
          rv_results$traceplot_fsr[[row_index_select]] <- fsr_lus_trace_all_mk[[row_index_select]]
        }else if(adopt_model == "Multi-sequence"){
          dt_result_fsr$Slope[row_index_select] <- fsr_multi_median_slope_all_mk[row_index_select]
          dt_result_fsr$Intercept[row_index_select] <- fsr_multi_median_intercept_all_mk[row_index_select]
          dt_result_fsr$Sigma[row_index_select] <- fsr_multi_median_sigma_all_mk[row_index_select]
          dt_result_fsr$MinUS[row_index_select] <- fsr_multi_min_us_all_mk[row_index_select]
          dt_result_fsr$Rhat_Slope[row_index_select] <- fsr_multi_rhat_slope_all_mk[row_index_select]
          dt_result_fsr$Rhat_Intercept[row_index_select] <- fsr_multi_rhat_intercept_all_mk[row_index_select]
          dt_result_fsr$Rhat_Sigma[row_index_select] <- fsr_multi_rhat_sigma_all_mk[row_index_select]
          dt_result_fsr$Rhat_Likelihood[row_index_select] <- fsr_multi_rhat_lp_all_mk[row_index_select]
          rv_results$dt_result_fsr <- dt_result_fsr
          rv_results$traceplot_fsr[[row_index_select]] <- fsr_multi_trace_all_mk[[row_index_select]]
        }
      })

      # Adopt model for DSR (Update the result table and the trace plot of the best model)
      observeEvent(input$act_adopt_model_dsr, {
        row_index_select <- input$dt_result_dsr_rows_selected
        dt_result_dsr <- rv_results$dt_result_dsr
        adopt_model <- input$select_adopt_model_dsr
        dt_result_dsr$Model[row_index_select] <- adopt_model
        if(adopt_model == "Allele"){
          dt_result_dsr$Slope[row_index_select] <- dsr_allele_median_slope_all_mk[row_index_select]
          dt_result_dsr$Intercept[row_index_select] <- dsr_allele_median_intercept_all_mk[row_index_select]
          dt_result_dsr$Sigma[row_index_select] <- dsr_allele_median_sigma_all_mk[row_index_select]
          dt_result_dsr$MinUS[row_index_select] <- NA
          dt_result_dsr$Rhat_Slope[row_index_select] <- dsr_allele_rhat_slope_all_mk[row_index_select]
          dt_result_dsr$Rhat_Intercept[row_index_select] <- dsr_allele_rhat_intercept_all_mk[row_index_select]
          dt_result_dsr$Rhat_Sigma[row_index_select] <- dsr_allele_rhat_sigma_all_mk[row_index_select]
          dt_result_dsr$Rhat_Likelihood[row_index_select] <- dsr_allele_rhat_lp_all_mk[row_index_select]
          rv_results$dt_result_dsr <- dt_result_dsr
          rv_results$traceplot_dsr[[row_index_select]] <- dsr_allele_trace_all_mk[[row_index_select]]
        }else if(adopt_model == "LUS"){
          dt_result_dsr$Slope[row_index_select] <- dsr_lus_median_slope_all_mk[row_index_select]
          dt_result_dsr$Intercept[row_index_select] <- dsr_lus_median_intercept_all_mk[row_index_select]
          dt_result_dsr$Sigma[row_index_select] <- dsr_lus_median_sigma_all_mk[row_index_select]
          dt_result_dsr$MinUS[row_index_select] <- NA
          dt_result_dsr$Rhat_Slope[row_index_select] <- dsr_lus_rhat_slope_all_mk[row_index_select]
          dt_result_dsr$Rhat_Intercept[row_index_select] <- dsr_lus_rhat_intercept_all_mk[row_index_select]
          dt_result_dsr$Rhat_Sigma[row_index_select] <- dsr_lus_rhat_sigma_all_mk[row_index_select]
          dt_result_dsr$Rhat_Likelihood[row_index_select] <- dsr_lus_rhat_lp_all_mk[row_index_select]
          rv_results$dt_result_dsr <- dt_result_dsr
          rv_results$traceplot_dsr[[row_index_select]] <- dsr_lus_trace_all_mk[[row_index_select]]
        }else if(adopt_model == "Multi-sequence"){
          dt_result_dsr$Slope[row_index_select] <- dsr_multi_median_slope_all_mk[row_index_select]
          dt_result_dsr$Intercept[row_index_select] <- dsr_multi_median_intercept_all_mk[row_index_select]
          dt_result_dsr$Sigma[row_index_select] <- dsr_multi_median_sigma_all_mk[row_index_select]
          dt_result_dsr$MinUS[row_index_select] <- dsr_multi_min_us_all_mk[row_index_select]
          dt_result_dsr$Rhat_Slope[row_index_select] <- dsr_multi_rhat_slope_all_mk[row_index_select]
          dt_result_dsr$Rhat_Intercept[row_index_select] <- dsr_multi_rhat_intercept_all_mk[row_index_select]
          dt_result_dsr$Rhat_Sigma[row_index_select] <- dsr_multi_rhat_sigma_all_mk[row_index_select]
          dt_result_dsr$Rhat_Likelihood[row_index_select] <- dsr_multi_rhat_lp_all_mk[row_index_select]
          rv_results$dt_result_dsr <- dt_result_dsr
          rv_results$traceplot_dsr[[row_index_select]] <- dsr_multi_trace_all_mk[[row_index_select]]
        }
      })

      #### Trace plot ####

      # Define reactive values for displaying the trace plot
      rv_traceplot <- reactiveValues()
      rv_traceplot$ae <- NULL
      rv_traceplot$bsr <- NULL
      rv_traceplot$fsr <- NULL
      rv_traceplot$dsr <- NULL
      rv_traceplot$m2sr <- NULL

      # Show the pop-up window to display the trace plot
      observeEvent(input$act_traceplot, {
        tabset_est_par_result <- input$tabset_est_par_result

        # AE
        if(tabset_est_par_result == "AE"){
          if(!is.null(ae_trace_all_mk)){
            row_index_select <- input$dt_result_ae_rows_selected
            if(isTruthy(row_index_select)){
              rv_traceplot$ae <- ae_trace_all_mk[[row_index_select]]

              showModal(modalDialog(
                title = "Trace plot",
                plotOutput(session$ns("traceplot_ae")),
                footer = tagList(
                  modalButton("Close")
                ),
                size = "l"
              ))

            }else{
              showModal(modalDialog(title = "Error", "Select a row!", easyClose = TRUE, footer = NULL))
            }
          }
        }

        # BSR
        else if(tabset_est_par_result == "BSR"){
          if(!is.null(bsr_best_trace_all_mk)){
            row_index_select <- input$dt_result_bsr_rows_selected
            if(isTruthy(row_index_select)){
              isolate(dt_result_bsr <- rv_results$dt_result_bsr)
              if(dt_result_bsr$Model[row_index_select] != "Not consider"){
                rv_traceplot$bsr <- rv_results$traceplot_bsr[[row_index_select]]

                showModal(modalDialog(
                  title = "Trace plot",
                  plotOutput(session$ns("traceplot_bsr")),
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
        }

        # FSR
        else if(tabset_est_par_result == "FSR"){
          if(!is.null(fsr_best_trace_all_mk)){
            row_index_select <- input$dt_result_fsr_rows_selected
            if(isTruthy(row_index_select)){
              isolate(dt_result_fsr <- rv_results$dt_result_fsr)
              if(dt_result_fsr$Model[row_index_select] != "Not consider"){
                rv_traceplot$fsr <- rv_results$traceplot_fsr[[row_index_select]]

                showModal(modalDialog(
                  title = "Trace plot",
                  plotOutput(session$ns("traceplot_fsr")),
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
        }

        # DSR
        else if(tabset_est_par_result == "DSR"){
          if(!is.null(dsr_best_trace_all_mk)){
            row_index_select <- input$dt_result_dsr_rows_selected
            if(isTruthy(row_index_select)){
              isolate(dt_result_dsr <- rv_results$dt_result_dsr)
              if(dt_result_dsr$Model[row_index_select] != "Not consider"){
                rv_traceplot$dsr <- rv_results$traceplot_dsr[[row_index_select]]

                showModal(modalDialog(
                  title = "Trace plot",
                  plotOutput(session$ns("traceplot_dsr")),
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
        }

        # M2SR
        else if(tabset_est_par_result == "M2SR"){
          if(!is.null(m2sr_best_trace_all_mk)){
            row_index_select <- input$dt_result_m2sr_rows_selected
            if(isTruthy(row_index_select)){
              isolate(dt_result_m2sr <- rv_results$dt_result_m2sr)
              if(dt_result_m2sr$Model[row_index_select] != "Not consider"){
                rv_traceplot$m2sr <- m2sr_best_trace_all_mk[[row_index_select]]

                showModal(modalDialog(
                  title = "Trace plot",
                  plotOutput(session$ns("traceplot_m2sr")),
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
        }
      })

      # Display the trace plot for AE
      output$traceplot_ae <- renderPlot({
        rv_traceplot$ae
      })

      # Display the trace plot for BSR
      output$traceplot_bsr <- renderPlot({
        rv_traceplot$bsr
      })

      # Display the trace plot for FSR
      output$traceplot_fsr <- renderPlot({
        rv_traceplot$fsr
      })

      # Display the trace plot for DSR
      output$traceplot_dsr <- renderPlot({
        rv_traceplot$dsr
      })

      # Display the trace plot for M2SR
      output$traceplot_m2sr <- renderPlot({
        rv_traceplot$m2sr
      })

      #### Data plot ####

      # Define reactive values for displaying the data plot
      rv_dataplot <- reactiveValues()
      rv_dataplot$bsr_xval <- NULL
      rv_dataplot$bsr_yval <- NULL
      rv_dataplot$fsr_xval <- NULL
      rv_dataplot$fsr_yval <- NULL
      rv_dataplot$dsr_xval <- NULL
      rv_dataplot$dsr_yval <- NULL
      rv_dataplot$m2sr_xval <- NULL
      rv_dataplot$m2sr_yval <- NULL

      # Show the pop-up window to display the data plot
      observeEvent(input$act_dataplot, {
        tabset_est_par_result <- input$tabset_est_par_result

        # AE
        if(tabset_est_par_result == "AE"){
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
        }

        # BSR
        else if(tabset_est_par_result == "BSR"){
          if(!is.null(bsr_data_all_mk)){
            row_index_select <- input$dt_result_bsr_rows_selected
            if(isTruthy(row_index_select)){
              isolate(dt_result_bsr <- rv_results$dt_result_bsr)
              model_bsr <- dt_result_bsr$Model[row_index_select]
              if(model_bsr != "Not consider"){
                bsr_data_selected <- bsr_data_all_mk[[row_index_select]]
                if(model_bsr == "LUS"){
                  rv_dataplot$bsr_xval <- bsr_data_selected$bsr_lus
                }else if(model_bsr == "Multi-sequence"){
                  rv_dataplot$bsr_xval <- bsr_data_selected$bsr_ca
                }else{
                  rv_dataplot$bsr_xval <- bsr_data_selected$bsr_al
                }
                rv_dataplot$bsr_yval <- bsr_data_selected$bsr

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
        }

        # FSR
        else if(tabset_est_par_result == "FSR"){
          if(!is.null(fsr_data_all_mk)){
            row_index_select <- input$dt_result_fsr_rows_selected
            if(isTruthy(row_index_select)){
              isolate(dt_result_fsr <- rv_results$dt_result_fsr)
              model_fsr <- dt_result_fsr$Model[row_index_select]
              if(model_fsr != "Not consider"){
                fsr_data_selected <- fsr_data_all_mk[[row_index_select]]
                if(model_fsr == "LUS"){
                  rv_dataplot$fsr_xval <- fsr_data_selected$fsr_lus
                }else if(model_fsr == "Multi-sequence"){
                  rv_dataplot$fsr_xval <- fsr_data_selected$fsr_ca
                }else{
                  rv_dataplot$fsr_xval <- fsr_data_selected$fsr_al
                }
                rv_dataplot$fsr_yval <- fsr_data_selected$fsr

                showModal(modalDialog(
                  title = "Data plot",
                  plotOutput(session$ns("dataplot_fsr")),
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
        }

        # DSR
        else if(tabset_est_par_result == "DSR"){
          if(!is.null(dsr_data_all_mk)){
            row_index_select <- input$dt_result_dsr_rows_selected
            if(isTruthy(row_index_select)){
              isolate(dt_result_dsr <- rv_results$dt_result_dsr)
              model_dsr <- dt_result_dsr$Model[row_index_select]
              if(model_dsr != "Not consider"){
                dsr_data_selected <- dsr_data_all_mk[[row_index_select]]
                if(model_dsr == "LUS"){
                  rv_dataplot$dsr_xval <- dsr_data_selected$dsr_lus
                }else if(model_dsr == "Multi-sequence"){
                  rv_dataplot$dsr_xval <- dsr_data_selected$dsr_ca
                }else{
                  rv_dataplot$dsr_xval <- dsr_data_selected$dsr_al
                }
                rv_dataplot$dsr_yval <- dsr_data_selected$dsr

                showModal(modalDialog(
                  title = "Data plot",
                  plotOutput(session$ns("dataplot_dsr")),
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
        }

        # M2SR
        else if(tabset_est_par_result == "M2SR"){
          if(!is.null(m2sr_data_all_mk)){
            row_index_select <- input$dt_result_m2sr_rows_selected
            if(isTruthy(row_index_select)){
              isolate(dt_result_m2sr <- rv_results$dt_result_m2sr)
              model_m2sr <- dt_result_m2sr$Model[row_index_select]
              if(model_m2sr != "Not consider"){
                m2sr_data_selected <- m2sr_data_all_mk[[row_index_select]]
                rv_dataplot$m2sr_xval <- m2sr_data_selected$m2sr_al
                rv_dataplot$m2sr_yval <- m2sr_data_selected$m2sr

                showModal(modalDialog(
                  title = "Data plot",
                  plotOutput(session$ns("dataplot_m2sr")),
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
        }
      })

      # Display the data plot for AE
      output$dataplot_ae <- renderPlot({
        boxplot(ae_data)
      })

      # Display the data plot for BSR
      output$dataplot_bsr <- renderPlot({
        plot(rv_dataplot$bsr_xval, rv_dataplot$bsr_yval)
      })

      # Display the data plot for FSR
      output$dataplot_fsr <- renderPlot({
        plot(rv_dataplot$fsr_xval, rv_dataplot$fsr_yval)
      })

      # Display the data plot for DSR
      output$dataplot_dsr <- renderPlot({
        plot(rv_dataplot$dsr_xval, rv_dataplot$dsr_yval)
      })

      # Display the data plot for M2SR
      output$dataplot_m2sr <- renderPlot({
        plot(rv_dataplot$m2sr_xval, rv_dataplot$m2sr_yval)
      })

      #### Result table ####

      # Display the result for AE
      output$dt_result_ae <- renderDataTable(server = FALSE, {
        datatable(
          rv_results$dt_result_ae,
          selection = list(mode = "single", target = "row"),
          options = list(iDisplayLength = 50, ordering = FALSE, autowidth = TRUE),
          rownames = FALSE
        )
      })

      # Display the result for BSR
      output$dt_result_bsr <- renderDataTable(server = FALSE, {
        datatable(
          rv_results$dt_result_bsr,
          selection = list(mode = "single", target = "row"),
          options = list(iDisplayLength = 50, ordering = FALSE, autowidth = TRUE),
          rownames = FALSE
        )
      })

      # Display the result for FSR
      output$dt_result_fsr <- renderDataTable(server = FALSE, {
        datatable(
          rv_results$dt_result_fsr,
          selection = list(mode = "single", target = "row"),
          options = list(iDisplayLength = 50, ordering = FALSE, autowidth = TRUE),
          rownames = FALSE
        )
      })

      # Display the result for DSR
      output$dt_result_dsr <- renderDataTable(server = FALSE, {
        datatable(
          rv_results$dt_result_dsr,
          selection = list(mode = "single", target = "row"),
          options = list(iDisplayLength = 50, ordering = FALSE, autowidth = TRUE),
          rownames = FALSE
        )
      })

      # Display the result for M2SR
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
