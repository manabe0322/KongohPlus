est_par_result_ui <- function(id){
  ns <- NS(id)

  tabPanel("Result",
           br(),
           sidebarLayout(
             sidebarPanel(width = 2,
                          actionButton(ns("act_traceplot"), label = "Trace plot"),
                          br(),
                          br(),
                          br(),
                          actionButton(ns("act_dataplot"), label = "Data plot"),
                          br(),
                          br(),
                          br(),
                          disabled(actionButton(ns("act_alcor"), label = "Corrected allele")),
                          br(),
                          br(),
                          br(),
                          actionButton(ns("act_aic_model"), label = "AIC"),
                          br(),
                          br(),
                          br(),
                          actionButton(ns("act_save_par"), label = "Save", class = "btn btn-primary btn-lg")
             ),
             mainPanel(width = 10,
                       tabsetPanel(id = ns("tabset_est_par_result"),
                                   type = "tabs",
                                   tabPanel("AE",
                                            br(),
                                            dataTableOutput(ns("dt_result_ae")),
                                            br()
                                   ),
                                   tabPanel("BSR",
                                            br(),
                                            dataTableOutput(ns("dt_result_bsr")),
                                            br()
                                   ),
                                   tabPanel("FSR",
                                            br(),
                                            dataTableOutput(ns("dt_result_fsr")),
                                            br()
                                   ),
                                   tabPanel("DSR",
                                            br(),
                                            dataTableOutput(ns("dt_result_dsr")),
                                            br()
                                   ),
                                   tabPanel("M2SR",
                                            br(),
                                            dataTableOutput(ns("dt_result_m2sr")),
                                            br()
                                   ),
                                   tabPanel("Gamma",
                                            br(),
                                            dataTableOutput(ns("dt_result_gamma")),
                                            br()
                                   )
                       )
             )
           )
  )
}

est_par_result_server <- function(id, path_pack, rv_output){
  moduleServer(
    id,
    function(input, output, session){
      results_est_par <- rv_output$results_est_par

      if(is.null(results_est_par)){
        markers_par <- NULL
        model_bsr_all_mk <- NULL
        model_fsr_all_mk <- NULL
        model_dsr_all_mk <- NULL

        ae_data <- NULL
        ae_median_myu_all_mk <- NULL
        ae_median_sigma_all_mk <- NULL
        ae_rhat_myu_all_mk <- NULL
        ae_rhat_sigma_all_mk <- NULL
        ae_rhat_lp_all_mk <- NULL
        ae_trace_all_mk <- NULL
        dt_result_ae <- NULL

        bsr_data_all_mk <- NULL
        dt_alcor_bsr <- NULL
        aics_bsr_all_mk <- NULL
        bsr_allele_median_slope_all_mk <- NULL
        bsr_allele_median_intercept_all_mk <- NULL
        bsr_allele_median_sigma_all_mk <- NULL
        bsr_allele_rhat_slope_all_mk <- NULL
        bsr_allele_rhat_intercept_all_mk <- NULL
        bsr_allele_rhat_sigma_all_mk <- NULL
        bsr_allele_rhat_lp_all_mk <- NULL
        bsr_lus_median_slope_all_mk <- NULL
        bsr_lus_median_intercept_all_mk <- NULL
        bsr_lus_median_sigma_all_mk <- NULL
        bsr_lus_rhat_slope_all_mk <- NULL
        bsr_lus_rhat_intercept_all_mk <- NULL
        bsr_lus_rhat_sigma_all_mk <- NULL
        bsr_lus_rhat_lp_all_mk <- NULL
        bsr_multi_median_slope_all_mk <- NULL
        bsr_multi_median_intercept_all_mk <- NULL
        bsr_multi_median_sigma_all_mk <- NULL
        bsr_multi_min_us_all_mk <- NULL
        bsr_multi_rhat_slope_all_mk <- NULL
        bsr_multi_rhat_intercept_all_mk <- NULL
        bsr_multi_rhat_sigma_all_mk <- NULL
        bsr_multi_rhat_lp_all_mk <- NULL
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
        dt_result_bsr <- NULL

        fsr_data_all_mk <- NULL
        dt_alcor_fsr <- NULL
        aics_fsr_all_mk <- NULL
        fsr_allele_median_slope_all_mk <- NULL
        fsr_allele_median_intercept_all_mk <- NULL
        fsr_allele_median_sigma_all_mk <- NULL
        fsr_allele_rhat_slope_all_mk <- NULL
        fsr_allele_rhat_intercept_all_mk <- NULL
        fsr_allele_rhat_sigma_all_mk <- NULL
        fsr_allele_rhat_lp_all_mk <- NULL
        fsr_lus_median_slope_all_mk <- NULL
        fsr_lus_median_intercept_all_mk <- NULL
        fsr_lus_median_sigma_all_mk <- NULL
        fsr_lus_rhat_slope_all_mk <- NULL
        fsr_lus_rhat_intercept_all_mk <- NULL
        fsr_lus_rhat_sigma_all_mk <- NULL
        fsr_lus_rhat_lp_all_mk <- NULL
        fsr_multi_median_slope_all_mk <- NULL
        fsr_multi_median_intercept_all_mk <- NULL
        fsr_multi_median_sigma_all_mk <- NULL
        fsr_multi_min_us_all_mk <- NULL
        fsr_multi_rhat_slope_all_mk <- NULL
        fsr_multi_rhat_intercept_all_mk <- NULL
        fsr_multi_rhat_sigma_all_mk <- NULL
        fsr_multi_rhat_lp_all_mk <- NULL
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
        dt_result_fsr <- NULL

        dsr_data_all_mk <- NULL
        dt_alcor_dsr <- NULL
        aics_dsr_all_mk <- NULL
        dsr_allele_median_slope_all_mk <- NULL
        dsr_allele_median_intercept_all_mk <- NULL
        dsr_allele_median_sigma_all_mk <- NULL
        dsr_allele_rhat_slope_all_mk <- NULL
        dsr_allele_rhat_intercept_all_mk <- NULL
        dsr_allele_rhat_sigma_all_mk <- NULL
        dsr_allele_rhat_lp_all_mk <- NULL
        dsr_lus_median_slope_all_mk <- NULL
        dsr_lus_median_intercept_all_mk <- NULL
        dsr_lus_median_sigma_all_mk <- NULL
        dsr_lus_rhat_slope_all_mk <- NULL
        dsr_lus_rhat_intercept_all_mk <- NULL
        dsr_lus_rhat_sigma_all_mk <- NULL
        dsr_lus_rhat_lp_all_mk <- NULL
        dsr_multi_median_slope_all_mk <- NULL
        dsr_multi_median_intercept_all_mk <- NULL
        dsr_multi_median_sigma_all_mk <- NULL
        dsr_multi_min_us_all_mk <- NULL
        dsr_multi_rhat_slope_all_mk <- NULL
        dsr_multi_rhat_intercept_all_mk <- NULL
        dsr_multi_rhat_sigma_all_mk <- NULL
        dsr_multi_rhat_lp_all_mk <- NULL
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
        dt_result_dsr <- NULL

        m2sr_data_all_mk <- NULL
        m2sr_best_model_all_mk <- NULL
        m2sr_best_median_myu_all_mk <- NULL
        m2sr_best_median_sigma_all_mk <- NULL
        m2sr_best_rhat_myu_all_mk <- NULL
        m2sr_best_rhat_sigma_all_mk <- NULL
        m2sr_best_rhat_lp_all_mk <- NULL
        m2sr_best_trace_all_mk <- NULL
        dt_result_m2sr <- NULL

        gamma_median_beta <- NULL
        gamma_rhat_beta <- NULL
        gamma_rhat_lp <- NULL
        gamma_traceplot <- NULL
        dt_result_gamma <- NULL

        dt_par <- NULL
      }else{
        markers_par <- results_est_par$markers_par
        model_bsr_all_mk <- results_est_par$model_bsr_all_mk
        model_fsr_all_mk <- results_est_par$model_fsr_all_mk
        model_dsr_all_mk <- results_est_par$model_dsr_all_mk

        ae_data <- results_est_par$ae_data
        ae_median_myu_all_mk <- results_est_par$ae_median_myu_all_mk
        ae_median_sigma_all_mk <- results_est_par$ae_median_sigma_all_mk
        ae_rhat_myu_all_mk <- results_est_par$ae_rhat_myu_all_mk
        ae_rhat_sigma_all_mk <- results_est_par$ae_rhat_sigma_all_mk
        ae_rhat_lp_all_mk <- results_est_par$ae_rhat_lp_all_mk
        ae_trace_all_mk <- results_est_par$ae_trace_all_mk
        dt_result_ae <- data.table(Marker = markers_par,
                                   Myu = signif(ae_median_myu_all_mk, 3),
                                   Sigma = signif(ae_median_sigma_all_mk, 3),
                                   Rhat_Myu = signif(ae_rhat_myu_all_mk, 3),
                                   Rhat_Sigma = signif(ae_rhat_sigma_all_mk, 3),
                                   Rhat_Likelihood = signif(ae_rhat_lp_all_mk, 3)
        )

        bsr_data_all_mk <- results_est_par$bsr_data_all_mk
        dt_alcor_bsr <- results_est_par$dt_alcor_bsr
        aics_bsr_all_mk <- results_est_par$aics_bsr_all_mk
        bsr_allele_median_slope_all_mk <- results_est_par$bsr_allele_median_slope_all_mk
        bsr_allele_median_intercept_all_mk <- results_est_par$bsr_allele_median_intercept_all_mk
        bsr_allele_median_sigma_all_mk <- results_est_par$bsr_allele_median_sigma_all_mk
        bsr_allele_rhat_slope_all_mk <- results_est_par$bsr_allele_rhat_slope_all_mk
        bsr_allele_rhat_intercept_all_mk <- results_est_par$bsr_allele_rhat_intercept_all_mk
        bsr_allele_rhat_sigma_all_mk <- results_est_par$bsr_allele_rhat_sigma_all_mk
        bsr_allele_rhat_lp_all_mk <- results_est_par$bsr_allele_rhat_lp_all_mk
        bsr_lus_median_slope_all_mk <- results_est_par$bsr_lus_median_slope_all_mk
        bsr_lus_median_intercept_all_mk <- results_est_par$bsr_lus_median_intercept_all_mk
        bsr_lus_median_sigma_all_mk <- results_est_par$bsr_lus_median_sigma_all_mk
        bsr_lus_rhat_slope_all_mk <- results_est_par$bsr_lus_rhat_slope_all_mk
        bsr_lus_rhat_intercept_all_mk <- results_est_par$bsr_lus_rhat_intercept_all_mk
        bsr_lus_rhat_sigma_all_mk <- results_est_par$bsr_lus_rhat_sigma_all_mk
        bsr_lus_rhat_lp_all_mk <- results_est_par$bsr_lus_rhat_lp_all_mk
        bsr_multi_median_slope_all_mk <- results_est_par$bsr_multi_median_slope_all_mk
        bsr_multi_median_intercept_all_mk <- results_est_par$bsr_multi_median_intercept_all_mk
        bsr_multi_median_sigma_all_mk <- results_est_par$bsr_multi_median_sigma_all_mk
        bsr_multi_min_us_all_mk <- results_est_par$bsr_multi_min_us_all_mk
        bsr_multi_rhat_slope_all_mk <- results_est_par$bsr_multi_rhat_slope_all_mk
        bsr_multi_rhat_intercept_all_mk <- results_est_par$bsr_multi_rhat_intercept_all_mk
        bsr_multi_rhat_sigma_all_mk <- results_est_par$bsr_multi_rhat_sigma_all_mk
        bsr_multi_rhat_lp_all_mk <- results_est_par$bsr_multi_rhat_lp_all_mk
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
        dt_result_bsr <- data.table(Marker = markers_par,
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

        fsr_data_all_mk <- results_est_par$fsr_data_all_mk
        dt_alcor_fsr <- results_est_par$dt_alcor_fsr
        aics_fsr_all_mk <- results_est_par$aics_fsr_all_mk
        fsr_allele_median_slope_all_mk <- results_est_par$fsr_allele_median_slope_all_mk
        fsr_allele_median_intercept_all_mk <- results_est_par$fsr_allele_median_intercept_all_mk
        fsr_allele_median_sigma_all_mk <- results_est_par$fsr_allele_median_sigma_all_mk
        fsr_allele_rhat_slope_all_mk <- results_est_par$fsr_allele_rhat_slope_all_mk
        fsr_allele_rhat_intercept_all_mk <- results_est_par$fsr_allele_rhat_intercept_all_mk
        fsr_allele_rhat_sigma_all_mk <- results_est_par$fsr_allele_rhat_sigma_all_mk
        fsr_allele_rhat_lp_all_mk <- results_est_par$fsr_allele_rhat_lp_all_mk
        fsr_lus_median_slope_all_mk <- results_est_par$fsr_lus_median_slope_all_mk
        fsr_lus_median_intercept_all_mk <- results_est_par$fsr_lus_median_intercept_all_mk
        fsr_lus_median_sigma_all_mk <- results_est_par$fsr_lus_median_sigma_all_mk
        fsr_lus_rhat_slope_all_mk <- results_est_par$fsr_lus_rhat_slope_all_mk
        fsr_lus_rhat_intercept_all_mk <- results_est_par$fsr_lus_rhat_intercept_all_mk
        fsr_lus_rhat_sigma_all_mk <- results_est_par$fsr_lus_rhat_sigma_all_mk
        fsr_lus_rhat_lp_all_mk <- results_est_par$fsr_lus_rhat_lp_all_mk
        fsr_multi_median_slope_all_mk <- results_est_par$fsr_multi_median_slope_all_mk
        fsr_multi_median_intercept_all_mk <- results_est_par$fsr_multi_median_intercept_all_mk
        fsr_multi_median_sigma_all_mk <- results_est_par$fsr_multi_median_sigma_all_mk
        fsr_multi_min_us_all_mk <- results_est_par$fsr_multi_min_us_all_mk
        fsr_multi_rhat_slope_all_mk <- results_est_par$fsr_multi_rhat_slope_all_mk
        fsr_multi_rhat_intercept_all_mk <- results_est_par$fsr_multi_rhat_intercept_all_mk
        fsr_multi_rhat_sigma_all_mk <- results_est_par$fsr_multi_rhat_sigma_all_mk
        fsr_multi_rhat_lp_all_mk <- results_est_par$fsr_multi_rhat_lp_all_mk
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
        dt_result_fsr <- data.table(Marker = markers_par,
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

        dsr_data_all_mk <- results_est_par$dsr_data_all_mk
        dt_alcor_dsr <- results_est_par$dt_alcor_dsr
        aics_dsr_all_mk <- results_est_par$aics_dsr_all_mk
        dsr_allele_median_slope_all_mk <- results_est_par$dsr_allele_median_slope_all_mk
        dsr_allele_median_intercept_all_mk <- results_est_par$dsr_allele_median_intercept_all_mk
        dsr_allele_median_sigma_all_mk <- results_est_par$dsr_allele_median_sigma_all_mk
        dsr_allele_rhat_slope_all_mk <- results_est_par$dsr_allele_rhat_slope_all_mk
        dsr_allele_rhat_intercept_all_mk <- results_est_par$dsr_allele_rhat_intercept_all_mk
        dsr_allele_rhat_sigma_all_mk <- results_est_par$dsr_allele_rhat_sigma_all_mk
        dsr_allele_rhat_lp_all_mk <- results_est_par$dsr_allele_rhat_lp_all_mk
        dsr_lus_median_slope_all_mk <- results_est_par$dsr_lus_median_slope_all_mk
        dsr_lus_median_intercept_all_mk <- results_est_par$dsr_lus_median_intercept_all_mk
        dsr_lus_median_sigma_all_mk <- results_est_par$dsr_lus_median_sigma_all_mk
        dsr_lus_rhat_slope_all_mk <- results_est_par$dsr_lus_rhat_slope_all_mk
        dsr_lus_rhat_intercept_all_mk <- results_est_par$dsr_lus_rhat_intercept_all_mk
        dsr_lus_rhat_sigma_all_mk <- results_est_par$dsr_lus_rhat_sigma_all_mk
        dsr_lus_rhat_lp_all_mk <- results_est_par$dsr_lus_rhat_lp_all_mk
        dsr_multi_median_slope_all_mk <- results_est_par$dsr_multi_median_slope_all_mk
        dsr_multi_median_intercept_all_mk <- results_est_par$dsr_multi_median_intercept_all_mk
        dsr_multi_median_sigma_all_mk <- results_est_par$dsr_multi_median_sigma_all_mk
        dsr_multi_min_us_all_mk <- results_est_par$dsr_multi_min_us_all_mk
        dsr_multi_rhat_slope_all_mk <- results_est_par$dsr_multi_rhat_slope_all_mk
        dsr_multi_rhat_intercept_all_mk <- results_est_par$dsr_multi_rhat_intercept_all_mk
        dsr_multi_rhat_sigma_all_mk <- results_est_par$dsr_multi_rhat_sigma_all_mk
        dsr_multi_rhat_lp_all_mk <- results_est_par$dsr_multi_rhat_lp_all_mk
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
        dt_result_dsr <- data.table(Marker = markers_par,Model = dsr_best_model_all_mk,
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

        m2sr_data_all_mk <- results_est_par$m2sr_data_all_mk
        m2sr_best_model_all_mk <- results_est_par$m2sr_best_model_all_mk
        m2sr_best_median_myu_all_mk <- results_est_par$m2sr_best_median_myu_all_mk
        m2sr_best_median_sigma_all_mk <- results_est_par$m2sr_best_median_sigma_all_mk
        m2sr_best_rhat_myu_all_mk <- results_est_par$m2sr_best_rhat_myu_all_mk
        m2sr_best_rhat_sigma_all_mk <- results_est_par$m2sr_best_rhat_sigma_all_mk
        m2sr_best_rhat_lp_all_mk <- results_est_par$m2sr_best_rhat_lp_all_mk
        m2sr_best_trace_all_mk <- results_est_par$m2sr_best_trace_all_mk
        dt_result_m2sr <- data.table(Marker = markers_par,
                                     Model = m2sr_best_model_all_mk,
                                     Myu = signif(m2sr_best_median_myu_all_mk, 3),
                                     Sigma = signif(m2sr_best_median_sigma_all_mk, 3),
                                     Rhat_Myu = signif(m2sr_best_rhat_myu_all_mk, 3),
                                     Rhat_Sigma = signif(m2sr_best_rhat_sigma_all_mk, 3),
                                     Rhat_Likelihood = signif(m2sr_best_rhat_lp_all_mk, 3)
        )

        gamma_median_beta <- results_est_par$gamma_median_beta
        gamma_rhat_beta <- results_est_par$gamma_rhat_beta
        gamma_rhat_lp <- results_est_par$gamma_rhat_lp
        gamma_traceplot <- results_est_par$gamma_traceplot
        dt_result_gamma <- data.table(Beta = signif(gamma_median_beta, 3),
                                      Rhat_Beta = signif(gamma_rhat_beta, 3),
                                      Rhat_Likelihood = signif(gamma_rhat_lp, 3)
        )

        n_mk <- length(markers_par)
        dt_par <- data.table(Marker = c(rep(markers_par, 5), NA),
                             Parameter = c(rep("AE", n_mk), rep("BSR", n_mk), rep("FSR", n_mk), rep("DSR", n_mk), rep("M2SR", n_mk), "Gamma"),
                             Model = c(rep("Log-normal", n_mk), bsr_best_model_all_mk, fsr_best_model_all_mk, dsr_best_model_all_mk, m2sr_best_model_all_mk, "Gamma"),
                             Myu = c(ae_median_myu_all_mk, bsr_best_median_myu_all_mk, fsr_best_median_myu_all_mk, dsr_best_median_myu_all_mk, m2sr_best_median_myu_all_mk, NA),
                             Slope = c(rep(NA, n_mk), bsr_best_median_slope_all_mk, fsr_best_median_slope_all_mk, dsr_best_median_slope_all_mk, rep(NA, n_mk), NA),
                             Intercept = c(rep(NA, n_mk), bsr_best_median_intercept_all_mk, fsr_best_median_intercept_all_mk, dsr_best_median_intercept_all_mk, rep(NA, n_mk), NA),
                             Sigma = c(ae_median_sigma_all_mk, bsr_best_median_sigma_all_mk, fsr_best_median_sigma_all_mk, dsr_best_median_sigma_all_mk, m2sr_best_median_sigma_all_mk, NA),
                             MinUS = c(rep(NA, n_mk), bsr_best_min_us_all_mk, fsr_best_min_us_all_mk, dsr_best_min_us_all_mk, rep(NA, n_mk), NA),
                             Beta = c(rep(NA, 5 * n_mk), gamma_median_beta),
                             Rhat_Myu = c(ae_rhat_myu_all_mk, bsr_best_rhat_myu_all_mk, fsr_best_rhat_myu_all_mk, dsr_best_rhat_myu_all_mk, m2sr_best_rhat_myu_all_mk, NA),
                             Rhat_Slope = c(rep(NA, n_mk), bsr_best_rhat_slope_all_mk, fsr_best_rhat_slope_all_mk, dsr_best_rhat_slope_all_mk, rep(NA, n_mk), NA),
                             Rhat_Intercept = c(rep(NA, n_mk), bsr_best_median_intercept_all_mk, fsr_best_median_intercept_all_mk, dsr_best_median_intercept_all_mk, rep(NA, n_mk), NA),
                             Rhat_Sigma = c(ae_rhat_sigma_all_mk, bsr_best_median_sigma_all_mk, fsr_best_median_sigma_all_mk, dsr_best_median_sigma_all_mk, m2sr_best_median_sigma_all_mk, NA),
                             Rhat_Beta = c(rep(NA, 5 * n_mk), gamma_rhat_beta),
                             Rhat_Likelihood = c(ae_rhat_lp_all_mk, bsr_best_rhat_lp_all_mk, fsr_best_rhat_lp_all_mk, dsr_best_rhat_lp_all_mk, m2sr_best_rhat_lp_all_mk, gamma_rhat_lp)
        )
      }

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
              if(dt_result_bsr$Model[row_index_select] != "Not consider"){
                rv_traceplot$bsr <- bsr_best_trace_all_mk[[row_index_select]]

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
              if(dt_result_fsr$Model[row_index_select] != "Not consider"){
                rv_traceplot$fsr <- fsr_best_trace_all_mk[[row_index_select]]

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
              if(dt_result_dsr$Model[row_index_select] != "Not consider"){
                rv_traceplot$dsr <- dsr_best_trace_all_mk[[row_index_select]]

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

        # Gamma
        else if(tabset_est_par_result == "Gamma"){
          if(!is.null(gamma_traceplot)){
            showModal(modalDialog(
              title = "Trace plot",
              plotOutput(session$ns("traceplot_gamma")),
              footer = tagList(
                modalButton("Close")
              ),
              size = "l"
            ))
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

      # Display the trace plot for the gamma parameter
      output$traceplot_gamma <- renderPlot({
        gamma_traceplot
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
              adopt_model_bsr <- dt_result_bsr$Model[row_index_select]
              if(adopt_model_bsr != "Not consider"){
                bsr_data_selected <- bsr_data_all_mk[[row_index_select]]
                if(adopt_model_bsr == "LUS"){
                  rv_dataplot$bsr_xval <- bsr_data_selected$bsr_lus
                }else if(adopt_model_bsr == "Multi-sequence"){
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
              adopt_model_fsr <- dt_result_fsr$Model[row_index_select]
              if(adopt_model_fsr != "Not consider"){
                fsr_data_selected <- fsr_data_all_mk[[row_index_select]]
                if(adopt_model_fsr == "LUS"){
                  rv_dataplot$fsr_xval <- fsr_data_selected$fsr_lus
                }else if(adopt_model_fsr == "Multi-sequence"){
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
              adopt_model_dsr <- dt_result_dsr$Model[row_index_select]
              if(adopt_model_dsr != "Not consider"){
                dsr_data_selected <- dsr_data_all_mk[[row_index_select]]
                if(adopt_model_dsr == "LUS"){
                  rv_dataplot$dsr_xval <- dsr_data_selected$dsr_lus
                }else if(adopt_model_dsr == "Multi-sequence"){
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
              adopt_model_m2sr <- dt_result_m2sr$Model[row_index_select]
              if(adopt_model_m2sr != "Not consider"){
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

      #### Show allele repeat corrections ####

      rv_dt_alcor <- reactiveVal(NULL)

      # Change the status of the 'Corrected allele' button
      observeEvent(input$tabset_est_par_result, {
        if(is.element(input$tabset_est_par_result, c("AE", "M2SR", "Gamma"))){
          disable("act_alcor")
        }else{
          enable("act_alcor")
        }
      })

      observeEvent(input$act_alcor, {
        bool_make_dt <- FALSE
        tabset_est_par_result <- input$tabset_est_par_result

        # BSR
        if(tabset_est_par_result == "BSR"){
          if(!is.null(dt_alcor_bsr)){
            row_index_select <- input$dt_result_bsr_rows_selected
            if(isTruthy(row_index_select)){
              adopt_model_bsr <- dt_result_bsr$Model[row_index_select]
              if(is.element(adopt_model_bsr, c("LUS", "Multi-sequence"))){
                bool_make_dt <- TRUE
                mk <- markers_par[row_index_select]
                setkey(dt_alcor_bsr, Marker)
                dt_alcor_bsr_mk <- dt_alcor_bsr[mk]
                alleles <- dt_alcor_bsr_mk[, Allele]
                alcors <- dt_alcor_bsr_mk[, Corrected]
              }else{
                showModal(modalDialog(title = "Information", "Information on corrected alleles are shown only when the adopted model is 'LUS' or 'Multi-sequence'.", easyClose = TRUE, footer = NULL))
              }
            }else{
              showModal(modalDialog(title = "Error", "Select a row!", easyClose = TRUE, footer = NULL))
            }
          }
        }

        # FSR
        else if(tabset_est_par_result == "fsr"){
          if(!is.null(dt_alcor_fsr)){
            row_index_select <- input$dt_result_fsr_rows_selected
            if(isTruthy(row_index_select)){
              adopt_model_fsr <- dt_result_fsr$Model[row_index_select]
              if(is.element(adopt_model_fsr, c("LUS", "Multi-sequence"))){
                bool_make_dt <- TRUE
                mk <- markers_par[row_index_select]
                setkey(dt_alcor_fsr, Marker)
                dt_alcor_fsr_mk <- dt_alcor_fsr[mk]
                alleles <- dt_alcor_fsr_mk[, Allele]
                alcors <- dt_alcor_fsr_mk[, Corrected]
              }else{
                showModal(modalDialog(title = "Information", "Information on corrected alleles are shown only when the adopted model is 'LUS' or 'Multi-sequence'.", easyClose = TRUE, footer = NULL))
              }
            }else{
              showModal(modalDialog(title = "Error", "Select a row!", easyClose = TRUE, footer = NULL))
            }
          }
        }

        # DSR
        else if(tabset_est_par_result == "dsr"){
          if(!is.null(dt_alcor_dsr)){
            row_index_select <- input$dt_result_dsr_rows_selected
            if(isTruthy(row_index_select)){
              adopt_model_dsr <- dt_result_dsr$Model[row_index_select]
              if(is.element(adopt_model_dsr, c("LUS", "Multi-sequence"))){
                bool_make_dt <- TRUE
                mk <- markers_par[row_index_select]
                setkey(dt_alcor_dsr, Marker)
                dt_alcor_dsr_mk <- dt_alcor_dsr[mk]
                alleles <- dt_alcor_dsr_mk[, Allele]
                alcors <- dt_alcor_dsr_mk[, Corrected]
              }else{
                showModal(modalDialog(title = "Information", "Information on corrected alleles are shown only when the adopted model is 'LUS' or 'Multi-sequence'.", easyClose = TRUE, footer = NULL))
              }
            }else{
              showModal(modalDialog(title = "Error", "Select a row!", easyClose = TRUE, footer = NULL))
            }
          }
        }

        if(bool_make_dt){
          dt_alcor <- data.table(Allele = alleles,
                                 Corrected = signif(alcors, 3))
          rv_dt_alcor(dt_alcor)

          showModal(modalDialog(
            title = "Corrected allele",
            dataTableOutput(session$ns("dt_alcor"))
          ))
        }else{
          rv_dt_alcor(NULL)
        }
      })

      output$dt_alcor <- renderDataTable(server = FALSE, {
        dt_alcor <- rv_dt_alcor()
        datatable(
          dt_alcor,
          selection = "none",
          options = list(dom = "t", iDisplayLength = nrow(dt_alcor), ordering = FALSE, autowidth = TRUE),
          rownames = FALSE
        )
      })

      #### Show AIC values and parameters ####

      # Define reactive values for displaying AIC values and parameters of each candidate model
      rv_dt_aic <- reactiveValues()
      rv_dt_aic$bsr <- NULL
      rv_dt_aic$fsr <- NULL
      rv_dt_aic$dsr <- NULL

      # Change the status of the 'AIC' button
      observeEvent(input$tabset_est_par_result, {
        if(is.element(input$tabset_est_par_result, c("AE", "M2SR", "Gamma"))){
          disable("act_aic_model")
        }else{
          enable("act_aic_model")
        }
      })

      # Show the pop-up window to show AIC values and parameters of candidate models
      observeEvent(input$act_aic_model, {
        tabset_est_par_result <- input$tabset_est_par_result

        # BSR
        if(tabset_est_par_result == "BSR"){
          if(!is.null(aics_bsr_all_mk)){
            row_index_select <- input$dt_result_bsr_rows_selected
            if(isTruthy(row_index_select)){
              aics <- aics_bsr_all_mk[[row_index_select]]
              model_bsr <- model_bsr_all_mk[row_index_select]
              adopt_model_bsr <- dt_result_bsr$Model[row_index_select]
              if(model_bsr == "Best"){
                models_cand <- names(aics)

                if(is.element("Allele", models_cand)){
                  bsr_slopes <- bsr_allele_median_slope_all_mk[row_index_select]
                  bsr_intercepts <- bsr_allele_median_intercept_all_mk[row_index_select]
                  bsr_sigmas <- bsr_allele_median_sigma_all_mk[row_index_select]
                  bsr_min_us <- NA
                  bsr_rhat_slopes <- bsr_allele_rhat_slope_all_mk[row_index_select]
                  bsr_rhat_intercepts <- bsr_allele_rhat_intercept_all_mk[row_index_select]
                  bsr_rhat_sigmas <- bsr_allele_rhat_sigma_all_mk[row_index_select]
                  bsr_rhat_lps <- bsr_allele_rhat_lp_all_mk[row_index_select]
                }

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

                rv_dt_aic$bsr <- data.table(Model = models_cand,
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
                  title = "AIC values and parameters of candidate models",
                  dataTableOutput(session$ns("dt_aic_bsr")),
                  br(),
                  h5(div(paste0("Adopted model: ", adopt_model_bsr), style = "color:#555555;font-weight:bold;")),
                  footer = tagList(
                    modalButton("Close")
                  ),
                  size = "l"
                ))
              }else if(adopt_model_bsr == "Not consider"){
                showModal(modalDialog(title = "Information", "The back stutter ratio is not considered in the selected marker.", easyClose = TRUE, footer = NULL))
              }else{
                showModal(modalDialog(title = "Information", "Information on AIC values are shown only when 'Best' is selected in the setting tab.", easyClose = TRUE, footer = NULL))
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
              model_fsr <- model_fsr_all_mk[row_index_select]
              adopt_model_fsr <- dt_result_fsr$Model[row_index_select]
              if(model_fsr == "Best"){
                models_cand <- names(aics)

                if(is.element("Allele", models_cand)){
                  fsr_slopes <- fsr_allele_median_slope_all_mk[row_index_select]
                  fsr_intercepts <- fsr_allele_median_intercept_all_mk[row_index_select]
                  fsr_sigmas <- fsr_allele_median_sigma_all_mk[row_index_select]
                  fsr_min_us <- NA
                  fsr_rhat_slopes <- fsr_allele_rhat_slope_all_mk[row_index_select]
                  fsr_rhat_intercepts <- fsr_allele_rhat_intercept_all_mk[row_index_select]
                  fsr_rhat_sigmas <- fsr_allele_rhat_sigma_all_mk[row_index_select]
                  fsr_rhat_lps <- fsr_allele_rhat_lp_all_mk[row_index_select]
                }

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

                rv_dt_aic$fsr <- data.table(Model = models_cand,
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
                  title = "AIC values and parameters of candidate models",
                  dataTableOutput(session$ns("dt_aic_fsr")),
                  br(),
                  h5(div(paste0("Adopted model: ", adopt_model_fsr), style = "color:#555555;font-weight:bold;")),
                  footer = tagList(
                    modalButton("Close")
                  ),
                  size = "l"
                ))
              }else if(adopt_model_fsr == "Not consider"){
                showModal(modalDialog(title = "Information", "The forward stutter ratio is not considered in the selected marker.", easyClose = TRUE, footer = NULL))
              }else{
                showModal(modalDialog(title = "Information", "Information on AIC values are shown only when 'Best' is selected in the setting tab.", easyClose = TRUE, footer = NULL))
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
              model_dsr <- model_dsr_all_mk[row_index_select]
              adopt_model_dsr <- dt_result_dsr$Model[row_index_select]
              if(model_dsr == "Best"){
                models_cand <- names(aics)

                if(is.element("Allele", models_cand)){
                  dsr_slopes <- dsr_allele_median_slope_all_mk[row_index_select]
                  dsr_intercepts <- dsr_allele_median_intercept_all_mk[row_index_select]
                  dsr_sigmas <- dsr_allele_median_sigma_all_mk[row_index_select]
                  dsr_min_us <- NA
                  dsr_rhat_slopes <- dsr_allele_rhat_slope_all_mk[row_index_select]
                  dsr_rhat_intercepts <- dsr_allele_rhat_intercept_all_mk[row_index_select]
                  dsr_rhat_sigmas <- dsr_allele_rhat_sigma_all_mk[row_index_select]
                  dsr_rhat_lps <- dsr_allele_rhat_lp_all_mk[row_index_select]
                }

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

                rv_dt_aic$dsr <- data.table(Model = models_cand,
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
                  title = "AIC values and parameters of candidate models",
                  dataTableOutput(session$ns("dt_aic_dsr")),
                  br(),
                  h5(div(paste0("Adopted model: ", adopt_model_dsr), style = "color:#555555;font-weight:bold;")),
                  footer = tagList(
                    modalButton("Close")
                  ),
                  size = "l"
                ))
              }else if(adopt_model_dsr == "Not consider"){
                showModal(modalDialog(title = "Information", "The double-back stutter ratio is not considered in the selected marker.", easyClose = TRUE, footer = NULL))
              }else{
                showModal(modalDialog(title = "Information", "Information on AIC values are shown only when 'Best' is selected in the setting tab.", easyClose = TRUE, footer = NULL))
              }
            }else{
              showModal(modalDialog(title = "Error", "Select a row!", easyClose = TRUE, footer = NULL))
            }
          }
        }
      })

      # Display information on AIC values and parameters of each back stutter model
      output$dt_aic_bsr <- renderDataTable(server = FALSE, {
        dt_aic_bsr <- rv_dt_aic$bsr
        datatable(
          dt_aic_bsr,
          selection = "none",
          options = list(dom = "t", iDisplayLength = nrow(dt_aic_bsr), ordering = FALSE, autowidth = TRUE),
          rownames = FALSE
        )
      })

      # Display information on AIC values and parameters of each forward stutter model
      output$dt_aic_fsr <- renderDataTable(server = FALSE, {
        dt_aic_fsr <- rv_dt_aic$fsr
        datatable(
          dt_aic_fsr,
          selection = "none",
          options = list(dom = "t", iDisplayLength = nrow(dt_aic_fsr), ordering = FALSE, autowidth = TRUE),
          rownames = FALSE
        )
      })

      # Display information on AIC values and parameters of each double-back stutter model
      output$dt_aic_dsr <- renderDataTable(server = FALSE, {
        dt_aic_dsr <- rv_dt_aic$dsr
        datatable(
          dt_aic_dsr,
          selection = "none",
          options = list(dom = "t", iDisplayLength = nrow(dt_aic_dsr), ordering = FALSE, autowidth = TRUE),
          rownames = FALSE
        )
      })

      #### Save parameters ####

      rv_par_names <- reactiveVal(NULL)

      observeEvent(input$act_save_par, {
        par_names <- list.files(paste0(path_pack, "/extdata/parameter"))
        rv_par_names(gsub(".{4}$", "", par_names))

        showModal(modalDialog(
          title = "Save parameters",
          textInput(session$ns("proj_par"), label = "Enter the project name"),
          footer = tagList(
            disabled(actionButton(session$ns("act_save_par_2"), "Save")),
            modalButton("Close")
          )
        ))
      })

      observeEvent(input$proj_par, {
        isolate(proj_par <- input$proj_par)
        if(nchar(proj_par) == 0){
          showFeedbackDanger(inputId = "proj_par", text = "Required")
          disable("act_save_par_2")
        }else if(is.element(proj_par, rv_par_names())){
          showFeedbackDanger(inputId = "proj_par", text = "This name is already registered.")
          disable("act_save_par_2")
        }else{
          hideFeedback("proj_par")
          enable("act_save_par_2")
        }
      })

      observeEvent(input$act_save_par_2, {
        isolate(proj_par <- input$proj_par)
        write.csv(dt_par, paste0(path_pack, "/extdata/parameter/", proj_par, ".csv"), row.names = FALSE)
        removeModal()
      })

      #### Result table ####

      # Display the result for AE
      output$dt_result_ae <- renderDataTable(server = FALSE, {
        datatable(
          dt_result_ae,
          selection = list(mode = "single", target = "row"),
          options = list(dom = "t", iDisplayLength = nrow(dt_result_ae), ordering = FALSE, autowidth = TRUE),
          rownames = FALSE
        )
      })

      # Display the result for BSR
      output$dt_result_bsr <- renderDataTable(server = FALSE, {
        datatable(
          dt_result_bsr,
          selection = list(mode = "single", target = "row"),
          options = list(dom = "t", iDisplayLength = nrow(dt_result_bsr), ordering = FALSE, autowidth = TRUE),
          rownames = FALSE
        )
      })

      # Display the result for FSR
      output$dt_result_fsr <- renderDataTable(server = FALSE, {
        datatable(
          dt_result_fsr,
          selection = list(mode = "single", target = "row"),
          options = list(dom = "t", iDisplayLength = nrow(dt_result_fsr), ordering = FALSE, autowidth = TRUE),
          rownames = FALSE
        )
      })

      # Display the result for DSR
      output$dt_result_dsr <- renderDataTable(server = FALSE, {
        datatable(
          dt_result_dsr,
          selection = list(mode = "single", target = "row"),
          options = list(dom = "t", iDisplayLength = nrow(dt_result_dsr), ordering = FALSE, autowidth = TRUE),
          rownames = FALSE
        )
      })

      # Display the result for M2SR
      output$dt_result_m2sr <- renderDataTable(server = FALSE, {
        datatable(
          dt_result_m2sr,
          selection = list(mode = "single", target = "row"),
          options = list(dom = "t", iDisplayLength = nrow(dt_result_m2sr), ordering = FALSE, autowidth = TRUE),
          rownames = FALSE
        )
      })

      # Display the result for the gamma parameter
      output$dt_result_gamma <- renderDataTable(server = FALSE, {
        datatable(
          dt_result_gamma,
          selection = list(mode = "single", target = "row"),
          options = list(dom = "t", iDisplayLength = nrow(dt_result_gamma), ordering = FALSE, autowidth = TRUE),
          rownames = FALSE
        )
      })
    }
  )
}
