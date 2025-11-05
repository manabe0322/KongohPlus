#' report_ui
#'
#' @description The function to create the ui module for generating a report
report_ui <- function(id){
  ns <- NS(id)

  tabPanel("Report",
           useWaiter(),
           textInput(ns("case_name"), label = h4("Case name")),
           textInput(ns("user_name"), label = h4("User name")),
           disabled(downloadButton(ns("download_report"), "Report"))
  )
}

#' report_server
#'
#' @description The function to create the server module for generating a report
report_server <- function(id, path_pack, rv_deconvo, rv_lr){
  moduleServer(
    id,
    function(input, output, session){
      rv_report <- reactiveValues()
      rv_report$dt_csp <- NULL
      rv_report$fn_csp <- NULL
      rv_report$fn_af <- NULL
      rv_report$dt_anameth <- NULL
      rv_report$fn_anameth <- NULL
      rv_report$markers_use <- NULL
      rv_report$data_hnc_list <- NULL
      rv_report$dt_ignore_mk <- NULL
      rv_report$hnc_all <- NULL
      rv_report$message_sr <- NULL
      rv_report$min_size <- NULL
      rv_report$max_size <- NULL
      rv_report$min_temp <- 1
      rv_report$max_temp <- 32000
      rv_report$min_rep_ae <- 0.2
      rv_report$max_rep_ae <- 5
      rv_report$min_deg <- 0.0001
      rv_report$stop_iter <- 500
      rv_report$max_iter <- 5000
      rv_report$num_converged <- 3
      rv_report$max_mle_runs <- 10
      rv_report$theta <- NULL
      rv_report$dt_ref <- NULL
      rv_report$fn_ref <- NULL
      rv_report$dt_lr_calc_record <- NULL
      rv_report$lr_list <- NULL
      rv_report$gt_prob_sum <- 0.99

      observe({
        req(rv_deconvo)
        data_deconvo <- rv_deconvo$data_deconvo

        if(!is.null(data_deconvo)){
          rv_report$dt_csp <- data_deconvo$dt_csp
          rv_report$fn_csp <- data_deconvo$fn_csp
          rv_report$fn_af <- data_deconvo$fn_af
          rv_report$dt_anameth <- data_deconvo$dt_anameth
          rv_report$fn_anameth <- data_deconvo$fn_anameth
          rv_report$markers_use <- data_deconvo$markers_use
          rv_report$data_hnc_list <- data_deconvo$data_hnc_list
          rv_report$dt_ignore_mk <- data_deconvo$dt_ignore_mk
          rv_report$hnc_all <- data_deconvo$hnc_all
          rv_report$message_sr <- data_deconvo$message_sr
          rv_report$min_size <- data_deconvo$min_size
          rv_report$max_size <- data_deconvo$max_size
          rv_report$min_temp <- data_deconvo$min_temp
          rv_report$max_temp <- data_deconvo$max_temp
          rv_report$min_rep_ae <- data_deconvo$min_rep_ae
          rv_report$max_rep_ae <- data_deconvo$max_rep_ae
          rv_report$min_deg <- data_deconvo$min_deg
          rv_report$stop_iter <- data_deconvo$stop_iter
          rv_report$max_iter <- data_deconvo$max_iter
          rv_report$num_converged <- data_deconvo$num_converged
          rv_report$max_mle_runs <- data_deconvo$max_mle_runs
          rv_report$theta <- data_deconvo$theta

          enable("download_report")
        }
      })

      observe({
        req(rv_lr)
        rv_report$dt_ref <- rv_lr$dt_ref
        rv_report$fn_ref <- rv_lr$fn_ref
        rv_report$dt_lr_calc_record <- rv_lr$dt_lr_calc_record
        rv_report$lr_list <- rv_lr$lr_list
      })

      output$download_report <- downloadHandler(
        filename = function() {
          paste0("report_", Sys.Date(), ".pdf")
        },
        content = function(file){
          waiter_show(html = spin_3k(), color = "white")

          report_tmp <- file.path(tempdir(), "report.Rmd")
          file.copy(paste0(path_pack, "/rmd/report.Rmd"), report_tmp, overwrite = TRUE)
          render(
            input = report_tmp,
            output_file = file,
            params = list(
              case_name = input$case_name,
              user_name = input$user_name,
              dt_csp = rv_report$dt_csp,
              fn_csp = rv_report$fn_csp,
              fn_af = rv_report$fn_af,
              dt_anameth = rv_report$dt_anameth,
              fn_anameth = rv_report$fn_anameth,
              markers_use = rv_report$markers_use,
              data_hnc_list = rv_report$data_hnc_list,
              dt_ignore_mk = rv_report$dt_ignore_mk,
              hnc_all = rv_report$hnc_all,
              message_sr = rv_report$message_sr,
              min_size = rv_report$min_size,
              max_size = rv_report$max_size,
              min_temp = rv_report$min_temp,
              max_temp = rv_report$max_temp,
              min_rep_ae = rv_report$min_rep_ae,
              max_rep_ae = rv_report$max_rep_ae,
              min_deg = rv_report$min_deg,
              stop_iter = rv_report$stop_iter,
              max_iter = rv_report$max_iter,
              num_converged = rv_report$num_converged,
              max_mle_runs = rv_report$max_mle_runs,
              theta = rv_report$theta,
              dt_ref = rv_report$dt_ref,
              fn_ref = rv_report$fn_ref,
              dt_lr_calc_record = rv_report$dt_lr_calc_record,
              lr_list = rv_report$lr_list,
              gt_prob_sum = rv_report$gt_prob_sum
            ),
            envir = new.env(parent = globalenv())
          )

          waiter_hide()
        }
      )
    }
  )
}
