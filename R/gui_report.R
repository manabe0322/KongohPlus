#' report_ui
#'
#' @description The function to create the ui module for generating a report
report_ui <- function(id){
  ns <- NS(id)

  tabPanel("Report",
           downloadButton(ns("download_report"), "Report")
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
      rv_report$dt_af <- NULL
      rv_report$fn_af <- NULL
      rv_report$fn_anameth <- NULL
      rv_report$data_hnc_list <- NULL
      rv_report$dt_ignore_mk <- NULL
      rv_report$hnc_all <- NULL
      rv_report$theta <- NULL
      rv_report$dt_lr_calc_record <- NULL
      rv_report$lr_list <- NULL

      observe({
        req(rv_deconvo)
        data_deconvo <- rv_deconvo$data_deconvo

        if(!is.null(data_deconvo)){
          rv_report$dt_csp <- data_deconvo$dt_csp
          rv_report$fn_csp <- data_deconvo$fn_csp
          rv_report$dt_af <- data_deconvo$dt_af
          rv_report$fn_af <- data_deconvo$fn_af
          rv_report$fn_anameth <- data_deconvo$fn_anameth
          rv_report$data_hnc_list <- data_deconvo$data_hnc_list
          rv_report$dt_ignore_mk <- data_deconvo$dt_ignore_mk
          rv_report$hnc_all <- data_deconvo$hnc_all
          rv_report$theta <- data_deconvo$theta
        }
      })

      observe({
        req(rv_lr)

        rv_report$dt_lr_calc_record <- rv_lr$dt_lr_calc_record
        rv_report$lr_list <- rv_lr$lr_list
      })

      output$download_report <- downloadHandler(
        filename = function() {
          paste0("report_", Sys.Date(), ".pdf")
        },
        content = function(file){
          report_tmp <- file.path(tempdir(), "report.Rmd")
          file.copy(paste0(path_pack, "/rmd/report.Rmd"), report_tmp, overwrite = TRUE)
          render(
            input = report_tmp,
            output_file = file,
            params = list(
              dt_csp = rv_report$dt_csp,
              fn_csp = rv_report$fn_csp,
              dt_af = rv_report$dt_af,
              fn_af = rv_report$fn_af,
              fn_anameth = rv_report$fn_anameth,
              data_hnc_list = rv_report$data_hnc_list,
              dt_ignore_mk = rv_report$dt_ignore_mk,
              hnc_all = rv_report$hnc_all,
              theta = rv_report$theta,
              dt_lr_calc_record = rv_report$dt_lr_calc_record,
              lr_list = rv_report$lr_list
            ),
            envir = new.env(parent = globalenv())
          )
        }
      )
    }
  )
}
