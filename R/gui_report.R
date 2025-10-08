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
      rv_report$data_hnc_list <- NULL
      rv_report$hnc_all <- NULL
      rv_report$dt_lr_calc_record <- NULL

      observe({
        req(rv_deconvo)
        data_deconvo <- rv_deconvo$data_deconvo

        if(!is.null(data_deconvo)){
          rv_report$dt_csp <- data_deconvo$dt_csp
          rv_report$fn_csp <- data_deconvo$fn_csp
          rv_report$dt_af <- data_deconvo$dt_af
          rv_report$fn_af <- data_deconvo$fn_af
          rv_report$data_hnc_list <- data_deconvo$data_hnc_list
          rv_report$hnc_all <- data_deconvo$hnc_all
        }
      })

      observe({
        req(rv_lr)

        if(!is.null(data_lr)){
          rv_report$dt_lr_calc_record <- rv_lr$dt_lr_calc_record
        }
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
              data_hnc_list = rv_report$data_hnc_list,
              hnc_all = rv_report$hnc_all,
              dt_lr_calc_record = rv_report$dt_lr_calc_record
            ),
            envir = new.env(parent = globalenv())
          )
        }
      )
    }
  )
}
