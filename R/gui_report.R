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
report_server <- function(id, path_pack, rv_deconvo){
  moduleServer(
    id,
    function(input, output, session){
      data_deconvo <- rv_deconvo$data_deconvo
      dt_csp <- data_deconvo$dt_csp

      output$download_report <- downloadHandler(
        filename = "report.html",
        content = function(file){
          report_tmp <- file.path(tempdir(), "report.Rmd")
          file.copy(paste0(path_pack, "/rmd/report.Rmd"), report_tmp, overwrite = TRUE)
          params <- list(dt_csp = dt_csp)
          render(report_tmp, output_file = file, params = params, envir = new.env(parent = globalenv()))
        }
      )
    }
  )
}
