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
report_server <- function(id){
  moduleServer(
    id,
    function(input, output, session){
      #output$download_report <- downloadHandler(
      #  filename = "report.html",
      #  content = function(file){
      #    report_tmp <- file.path(tempdir(), "report.Rmd")
      #    file.copy(paste0(path_pack, "/rmd/report.Rmd"), report_tmp, overwrite = TRUE)
      #    params <- list(dt_csp = dt_csp, hnc = hnc, markers_use = markers_use, dt_mpdeg = dt_mpdeg, dt_mcmc = dt_mcmc, dt_pg_list = dt_pg_list, dt_log10like = dt_log10like)
      #    render(report_tmp, output_file = file, params = params, envir = new.env(parent = globalenv()))
      #  }
      #)
    }
  )
}
