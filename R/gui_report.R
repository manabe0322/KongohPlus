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

      observe({
        req(rv_deconvo)
        data_deconvo <- rv_deconvo$data_deconvo

        if(!is.null(data_deconvo)){
          rv_report$dt_csp <- data_deconvo$dt_csp

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
              dt_csp = rv_report$dt_csp
            ),
            envir = new.env(parent = globalenv())
          )
        }
      )
    }
  )
}
