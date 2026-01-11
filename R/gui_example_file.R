#' example_file_ui
#'
#' @description The function to create the ui module for the example files
#' @param id The ui ID
example_file_ui <- function(id){
  ns <- NS(id)

  tabPanel(
    "Example files",
    titlePanel("Example files"),
    br(),
    downloadButton(ns("download_csp"), "Crime stain profile"),
    br(),
    br(),
    downloadButton(ns("download_ref"), "Reference profile"),
    br(),
    br(),
    downloadButton(ns("download_af"), "Allele frequencies")
  )
}

#' example_file
#'
#' @description The function to create the server module for the example files
#' @param id The server ID
#' @param path_pack The path of the package KongohPlus
example_file_server <- function(id, path_pack){
  moduleServer(
    id,
    function(input, output, session){

      # --- Download the example file for crime stain profiles --- #

      # CSP: PROVEDIt, GlobalFiler, RD14-0003-31_32-1;1-M1a-0.25GF-Q1.2

      output$download_csp <- downloadHandler(
        filename = "example_csp.csv",
        content = function(file){
          csvfile <- read.csv(paste0(path_pack, "/extdata/examples/example_csp.csv"))
          csvfile[is.na(csvfile)] <- ""
          write.csv(csvfile, file, row.names = FALSE)
        }
      )

      # --- Download the example file for reference profiles --- #

      # Ref: PROVEDIt, GlobalFiler, 31 and 32

      output$download_ref <- downloadHandler(
        filename = "example_ref.csv",
        content = function(file){
          csvfile <- read.csv(paste0(path_pack, "/extdata/examples/example_ref.csv"))
          csvfile[is.na(csvfile)] <- ""
          write.csv(csvfile, file, row.names = FALSE)
        }
      )

      # --- Download the example file for allele frequencies --- #

      # Af: GlobalFiler, Japanese

      output$download_af <- downloadHandler(
        filename = "Allele-count_GF_Japanese.csv",
        content = function(file){
          csvfile <- read.csv(paste0(path_pack, "/extdata/examples/Allele-count_GF_Japanese.csv"))
          csvfile[is.na(csvfile)] <- ""
          write.csv(csvfile, file, row.names = FALSE)
        }
      )
    }
  )
}
