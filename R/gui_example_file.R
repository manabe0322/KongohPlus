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
    downloadButton(ns("download_af"), "Allele frequencies"),
    br(),
    br(),

    tags$hr(),
    h4("References"),

    tags$p(
      tags$strong("Crime stain and reference profiles: "),
      "L.E. Alfonse, A.D. Garrett, D.S. Lun, K.R. Duffy, C.M. Grgicak, A large-scale dataset of single and mixed-source short tandem repeat profiles to inform human identification strategies: PROVEDIt, Forensic Sci. Int. Genet. 32 (2018) 62-70."
    ),

    tags$p(
      tags$strong("Allele frequencies: "),
      "C.R. Steffen, M.D. Coble, K.B. Gettings, P.M. Vallone, Corrigendum to 'U.S. Population Data for 29 Autosomal STR Loci' [Forensic Sci. Int. Genet. 7 (2013) e82-e83], Forensic Sci. Int. Genet. 31 (2017) e36-e40."
    )
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

      # Af: GlobalFiler, Caucasian

      output$download_af <- downloadHandler(
        filename = "example_af.csv",
        content = function(file){
          csvfile <- read.csv(paste0(path_pack, "/extdata/examples/example_af.csv"))
          csvfile[is.na(csvfile)] <- ""
          write.csv(csvfile, file, row.names = FALSE)
        }
      )
    }
  )
}
