#' KongohPlus
#'
#' @description Main window
#' @usage KongohPlus()
#' @export
KongohPlus <- function(bool_compile = TRUE){
  ver_soft <- packageVersion("KongohPlus")
  path_pack <- path.package("KongohPlus", quiet = FALSE)
  options(shiny.maxRequestSize = 500 * 1024^2)

  # Install and activate pandoc
  pandoc_install()
  pandoc_activate()

  ui <- fluidPage(theme = shinytheme("cerulean"),
                  useShinyjs(),
                  useShinyFeedback(),
                  tags$style(type = "text/css", "body{padding-top: 70px;}"),

                  navbarPage(title = paste0("KongohPlus ver. ", ver_soft),
                             id = "navbar",
                             position = c("fixed-top"),

                             file_ui("file"),

                             deconvo_ui("deconvo"),

                             lr_ui("lr"),

                             tabPanel("Project"),

                             navbarMenu("Tools",
                                        tabPanel("Example files"),
                                        tabPanel("Support manual interpretation"),
                                        est_par_ui("est_par"),
                                        tabPanel("Validation mode")
                             ),

                             navbarMenu("Setting",
                                        tabPanel("Analysis method"),
                                        tabPanel("Allele repeat correction"),
                                        tabPanel("Typing kit")
                             ),

                             tabPanel("Example files"),

                             tabPanel("Manual")
                  )
  )

  server <- function(input, output, session){
    rv_file <- file_server("file", session)
    rv_deconvo <- deconvo_server("deconvo", path_pack, rv_file)
    lr_server("lr", path_pack, rv_file, rv_deconvo)

    est_par_server("est_par", path_pack)
  }

  shinyApp(ui, server)
}
