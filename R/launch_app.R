#' KongohPlus
#'
#' @description Main window
#' @usage KongohPlus()
#' @export
KongohPlus <- function(bool_compile = TRUE){
  ver_soft <- packageVersion("KongohPlus")
  path_pack <- path.package("KongohPlus", quiet = FALSE)
  options(shiny.maxRequestSize = 500 * 1024^2)

  # Compile stan codes
  if(bool_compile){
    compile_stan_code(path_pack)
  }

  ui <- fluidPage(useShinyjs(),
                  useShinyFeedback(),
                  theme = shinytheme("cerulean"),
                  tags$style(type = "text/css", "body{padding-top: 70px;}"),

                  navbarPage(title = paste0("KongohPlus ver. ", ver_soft),
                             id = "navbar",
                             position = c("fixed-top"),

                             load_ui("load"),

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
    rv_file <- load_server("load", session)
    observe({
      deconvo_server("deconvo", path_pack, rv_file)
      lr_server("lr", path_pack, rv_file)
    })

    est_par_server("est_par", path_pack)
  }

  shinyApp(ui, server)
}
