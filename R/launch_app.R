#' KongohPlus
#'
#' @description Main window
#' @usage KongohPlus()
#' @export
KongohPlus <- function(){
  ver_soft <- packageVersion("KongohPlus")
  path_pack <- path.package("KongohPlus", quiet = FALSE)
  options(shiny.maxRequestSize = 500 * 1024^2)

  ui <- fluidPage(useShinyjs(),
                  theme = shinytheme("cerulean"),
                  navbarPage(title = paste0("KongohPlus ver. ", ver_soft),
                             id = "navbar",
                             position = c("fixed-top"),
                             tags$style(type = "text/css", "body{padding-top: 70px;}"),

                             tabPanel("Load"),

                             tabPanel("Deconvolution"),

                             tabPanel("Likelihood ratio"),

                             tabPanel("Project"),

                             navbarMenu("Tools",
                                        tabPanel("Example files"),
                                        tabPanel("Support manual interpretation"),
                                        tabPanel("Estimate parameters"),
                                        tabPanel("Validation mode")
                             ),

                             navbarMenu("Setting",
                                        tabPanel("Analysis method"),
                                        tabPanel("Allele repeat correction"),
                                        tabPanel("Typing kit")
                             ),

                             tabPanel("Manual")
                  )
  )

  server <- function(input, output, session){

  }

  shinyApp(ui, server)
}
