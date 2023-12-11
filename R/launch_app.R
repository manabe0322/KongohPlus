#' KongohPlus
#'
#' @description Main window
#' @usage KongohPlus()
#' @export
KongohPlus <- function(){
  ver_soft <- "0.1.0"
  #ver_soft <- packageVersion("KongohPlus")
  path_pack <- "D:/RStudio_GitHub/KongohPlus/inst"
  #path_pack <- "C:/Users/manab/Documents/RStudio/KongohPlus/inst"
  #path_pack <- path.package("KongohPlus", quiet = FALSE)
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
    est_par_server("est_par", path_pack)
  }

  shinyApp(ui, server)
}
