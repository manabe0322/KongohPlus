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

  stan_models <- list.files("D:/RStudio_GitHub/KongohPlus/data")
  #stan_models <- list.files(paste0(path_pack, "/data"))
  if(!is.element("model_ae.rds", stan_models)){
    model_ae <- cmdstan_model(paste0(path_pack, "/stan/model_ae.stan"))
    saveRDS(model_ae, "D:/RStudio_GitHub/KongohPlus/data/model_ae.rds")
    #saveRDS(model_ae, paste0(path_pack, "/data/model_ae.rds"))
  }

  ui <- fluidPage(useShinyjs(),
                  theme = shinytheme("cerulean"),
                  tags$style(type = "text/css", "body{padding-top: 70px;}"),

                  navbarPage(title = paste0("KongohPlus ver. ", ver_soft),
                             id = "navbar",
                             position = c("fixed-top"),


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
