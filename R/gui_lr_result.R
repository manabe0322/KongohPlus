#' lr_result_ui
#'
#' @description The function to create the ui module for results of the likelihood ratio
lr_result_ui <- function(id){
  ns <- NS(id)

  tabPanel("Result",
           br(),
           fluidRow(
             column(6,
                    wellPanel(
                      h4("Probabilistic genotyping")
                    )
             ),
             column(6,
                    wellPanel(
                      h4("Likelihood ratio")
                    )
             )
           )
  )
}

#' lr_result_server
#'
#' @description The function to create the server module for results of the likelihood ratio
lr_result_server <- function(id){

}
