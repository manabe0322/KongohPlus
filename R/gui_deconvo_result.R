#' deconvo_result_ui
#'
#' @description The function to create the ui module for results of deconvolution
deconvo_result_ui <- function(id){
  ns <- NS(id)

  tabPanel("Result",
           br(),
           fluidRow(
             column(3,
                    wellPanel(
                      h4("Estimated parameter")
                    )
             ),
             column(6,
                    wellPanel(
                      h4("Probabilistic genotyping")
                    )
             ),
             column(3,
                    wellPanel(
                      h4("Likelihood")
                    )
             )
           )
  )
}

#' deconvo_result_server
#'
#' @description The function to create the server module for results of deconvolution
deconvo_result_server <- function(id, rv_deconvo_result){
  moduleServer(
    id,
    function(input, output, session){
      results <- rv_deconvo_result
      if(!is.null(results)){
        dt_pg_list <- results$dt_pg_list
        dt_likelihood <- results$dt_likelihood
        dt_mpdeg <- results$dt_mpdeg
        tph <- results$tph
        rhat_lp <- results$rhat_lp
        traceplot <- results$traceplot
      }
    }
  )
}
