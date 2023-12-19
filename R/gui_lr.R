#' lr_ui
#'
#' @description The function to create the ui module for the likelihood ratio
lr_ui <- function(id){
  ns <- NS(id)

  tabPanel("Likelihood ratio",
           sidebarLayout(
             sidebarPanel(width = 3,
                          numericInput(ns("noc_lr"), label = h4("Number of contributors"), value = NULL),
                          uiOutput(ns("hp")),
                          uiOutput(ns("hd"))
             ),
             mainPanel(width = 9,
                       selectInput(ns("analysis_method_lr"), label = h4("Analysis method"), choices = LETTERS[1:3])
             )
           )
  )
}

#' lr_server
#'
#' @description The function to create the server module for the likelihood ratio
lr_server <- function(id, rv_file){
  moduleServer(
    id,
    function(input, output, session){
      dt_ref <- rv_file$dt_ref

      rv_setting_lr <- reactiveValues()
      rv_setting_lr$refs <- NULL

      if(!is.null(dt_ref)){
        rv_setting_lr$refs <- unique(dt_ref[, SampleName])
      }

      output$hp <- renderUI({
        checkboxGroupInput(session$ns("hp"), label = h4("Prosecutor hypothesis"), choices = rv_setting_lr$refs)
      })

      output$hd <- renderUI({
        checkboxGroupInput(session$ns("hd"), label = h4("Defense hypothesis"), choices = rv_setting_lr$refs)
      })
    }
  )
}
