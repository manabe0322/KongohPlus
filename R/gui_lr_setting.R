#' lr_setting_ui
#'
#' @description The function to create the ui module for setting conditions of the likelihood ratio
lr_setting_ui <- function(id){
  ns <- NS(id)

  tabPanel("Likelihood ratio",
           br(),
           sidebarLayout(
             sidebarPanel(width = 3,
                          numericInput(ns("noc_lr"), label = h4("Number of contributors"), value = NULL),
                          uiOutput(ns("hp")),
                          uiOutput(ns("hd")),
                          br(),
                          actionButton(ns("act_analyze_lr"), label = "Analyze", class = "btn btn-primary btn-lg")
             ),
             mainPanel(width = 9,
                       uiOutput(ns("analysis_method_lr"))
             )
           )
  )
}

#' lr_setting_server
#'
#' @description The function to create the server module for setting conditions of the likelihood ratio
lr_setting_server <- function(id, path_pack, rv_file){
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

      analysis_methods <- list.files(paste0(path_pack, "/extdata/analysis_method"))

      output$analysis_method_lr <- renderUI({
        selectInput(session$ns("analysis_method_lr"), label = h4("Analysis method"), choices = analysis_methods, selected = analysis_methods[1])
      })
    }
  )
}
