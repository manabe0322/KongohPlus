#' deconvo_ui
#'
#' @description The function to create the ui module for deconvolution
deconvo_ui <- function(id){
  ns <- NS(id)

  tabPanel("Deconvolution",
           sidebarLayout(
             sidebarPanel(width = 3,
                          numericInput(ns("noc_deconvo"), label = h4("Number of contributors"), value = NULL),
                          br(),
                          actionButton(ns("act_analyze_deconvo"), label = "Analyze", class = "btn btn-primary btn-lg")
             ),
             mainPanel(width = 9,
                       uiOutput(ns("analysis_method_deconvo")),
                       br(),
                       fluidRow(
                         column(4,
                                dataTableOutput(ns("dt_at")),
                                br()
                         ),
                         column(8,
                                h4("Kit"),
                                textOutput(ns("kit")),
                                br(),
                                h4("Allele probability"),
                                textOutput(ns("alprob")),
                                br(),
                                h4("Theta"),
                                textOutput(ns("theta")),
                                br(),
                                h4("Parameters"),
                                actionButton(ns("act_ae"), label = "AE"),
                                actionButton(ns("act_bsr"), label = "BSR"),
                                actionButton(ns("act_fsr"), label = "FSR"),
                                actionButton(ns("act_dsr"), label = "DSR"),
                                actionButton(ns("act_m2sr"), label = "M2SR"),
                                br(),
                                br(),
                                h4("Conditions to remove unrealistic genotype combinations")
                         )
                       )
             )
           )
  )
}

#' deconvo_server
#'
#' @description The function to create the server module for deconvolution
deconvo_server <- function(id, path_pack){
  moduleServer(
    id,
    function(input, output, session){
      analysis_methods <- list.files(paste0(path_pack, "/extdata/analysis_method"))

      output$analysis_method_deconvo <- renderUI({
        selectInput(session$ns("analysis_method_deconvo"), label = h4("Analysis method"), choices = analysis_methods, selected = analysis_methods[1])
      })

      observeEvent(input$analysis_method_deconvo, {
        isolate(analysis_method <- input$analysis_method_deconvo)
        dt_analysis_method <- fread(paste0(path_pack, "/extdata/analysis_method/", analysis_method))
        names <- dt_analysis_method[, Name]
        values <- dt_analysis_method[, Value]

        index_at <- grep("AT_", names)
        dt_at <- data.table(Marker = gsub("AT_", "", names[index_at]),
                            AT = values[index_at])

        output$dt_at <- renderDataTable(server = FALSE, {
          datatable(
            dt_at,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_at), ordering = FALSE, autowidth = TRUE),
            rownames = FALSE
          )
        })

        output$kit <- renderText({values[names == "Kit"]})
        output$alprob <- renderText({values[names == "Allele probability"]})
        output$theta <- renderText({values[names == "Theta"]})
      })
    }
  )
}
