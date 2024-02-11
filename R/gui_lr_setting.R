#' lr_setting_ui
#'
#' @description The function to create the ui module for setting conditions of the likelihood ratio
lr_setting_ui <- function(id){
  ns <- NS(id)

  tabPanel("Likelihood ratio",
           br(),
           sidebarLayout(
             sidebarPanel(width = 3,
                          numericInput(ns("noc_lr"), label = h4("Number of contributors"), value = 1),
                          uiOutput(ns("hp")),
                          uiOutput(ns("hd")),
                          br(),
                          uiOutput(ns("anameth_lr")),
                          br(),
                          actionButton(ns("act_analyze_lr"), label = "Analyze", class = "btn btn-primary btn-lg")
             ),
             mainPanel(width = 9,
                       fluidRow(
                         column(4,
                                dataTableOutput(ns("dt_at_lr")),
                                br()
                         ),
                         column(8,
                                h4("Kit")
                         )
                       )
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

      rv_lr_setting <- reactiveValues()
      rv_lr_setting$refs <- NULL
      rv_lr_setting$dt_anameth <- NULL

      if(!is.null(dt_ref)){
        rv_lr_setting$refs <- unique(dt_ref[, SampleName])
      }

      observeEvent(input$noc_lr, {
        noc_lr <- input$noc_lr
        if(!is.integer(noc_lr) || noc_lr < 1 || noc_lr > 4){
          showFeedbackDanger(inputId = "noc_lr", text = "Only 1, 2, 3, and 4 is allowed.")
          disable("act_analyze_lr")
        }else{
          hideFeedback("noc_lr")
          enable("act_analyze_lr")
        }
      })

      output$hp <- renderUI({
        checkboxGroupInput(session$ns("hp"), label = h4("Prosecutor hypothesis"), choices = rv_lr_setting$refs)
      })

      output$hd <- renderUI({
        checkboxGroupInput(session$ns("hd"), label = h4("Defense hypothesis"), choices = rv_lr_setting$refs)
      })

      anameth <- list.files(paste0(path_pack, "/extdata/anameth"))

      output$anameth_lr <- renderUI({
        anameth_all <- list.files(paste0(path_pack, "/extdata/anameth"))
        anameth_all <- gsub(".{4}$", "", anameth_all)
        selectInput(session$ns("anameth_lr"), label = h4("Analysis method"), choices = anameth_all, selected = anameth_all[1])
      })

      observeEvent(input$anameth_lr, {
        # Load the selected analysis method
        anameth <- input$anameth_lr
        rv_lr_setting$dt_anameth <- dt_anameth <- fread(paste0(path_pack, "/extdata/anameth/", anameth, ".csv"))

        # Extract data per column
        folders <- dt_anameth[, Folder]
        names <- dt_anameth[, Name]

        rv_deconvo_setting$dt_at <- dt_at <- fread(paste0(path_pack, "/extdata/at/", names[folders == "at"], ".csv"))

        # Display analytical thresholds
        output$dt_at_lr <- renderDataTable(server = FALSE, {
          datatable(
            dt_at,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_at), ordering = FALSE, autowidth = TRUE),
            rownames = FALSE
          )
        })
      })
    }
  )
}
