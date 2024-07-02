#' anameth_ui
#'
#' @description The function to create the ui module for the analysis method
anameth_ui <- function(id){
  ns <- NS(id)

  tabPanel("Analysis method",
           titlePanel("Analysis method"),
           br(),
           sidebarLayout(
             sidebarPanel(width = 2,
                          actionButton(ns("act_anameth_new"), "New"),
                          br(),
                          br(),
                          actionButton(ns("act_anameth_edit"), "Edit"),
                          br(),
                          br(),
                          actionButton(ns("act_anameth_delete"), "Delete"),
                          br(),
                          br(),
                          actionButton(ns("act_anameth_copy"), "Copy")
             ),
             mainPanel(width = 10,
                       tags$head(
                         tags$style(HTML("table {table-layout: fixed;")),
                         tags$style(".datatables .display {margin-left: 0;}")
                       ),
                       dataTableOutput(ns("dt_info_anameth"))
             )
           )
  )
}

#' anameth_server
#'
#' @description The function to create the server module for the analysis method
anameth_server <- function(id, path_pack){
  moduleServer(
    id,
    function(input, output, session){

      ######################################
      # Define the initial reactive values #
      ######################################

      rv_anameth <- reactiveValues()
      rv_anameth$dt_info_anameth <- create_dt_info_anameth(path_pack)

      ##############################
      # Make a new analysis method #
      ##############################

      observeEvent(input$act_anameth_new, {
        showModal(modalDialog(
          title = "Make a new analysis method",
          textInput(session$ns("input_name_new_anameth"), label = "Enter a name of the new analysis method", value = NULL),
          tabsetPanel(id = session$ns("tabset_new_anameth"),
                      type = "tabs",
                      tabPanel("Kit & Parameter"),
                      tabPanel("AT"),
                      tabPanel("Allele probability"),
                      tabPanel("Drop-in")
          ),
          footer = tagList(
            actionButton(session$ns("act_new_anameth_save"), "Save"),
            modalButton("Cancel")
          ),
          size = "l"
        ))
      }, ignoreInit = TRUE)

      output$dt_info_anameth <- renderDataTable({
        dt_info_anameth <- rv_anameth$dt_info_anameth
        datatable(
          dt_info_anameth,
          colnames = c("Analysis method", "Kit", "Parameter file", "AT", "Allele probability", "Drop-in"),
          selection = list(mode = "single", target = "row"),
          extensions = c("FixedColumns", "FixedHeader"),
          options = list(dom = "t", iDisplayLength = nrow(dt_info_anameth), ordering = FALSE, autoWidth = TRUE, scrollX = TRUE, scrollY = "650px",
                         fixedHeader = TRUE, fixedColumns = list(leftColumns = 1, rightColumns = 0),
                         columnDefs = list(list(targets = c(0, 2), width = "300px"),
                                           list(targets = c(1, 3, 4, 5), width = "200px"))),
          rownames = FALSE
        )
      })
    }
  )
}
