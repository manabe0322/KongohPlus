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

      numericColumn <- function(len, col, ...){
        inputs <- character(len)
        for(i in seq_len(len)){
          inputs[i] <- as.character(
            numericInput(paste0("numeric_", col, "_", i), label = NULL, ...)
          )
        }
        return(inputs)
      }

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
          title = "Create a new analysis method",
          textInput(session$ns("anameth_new_input_name"), label = "Enter a name of the new analysis method", value = NULL),
          tabsetPanel(id = session$ns("tabset_new_anameth"),
                      type = "tabs",
                      tabPanel("Parameter & Kit",
                               br(),
                               uiOutput(session$ns("anameth_new_param")),
                               br(),
                               h4("Kit"),
                               textOutput(session$ns("anameth_new_kit"))
                      ),
                      tabPanel("AT",
                               br(),
                               dataTableOutput(session$ns("anameth_new_at"))
                      ),
                      tabPanel("Allele probability",
                               br()
                      ),
                      tabPanel("Allele repeat correction",
                               br()
                      ),
                      tabPanel("Drop-in",
                               br()
                      )
          ),
          footer = tagList(
            #actionButton(session$ns("act_new_anameth_save"), "Save"),
            modalButton("Cancel")
          ),
          size = "l"
        ))
      }, ignoreInit = TRUE)

      output$anameth_new_param <- renderUI({
        param_all <- list.files(paste0(path_pack, "/extdata/param"))
        param_all <- gsub(".{4}$", "", param_all)
        selectInput(session$ns("anameth_new_param"), label = h4("Parameter data"), choices = param_all)
      })

      observeEvent(input$anameth_new_param, {
        anameth_new_param <- input$anameth_new_param
        if(is.null(anameth_new_param)){
          output$anameth_new_kit <- renderText({"NA"})
        }else{
          # Set the kit name
          anameth_new_param <- input$anameth_new_param
          dt_par_select <- fread(paste0(path_pack, "/extdata/param/", anameth_new_param, ".csv"))
          kit_par <- dt_par_select$Kit[1]
          output$anameth_new_kit <- renderText({paste0(kit_par)})

          dt_kit_select <- fread(paste0(path_pack, "/extdata/kit/", kit_par, ".csv"))
          markers <- unique(dt_kit_select$Marker[dt_kit_select$Chr == "Auto"])
          n_mk <- length(markers)

          # Create a data.table for the AT
          output$anameth_new_at <- renderDataTable(server = FALSE, {
            dt_anameth_new_at <- data.table(Marker = markers)
            dt_anameth_new_at$AT <- numericColumn(n_mk, 2, value = 100, width = "75%")
            datatable(
              dt_anameth_new_at,
              selection = "none",
              options = list(dom = "t", iDisplayLength = nrow(dt_anameth_new_at), ordering = FALSE, autoWidth = TRUE),
              rownames = FALSE,
              escape = FALSE
            )
          })
        }
      })





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
