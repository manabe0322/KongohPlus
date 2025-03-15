#' kit_ui
#'
#' @description The function to create the ui module for the kit
kit_ui <- function(id){
  ns <- NS(id)

  tabPanel("Kit",
           titlePanel("Kit"),
           br(),
           sidebarLayout(
             sidebarPanel(width = 3,
                          textInput(ns("name_kit"), label = h4("Kit name"), value = ""),
                          br(),
                          h4("Add marker"),
                          textInput(ns("marker_add"), label = "Enter the marker name", value = ""),
                          br(),
                          actionButton(ns("act_add_marker"), label = "Add"),
                          br(),
                          br(),
                          h4("Delete marker"),
                          uiOutput(ns("ui_marker_delete")),
                          br(),
                          actionButton(ns("act_delete_marker"), label = "Delete"),
                          br(),
                          br(),
                          actionButton(ns("act_register_kit"), label = "Register", class = "btn btn-primary btn-lg")
             ),
             mainPanel(width = 9,
                       rHandsontableOutput(ns("dt_info_kit"))
             )
           )
  )
}

#' kit_server
#'
#' @description The function to create the server module for the kit
kit_server <- function(id, path_pack){
  moduleServer(
    id,
    function(input, output, session){
      rv_kit <- reactiveValues()
      rv_kit$dt_info_kit <- data.table("Marker" = character(0),
                                       "Repeat unit length" = integer(0),
                                       "Chromosome" = character(0))

      observeEvent(input$act_add_marker, {
        marker_add <- input$marker_add
        if(marker_add != ""){
          dt_info_kit <- rv_kit$dt_info_kit
          dt_tmp <- data.table("Marker" = marker_add,
                               "Repeat unit length" = NA,
                               "Chromosome" = "Auto")
          rv_kit$dt_info_kit <- rbind(dt_info_kit, dt_tmp)
        }
      })

      output$ui_marker_delete <- renderUI({
        markers <- rv_kit$dt_info_kit$Marker
        selectInput(session$ns("marker_delete"), label = "", choices = markers, selected = markers[1], multiple = FALSE, selectize = FALSE)
      })

      observeEvent(input$act_delete_marker, ignoreInit = TRUE, {
        marker_delete <- input$marker_delete
        markers <- rv_kit$dt_info_kit$Marker
        index_delete <- which(markers == marker_delete)
        if(length(index_delete) == 1){
          rv_kit$dt_info_kit <- rv_kit$dt_info_kit[-index_delete]
        }
      })

      output$dt_info_kit <- renderRHandsontable({
        rhandsontable(rv_kit$dt_info_kit, rowHeaders = FALSE, overflow = "visible") %>%
          hot_table(fillHandle = FALSE, stretchH = "all", rowHeights = 30) %>%
          hot_col(c("Marker"), type = "text", readOnly = TRUE) %>%
          hot_col(c("Repeat unit length"), type = "numeric", format = "0") %>%
          hot_col(c("Chromosome"), type = "dropdown", source = c("Autosomal", "Sex"))
      })

      observeEvent(input$act_register_kit, {
        dt_info_kit <- rv_kit$dt_info_kit
      })
    }
  )
}
