#' size_ui
#'
#' @description The function to create the ui module for the allele sizes
size_ui <- function(id){
  ns <- NS(id)

  tabPanel("Allele size",
           titlePanel("Allele size"),
           br(),
           sidebarLayout(
             sidebarPanel(width = 3,
                          uiOutput(ns("select_kit_size")),
                          br(),
                          actionButton(ns("act_size_edit"), label = "Edit")
             ),
             mainPanel(width = 9,
                       dataTableOutput(ns("dt_info_size"))
             )
           )
  )
}

#' size_server
#'
#' @description The function to create the server module for the allele sizes
size_server <- function(id, path_pack){
  moduleServer(
    id,
    function(input, output, session){
      output$select_kit_size <- renderUI({
        kit_all <- list.files(paste0(path_pack, "/extdata/kit"))
        kit_all <- gsub(".{4}$", "", kit_all)
        selectInput(session$ns("select_kit_size"), label = h4("Kit"), choices = kit_all, selected = kit_all[1], selectize = FALSE)
      })
    }
  )
}
