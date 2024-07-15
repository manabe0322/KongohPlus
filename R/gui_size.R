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
                          uiOutput(ns("select_marker_size")),
                          br(),
                          actionButton(ns("act_size_edit"), label = "Edit")
             ),
             mainPanel(width = 9,
                       dataTableOutput(ns("dt_info_size"))
             )
           )
  )
}

#' create_init_dt_info_size
#'
#' @description The function to create the initial data.table for information on the allele sizes
create_init_dt_info_size <- function(path_pack){
  kit_all <- list.files(paste0(path_pack, "/extdata/kit"))
  if(length(kit_all) != 0){
    dt_kit <- fread(paste0(path_pack, "/extdata/kit/", kit_all[1]))
    dt_kit <- change_data_type_kit(dt_kit)
    marker <- dt_kit$Marker[1]
    dt_info_size <- create_dt_info_size(dt_kit, marker)
  }else{
    dt_info_size <- NULL
  }
  return(dt_info_size)
}

#' create_dt_info_size
#'
#' @description The function to create the data.table for information on the allele sizes
create_dt_info_size <- function(dt_kit, marker){
  index_mk <- which(dt_kit$Marker == marker)
  dt_info_size <- data.table(Marker = dt_kit$Marker[index_mk],
                             Dye = dt_kit$Dye[index_mk],
                             Length_of_the_repeat_unit = dt_kit$Length_of_the_repeat_unit[index_mk],
                             Chr = dt_kit$Chr[index_mk],
                             Allele = dt_kit$Allele[index_mk],
                             Size = dt_kit$Size[index_mk])
  return(dt_info_size)
}

#' size_server
#'
#' @description The function to create the server module for the allele sizes
size_server <- function(id, path_pack){
  moduleServer(
    id,
    function(input, output, session){

      ######################################
      # Define the initial reactive values #
      ######################################

      rv_size <- reactiveValues()
      rv_size$dt_info_size <- create_init_dt_info_size(path_pack)
      rv_size$dt_kit <- NULL

      output$select_kit_size <- renderUI({
        kit_all <- list.files(paste0(path_pack, "/extdata/kit"))
        kit_all <- gsub(".{4}$", "", kit_all)
        selectInput(session$ns("select_kit_size"), label = h4("Kit"), choices = kit_all, selected = kit_all[1], selectize = FALSE)
      })

      observeEvent(input$select_kit_size, {
        dt_kit <- fread(paste0(path_pack, "/extdata/kit/", input$select_kit_size, ".csv"))
        rv_size$dt_kit <- dt_kit <- change_data_type_kit(dt_kit)
        markers <- unique(dt_kit$Marker)

        output$select_marker_size <- renderUI({
          selectInput(session$ns("select_marker_size"), label = h4("Marker"), choices = markers, selected = markers[1], selectize = FALSE)
        })
      })

      observeEvent(input$select_marker_size, {
        rv_size$dt_info_size <- create_dt_info_size(rv_size$dt_kit, input$select_marker_size)
      })

      output$dt_info_size <- renderDataTable({
        dt_info_size <- rv_size$dt_info_size
        if(!is.null(dt_info_size)){
          datatable(
            dt_info_size,
            colnames = c("Marker", "Dye", "Length of the repeat unit", "Chromosome", "Allele", "Size"),
            selection = list(mode = "single", target = "row"),
            options = list(dom = "t", iDisplayLength = nrow(dt_info_size), ordering = FALSE, autoWidth = TRUE),
            rownames = FALSE
          ) %>%
            formatRound(columns = c("Size"), digits = 2)
        }
      })
    }
  )
}
