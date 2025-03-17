#' seq_data_ui
#'
#' @description The function to create the ui module for sequence data
seq_data_ui <- function(id){
  ns <- NS(id)

  tabPanel("Sequence data",
           titlePanel("Sequence data"),
           br(),
           sidebarLayout(
             sidebarPanel(width = 3,
                          fileInput(ns("file_seq"), label = h4("Load data"), accept = ".csv", buttonLabel = "Load"),
                          uiOutput(ns("select_seq")),
                          br(),
                          uiOutput(ns("select_marker_seq"))
             ),
             mainPanel(width = 9,
                       dataTableOutput(ns("dt_info_seq"))
             )
           )
  )
}

#' create_init_dt_info_seq
#'
#' @description The function to create the initial data.table for sequence data
create_init_dt_info_seq <- function(path_pack){
  seq_all <- list.files(paste0(path_pack, "/extdata/seq"))
  if(length(seq_all) != 0){
    dt_seq <- fread(paste0(path_pack, "/extdata/seq/", seq_all[1]))
    marker <- dt_seq$Marker[1]
    dt_info_seq <- create_dt_info_seq(dt_seq, marker)
  }else{
    dt_info_seq <- NULL
  }
  return(dt_info_seq)
}

#' create_dt_seq
#'
#' @description The function to create the data.table for sequence data
create_dt_info_seq <- function(dt_seq, marker){
  index_mk <- which(dt_seq$Marker == marker)
  dt_info_seq <- data.table(Allele = dt_seq$Allele[index_mk],
                            Repeat_Region = dt_seq[["Repeat Region"]][index_mk],
                            Flanking_Region = dt_seq[["Flanking Region"]][index_mk],
                            Count = dt_seq$Count[index_mk])

  # Change data type
  col_char <- c("Allele", "Repeat_Region", "Flanking_Region")
  dt_info_seq[, (col_char) := lapply(.SD, as.character), .SDcols = col_char]
  col_int <- c("Count")
  dt_info_seq[, (col_int) := lapply(.SD, as.integer), .SDcols = col_int]

  return(dt_info_seq)
}

#' seq_data_server
#'
#' @description The function to create the server module for sequence data
seq_data_server <- function(id, path_pack){
  moduleServer(
    id,
    function(input, output, session){
      rv_seq <- reactiveValues()
      rv_seq$seq_all <- gsub(".{4}$", "", list.files(paste0(path_pack, "/extdata/seq")))
      rv_seq$dt_info_seq <- create_init_dt_info_seq(path_pack)

      # Load the crime stain profile
      observeEvent(input$file_seq, {
        file_input <- input$file_seq
        if(!is.null(file_input)){
          dt_seq <- fread(file_input$datapath)
          fn_dt_seq <- input$file_seq$name
          write.csv(dt_seq, paste0(path_pack, "/extdata/seq/", fn_dt_seq), row.names = FALSE)

          marker <- dt_seq$Marker[1]

          seq_before <- rv_seq$seq_all
          seq_load <- gsub(".{4}$", "", fn_dt_seq)
          seq_after <- c(seq_before, seq_load)
          rv_seq$seq_all <- seq_after
          updateSelectInput(inputId = "select_seq", choices = seq_after, selected = seq_load)

          rv_seq$dt_info_seq <- create_dt_info_seq(dt_seq, marker)
        }
      })

      output$select_seq <- renderUI({
        seq_all <- rv_seq$seq_all
        selectInput(session$ns("select_seq"), label = h4("Select data"), choices = seq_all, selected = seq_all[1], selectize = FALSE)
      })

      observeEvent(input$select_seq, {
        dt_seq <- fread(paste0(path_pack, "/extdata/seq/", input$select_seq, ".csv"))
        markers_seq <- unique(dt_seq$Marker)

        output$select_marker_seq <- renderUI({
          selectInput(session$ns("select_marker_seq"), label = h4("Marker"), choices = markers_seq, selected = markers_seq[1], selectize = FALSE)
        })
      })

      observeEvent(input$select_marker_seq, {
        dt_seq <- fread(paste0(path_pack, "/extdata/seq/", input$select_seq, ".csv"))
        rv_seq$dt_info_seq <- create_dt_info_seq(dt_seq, input$select_marker_seq)
      })

      output$dt_info_seq <- renderDataTable({
        dt_info_seq <- rv_seq$dt_info_seq
        if(!is.null(dt_info_seq)){
          datatable(
            dt_info_seq,
            colnames = c("Allele", "Repeat Region", "Flanking Region", "Count"),
            selection = "none",
            options = list(iDisplayLength = 10, ordering = FALSE, autoWidth = TRUE,
                           columnDefs = list(
                             list(targets = 0, width = "60px"),
                             list(targets = 1, width = "400px"),
                             list(targets = 2, width = "200px"),
                             list(targets = 3, width = "60px")
                           )
                           ),
            rownames = FALSE
          )
        }
      })
    }
  )
}
