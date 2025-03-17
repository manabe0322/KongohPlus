#' kit_ui
#'
#' @description The function to create the ui module for the kit
kit_ui <- function(id){
  ns <- NS(id)

  tabPanel("1. Register kits",
           titlePanel("1. Register kits"),
           br(),
           sidebarLayout(
             sidebarPanel(width = 3,
                          h4("Kit name"),
                          textInput(ns("name_kit"), label = "Enter the kit name", value = ""),
                          br(),
                          h4("Add marker"),
                          textInput(ns("marker_add"), label = "Enter the marker name", value = ""),
                          actionButton(ns("act_add_marker"), label = "Add"),
                          br(),
                          br(),
                          h4("Delete marker"),
                          uiOutput(ns("ui_marker_delete")),
                          actionButton(ns("act_delete_marker"), label = "Delete"),
                          br(),
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
                                       "Repeat Unit Length" = character(0),
                                       "Chromosome" = character(0))

      observeEvent(input$act_add_marker, {
        marker_add <- input$marker_add
        if(marker_add != ""){
          dt_info_kit <- hot_to_r(input$dt_info_kit)
          if(!is.element(marker_add, dt_info_kit$Marker)){
            dt_tmp <- data.table("Marker" = marker_add,
                                 "Repeat Unit Length" = NA,
                                 "Chromosome" = "Autosomal")
            rv_kit$dt_info_kit <- rbind(dt_info_kit, dt_tmp)
          }
        }
      })

      output$ui_marker_delete <- renderUI({
        markers <- rv_kit$dt_info_kit$Marker
        selectInput(session$ns("marker_delete"), label = "Select a marker", choices = markers, selected = markers[1], multiple = FALSE, selectize = FALSE)
      })

      observeEvent(input$act_delete_marker, ignoreInit = TRUE, {
        marker_delete <- input$marker_delete
        markers <- rv_kit$dt_info_kit$Marker
        index_delete <- which(markers == marker_delete)
        if(length(index_delete) == 1){
          dt_info_kit <- hot_to_r(input$dt_info_kit)
          rv_kit$dt_info_kit <- dt_info_kit[-index_delete]
        }
      })

      output$dt_info_kit <- renderRHandsontable({
        rhandsontable(rv_kit$dt_info_kit, rowHeaders = FALSE, overflow = "visible") %>%
          hot_table(fillHandle = FALSE, stretchH = "all", rowHeights = 30) %>%
          hot_col(c("Marker"), type = "text", readOnly = TRUE) %>%
          hot_col(c("Repeat Unit Length"), type = "dropdown", source = as.character(c(2, 3, 4, 5, 6, 7))) %>%
          hot_col(c("Chromosome"), type = "dropdown", source = c("Autosomal", "Sex"))
      })

      observeEvent(input$act_register_kit, {
        name_kit <- input$name_kit
        dt_info_kit <- hot_to_r(input$dt_info_kit)
        col_integer <- c("Repeat Unit Length")
        dt_info_kit[, (col_integer) := lapply(.SD, as.integer), .SDcols = col_integer]
        index_auto <- which(dt_info_kit[["Chromosome"]] == "Autosomal")
        markers_auto <- dt_info_kit[["Marker"]][index_auto]
        rep_uni_len_auto <- as.numeric(dt_info_kit[["Repeat Unit Length"]][index_auto])

        tmp <- list.files(paste0(path_pack, "/extdata/kit"), pattern = "csv")
        kits_registerd <- str_replace(tmp, ".{4}$", "")
        if(name_kit == ""){
          showModal(modalDialog(title = "Error", "The kit name must be defined.", easyClose = TRUE, footer = NULL))
        }else if(is.element(name_kit, kits_registerd)){
          showModal(modalDialog(title = "Error", "The kit name is already registered.", easyClose = TRUE, footer = NULL))
        }else if(length(index_auto) == 0){
          showModal(modalDialog(title = "Error", "At least one autosomal marker must be defined.", easyClose = TRUE, footer = NULL))
        }else if(any(is.na(rep_uni_len_auto))){
          markers_error <- markers_auto[is.na(rep_uni_len_auto)]
          error_message <- "The 'Repeat Unit Length' of the following markers must be defined:<br>"
          for(i in 1:length(markers_error)){
            error_message <- paste0(error_message, "  - ", markers_error[i], "<br>")
          }
          showModal(modalDialog(title = "Error", HTML(error_message), easyClose = TRUE, footer = NULL))
        }else{
          n_mk <- nrow(dt_info_kit)
          dt_kit_right <- data.table("Allele" = rep("", n_mk),
                                     "Corrected Allele" = rep("", n_mk),
                                     "Size" = rep("", n_mk))
          dt_kit <- cbind(dt_info_kit, dt_kit_right)

          write.csv(dt_kit, paste0(path_pack, "/extdata/kit/", name_kit, ".csv"), row.names = FALSE)
          showModal(modalDialog(title = "Information", "The kit information was saved successfully.", easyClose = TRUE, footer = NULL))
        }
      })
    }
  )
}
