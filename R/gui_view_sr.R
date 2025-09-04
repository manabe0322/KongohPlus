#' view_sr_ui
#'
#' @description The function to create the ui module for viewing stutter ratios
view_sr_ui <- function(id){
  ns <- NS(id)

  tabPanel(
    "View stutter ratios",
    titlePanel("View stutter ratios"),
    br(),
    fluidRow(
      column(3,
             wellPanel(
               h4("Stutter data"),
               fluidRow(
                 column(3, br(), shinyFilesButton(ns("load_stutter_data"), label = "Load", title = "Select a file", multiple = FALSE)),
                 column(9, uiOutput(ns("fn_stutter_data")))
               ),
               uiOutput(ns("select_mk")),
               uiOutput(ns("select_st_type"))
             )
      ),
      column(4,
             wellPanel(
               dataTableOutput(ns("dt_average_sr_display"))
             )
      ),
      column(5,
             wellPanel(
               dataTableOutput(ns("dt_param_sr"))
             )
      )
    )
  )
}

#' view_sr_server
#'
#' @description The function to create the server module for viewing stutter ratios
view_sr_server <- function(id, path_pack){
  moduleServer(
    id,
    function(input, output, session){
      rv_view_sr <- reactiveValues()
      rv_view_sr$fn_stutter_data <- NULL
      rv_view_sr$markers_auto <- character(0)
      rv_view_sr$list_dt_average_sr <- NULL

      shinyFileChoose(input, "load_stutter_data", roots = c(documents = paste0(path_pack, "/extdata/stutter")), session = session, filetypes = c("RData"))

      observeEvent(input$load_stutter_data, {
        req(input$load_stutter_data)
        path_stutter_data <- parseFilePaths(roots = c(documents = paste0(path_pack, "/extdata/stutter")), input$load_stutter_data)$datapath

        if(length(path_stutter_data) > 0){
          rv_view_sr$fn_stutter_data <- basename(path_stutter_data)

          load(path_stutter_data)

          #### data_est_sr$list_dt_average_sr

          rv_view_sr$markers_auto <- names(list_dt_average_sr)
          rv_view_sr$list_dt_average_sr <- list_dt_average_sr
        }
      })

      output$fn_stutter_data <- renderUI({
        disabled(textInput(session$ns("fn_stutter_data"), label = "", value = rv_view_sr$fn_stutter_data))
      })

      output$select_mk <- renderUI({
        markers_auto <- rv_view_sr$markers_auto

        if(length(markers_auto) == 0){
          selectInput(session$ns("select_mk"), label = h4("Marker"), choices = markers_auto, selected = NULL, multiple = FALSE, selectize = FALSE)
        }else{
          selectInput(session$ns("select_mk"), label = h4("Marker"), choices = markers_auto, selected = markers_auto[1], multiple = FALSE, selectize = FALSE)
        }
      })

      output$select_st_type <- renderUI({
        fn_stutter_data <- rv_view_sr$fn_stutter_data

        if(is.null(fn_stutter_data)){
          selectInput(session$ns("select_st_type"), label = h4("Stutter type"), choices = character(0), selected = NULL, multiple = FALSE, selectize = FALSE)
        }else{
          selectInput(session$ns("select_st_type"), label = h4("Stutter type"), choices = c("Back stutter", "Forward stutter", "Double-back stutter", "Minus 2-nt stutter"), selected = "Back stutter", multiple = FALSE, selectize = FALSE)
        }
      })

      output$dt_average_sr_display <- renderDataTable({
        select_mk <- input$select_mk
        list_dt_average_sr <- rv_view_sr$list_dt_average_sr

        if(!is.null(list_dt_average_sr) && !is.null(select_mk)){
          dt_average_sr <- list_dt_average_sr[[select_mk]]
          dt_average_sr_display <- dt_average_sr[, !c("Corrected allele (BSR)", "Corrected allele (FSR)", "Corrected allele (DSR)", "Corrected allele (M2SR)"), with = FALSE]

          datatable(
            dt_average_sr_display,
            selection = list(mode = "single", target = "row"),
            options = list(dom = "t", iDisplayLength = nrow(dt_average_sr_display), ordering = FALSE, autoWidth = TRUE, scrollY = "600px", fixedHeader = TRUE),
            rownames = FALSE
          ) %>%
            formatSignif(columns = c("Average BSR", "Average FSR", "Average DSR", "Average M2SR"), digits = 4) %>%
            formatStyle(names(dt_average_sr_display), lineHeight = "80%")
        }
      })
    }
  )
}
