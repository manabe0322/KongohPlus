#' est_sr_result_ui
#'
#' @description The function to create the ui module to create the result tab for estimating average stutter ratios
est_sr_result_ui <- function(id){
  ns <- NS(id)

  tabPanel("Result",
           br(),
           fluidRow(
             column(3,
                    wellPanel(
                      uiOutput(ns("select_mk")),
                      uiOutput(ns("select_st_type"))
                    ),
                    br(),
                    actionButton(ns("act_saveas_est_sr_1"), label = "Save as", class = "btn btn-primary btn-lg")
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

#' est_sr_setting_server
#'
#' @description The function to create the server module to create the result tab for estimating average stutter ratios
est_sr_result_server <- function(id, path_pack, rv_est_sr){
  moduleServer(
    id,
    function(input, output, session){
      data_est_sr <- rv_est_sr$data_est_sr

      if(!is.null(data_est_sr)){
        rv_est_sr_result <- reactiveValues()
        rv_est_sr_result$markers_auto <- names(data_est_sr$list_dt_average_sr)
        rv_est_sr_result$list_dt_average_sr <- data_est_sr$list_dt_average_sr

        output$select_mk <- renderUI({
          markers_auto <- rv_est_sr_result$markers_auto

          if(length(markers_auto) == 0){
            selectInput(session$ns("select_mk"), label = h4("Marker"), choices = markers_auto, selected = NULL, multiple = FALSE, selectize = FALSE)
          }else{
            selectInput(session$ns("select_mk"), label = h4("Marker"), choices = markers_auto, selected = markers_auto[1], multiple = FALSE, selectize = FALSE)
          }
        })

        output$select_st_type <- renderUI({
          fn_stutter_data <- rv_est_sr_result$fn_stutter_data

          if(is.null(fn_stutter_data)){
            selectInput(session$ns("select_st_type"), label = h4("Stutter type"), choices = character(0), selected = NULL, multiple = FALSE, selectize = FALSE)
          }else{
            selectInput(session$ns("select_st_type"), label = h4("Stutter type"), choices = c("Back stutter", "Forward stutter", "Double-back stutter", "Minus 2-nt stutter"), selected = "Back stutter", multiple = FALSE, selectize = FALSE)
          }
        })

        observeEvent(input$act_saveas_est_sr_1, {
          showModal(modalDialog(
            title = "Save stutter data",
            textInput(session$ns("fn_est_sr"), label = "Enter the file name"),
            footer = tagList(
              disabled(actionButton(session$ns("act_saveas_est_sr_2"), "Save")),
              modalButton("Close")
            )
          ))
        })

        observeEvent(input$fn_est_sr, {
          fn_est_sr <- input$fn_est_sr
          sr_file_names <- list.files(paste0(path_pack, "/extdata/stutter"))
          if(is.element(paste0(fn_est_sr, ".RData"), sr_file_names)){
            showFeedbackDanger(inputId = "fn_est_sr", text = "A file with the same name already exists.")
            disable("act_saveas_est_sr_2")
          }else if(nchar(fn_est_sr) == 0){
            hideFeedback("fn_est_sr")
            disable("act_saveas_est_sr_2")
          }else{
            hideFeedback("fn_est_sr")
            enable("act_saveas_est_sr_2")
          }
        })

        observeEvent(input$act_saveas_est_sr_2, {
          fn_est_sr <- input$fn_est_sr
          save(list = c("data_est_sr"), file = paste0(path_pack, "/extdata/stutter/", fn_est_sr, ".RData"))
          removeModal()
        })

        output$dt_average_sr_display <- renderDataTable({
          select_mk <- input$select_mk
          list_dt_average_sr <- rv_est_sr_result$list_dt_average_sr

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
    }
  )
}
