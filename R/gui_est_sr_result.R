#' est_sr_result_ui
#'
#' @description The function to create the ui module to create the result tab for estimating average stutter ratios
est_sr_result_ui <- function(id){
  ns <- NS(id)

  tabPanel("Result",
           br(),
           fluidRow(
             column(2,
                    wellPanel(
                      uiOutput(ns("select_mk")),
                      uiOutput(ns("select_st_type"))
                    )
             ),
             column(5,
                    wellPanel(
                      dataTableOutput(ns("dt_ave_sr"))
                    )
             ),
             column(5,
                    wellPanel(),
                    wellPanel(
                      dataTableOutput(ns("dt_diag_sr"))
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

        dt_setting_est_sr <- rv_est_sr$dt_setting_est_sr

        dt_result_est_bsr <- data_est_sr$dt_result_est_bsr
        trace_bsr <- data_est_sr$trace_bsr
        dt_average_bsr <- data_est_sr$dt_average_bsr
        dt_result_est_fsr <- data_est_sr$dt_result_est_fsr
        trace_fsr <- data_est_sr$trace_fsr
        dt_average_fsr <- data_est_sr$dt_average_fsr
        dt_result_est_dsr <- data_est_sr$dt_result_est_dsr
        trace_dsr <- data_est_sr$trace_dsr
        dt_average_dsr <- data_est_sr$dt_average_dsr
        dt_result_est_m2sr <- data_est_sr$dt_result_est_m2sr
        trace_m2sr <- data_est_sr$trace_m2sr
        dt_average_m2sr <- data_est_sr$dt_average_m2sr

        output$select_mk <- renderUI({
          ns <- session$ns
          markers_auto <- dt_setting_est_sr$Marker
          selectInput(ns("select_mk"), label = h4("Marker"), choices = markers_auto, selected = markers_auto[1], multiple = FALSE, selectize = FALSE)
        })

        output$select_st_type <- renderUI({
          ns <- session$ns
          selectInput(ns("select_st_type"), label = h4("Stutter type"), choices = c("Back stutter", "Forward stutter", "Double-back stutter", "Minus 2-nt stutter"), selected = "Back stutter", multiple = FALSE, selectize = FALSE)
        })

        observeEvent(input$select_mk, {
          select_mk <- input$select_mk
          select_st_type <- input$select_st_type
          if(select_st_type == "Back stutter"){
            dt_ave_sr <- dt_average_bsr[[select_mk]]
          }else if(select_st_type == "Forward stutter"){
            dt_ave_sr <- dt_average_fsr[[select_mk]]
          }else if(select_st_type == "Double-back stutter"){
            dt_ave_sr <- dt_average_dsr[[select_mk]]
          }else if(select_st_type == "Minus 2-nt stutter"){
            dt_ave_sr <- dt_average_m2sr[[select_mk]]
          }

          output$dt_ave_sr <- renderDataTable(server = FALSE, {
            datatable(
              dt_ave_sr,
              selection = "none",
              options = list(iDisplayLength = 10, ordering = FALSE, autowidth = TRUE),
              rownames = FALSE,
              colnames = c("Allele", "Corrected", "Average SR")
            )
          })
        })
      }
    }
  )
}
