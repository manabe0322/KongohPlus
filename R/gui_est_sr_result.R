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
                      uiOutput(ns("select_st_type")),
                      h4("Model"),
                      textOutput(ns("meth_sr"))
                    ),
                    br(),
                    actionButton(ns("act_saveas_est_sr_1"), label = "Save as", class = "btn btn-primary btn-lg")
             ),
             column(5,
                    wellPanel(
                      dataTableOutput(ns("dt_ave_sr"))
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

#' full_join_average_sr
#'
#' @description The function to combine data.tables of average stutter ratios in each stutter type into one data.table (per marker)
#' @param list_dt_average_bsr A list of data.tables of average back stutter ratios in each marker
#' @param list_dt_average_fsr A list of data.tables of average forward stutter ratios in each marker
#' @param list_dt_average_dsr A list of data.tables of average double-back stutter ratios in each marker
#' @param list_dt_average_m2sr A list of data.tables of average minus 2-nt stutter ratios in each marker
full_join_average_sr <- function(list_dt_average_bsr, list_dt_average_fsr, list_dt_average_dsr, list_dt_average_m2sr){
  markers_auto <- names(list_dt_average_bsr)
  n_mk <- length(markers_auto)
  list_dt_average_sr <- vector(mode = "list", length = n_mk)
  names(list_dt_average_sr) <- markers_auto

  for(i in 1:n_mk){
    mk <- markers_auto[i]
    dt_ave_bsr <- list_dt_average_bsr[[mk]]
    names(dt_ave_bsr) <- c("Allele", "Corrected allele (BSR)", "Average BSR")
    dt_ave_fsr <- list_dt_average_fsr[[mk]]
    names(dt_ave_fsr) <- c("Allele", "Corrected allele (FSR)", "Average FSR")
    dt_ave_dsr <- list_dt_average_dsr[[mk]]
    names(dt_ave_dsr) <- c("Allele", "Corrected allele (DSR)", "Average DSR")
    dt_ave_m2sr <- list_dt_average_m2sr[[mk]]
    names(dt_ave_m2sr) <- c("Allele", "Corrected allele (M2SR)", "Average M2SR")

    dt_ave_sr <- copy(dt_ave_bsr)
    dt_ave_sr <- full_join(dt_ave_sr, dt_ave_fsr, by = c("Allele"))
    dt_ave_sr <- full_join(dt_ave_sr, dt_ave_dsr, by = c("Allele"))
    dt_ave_sr <- full_join(dt_ave_sr, dt_ave_m2sr, by = c("Allele"))

    dt_ave_sr[["Average BSR"]][is.na(dt_ave_sr[["Average BSR"]])] <- 0
    dt_ave_sr[["Average FSR"]][is.na(dt_ave_sr[["Average FSR"]])] <- 0
    dt_ave_sr[["Average DSR"]][is.na(dt_ave_sr[["Average DSR"]])] <- 0
    dt_ave_sr[["Average M2SR"]][is.na(dt_ave_sr[["Average M2SR"]])] <- 0

    list_dt_average_sr[[i]] <- dt_ave_sr
  }

  return(list_dt_average_sr)
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
        fn_exp_sr <- rv_est_sr$fn_exp_sr
        fn_gt_sr <- rv_est_sr$fn_gt_sr
        dt_setting_est_sr <- rv_est_sr$dt_setting_est_sr
        markers_auto <- dt_setting_est_sr$Marker
        name_kit_stutter <- rv_est_sr$name_kit_stutter
        min_slope <- rv_est_sr$min_slope
        max_slope <- rv_est_sr$max_slope
        min_intercept <- rv_est_sr$min_intercept
        max_intercept <- rv_est_sr$max_intercept
        min_myu <- rv_est_sr$min_myu
        max_myu <- rv_est_sr$max_myu
        min_sigma <- rv_est_sr$min_sigma
        max_sigma <- rv_est_sr$max_sigma
        initial_point_factor <- rv_est_sr$initial_point_factor

        dt_result_est_bsr <- data_est_sr$dt_result_est_bsr
        list_dt_average_bsr <- data_est_sr$list_dt_average_bsr
        dt_result_est_fsr <- data_est_sr$dt_result_est_fsr
        list_dt_average_fsr <- data_est_sr$list_dt_average_fsr
        dt_result_est_dsr <- data_est_sr$dt_result_est_dsr
        list_dt_average_dsr <- data_est_sr$list_dt_average_dsr
        dt_result_est_m2sr <- data_est_sr$dt_result_est_m2sr
        list_dt_average_m2sr <- data_est_sr$list_dt_average_m2sr

        output$select_mk <- renderUI({
          ns <- session$ns
          selectInput(ns("select_mk"), label = h4("Marker"), choices = markers_auto, selected = markers_auto[1], multiple = FALSE, selectize = FALSE)
        })

        output$select_st_type <- renderUI({
          ns <- session$ns
          selectInput(ns("select_st_type"), label = h4("Stutter type"), choices = c("Back stutter", "Forward stutter", "Double-back stutter", "Minus 2-nt stutter"), selected = "Back stutter", multiple = FALSE, selectize = FALSE)
        })

        observeEvent(c(input$select_mk,input$select_st_type), {
          select_mk <- input$select_mk
          select_st_type <- input$select_st_type
          index_mk <- which(markers_auto == select_mk)
          if(select_st_type == "Back stutter"){
            dt_ave_sr <- list_dt_average_bsr[[select_mk]]
            dt_result_est_sr <- dt_result_est_bsr
          }else if(select_st_type == "Forward stutter"){
            dt_ave_sr <- list_dt_average_fsr[[select_mk]]
            dt_result_est_sr <- dt_result_est_fsr
          }else if(select_st_type == "Double-back stutter"){
            dt_ave_sr <- list_dt_average_dsr[[select_mk]]
            dt_result_est_sr <- dt_result_est_dsr
          }else if(select_st_type == "Minus 2-nt stutter"){
            dt_ave_sr <- list_dt_average_m2sr[[select_mk]]
            dt_result_est_sr <- dt_result_est_m2sr
          }

          meth_sr <- dt_result_est_sr$Model[index_mk]
          output$meth_sr <- renderText({paste0(meth_sr)})

          output$dt_ave_sr <- renderDataTable(server = FALSE, {
            datatable(
              dt_ave_sr,
              selection = "none",
              options = list(iDisplayLength = 10, ordering = FALSE, autowidth = TRUE),
              rownames = FALSE,
              colnames = c("Allele", "Corrected", "Average SR")
            ) %>%
              formatSignif(columns = c("AlCor", "AveSR"), digits = 3)
          })

          if(meth_sr == "Not consider"){
            param <- character(0)
            val <- numeric(0)
          }else if(is.element(meth_sr, c("Uniform", "Uniform (Multiple loci)"))){
            param <- c("Myu", "Sigma")
            val <- c(dt_result_est_sr$Myu[index_mk], dt_result_est_sr$Sigma[index_mk])
          }else{
            param <- c("Slope", "Intercept", "Sigma")
            val <- c(dt_result_est_sr$Slope[index_mk], dt_result_est_sr$Intercept[index_mk], dt_result_est_sr$Sigma[index_mk])
          }
          dt_param_sr <- data.table(Parameter = param,
                                   EstVar = val)
          names(dt_param_sr) <- c("Parameter", "Estimated value")

          output$dt_param_sr <- renderDataTable(server = FALSE, {
            datatable(
              dt_param_sr,
              selection = "none",
              options = list(dom = "t", iDisplayLength = nrow(dt_param_sr), ordering = FALSE, autowidth = TRUE),
              rownames = FALSE
            ) %>%
              formatSignif(columns = c("Estimated value"), digits = 3)
          })
        })

        rv_sr_file_names <- reactiveVal(NULL)

        observeEvent(input$act_saveas_est_sr_1, {
          sr_file_names <- list.files(paste0(path_pack, "/extdata/stutter_file"))
          rv_sr_file_names(gsub(".{6}$", "", sr_file_names))

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
          if(nchar(fn_est_sr) == 0){
            showFeedbackDanger(inputId = "fn_est_sr", text = "Required")
            disable("act_saveas_est_sr_2")
          }else if(is.element(fn_est_sr, rv_sr_file_names())){
            showFeedbackDanger(inputId = "fn_est_sr", text = "This name is already registered.")
            disable("act_saveas_est_sr_2")
          }else{
            hideFeedback("fn_est_sr")
            enable("act_saveas_est_sr_2")
          }
        })

        observeEvent(input$act_saveas_est_sr_2, {
          fn_est_sr <- input$fn_est_sr
          list_dt_average_sr <- full_join_average_sr(list_dt_average_bsr, list_dt_average_fsr, list_dt_average_dsr, list_dt_average_m2sr)

          save(list = c("fn_exp_sr", "fn_gt_sr", "dt_setting_est_sr", "name_kit_stutter",
                        "min_slope", "max_slope", "min_intercept", "max_intercept", "min_myu", "max_myu", "min_sigma", "max_sigma", "initial_point_factor",
                        "dt_result_est_bsr", "dt_result_est_fsr", "dt_result_est_dsr", "dt_result_est_m2sr", "list_dt_average_sr"),
               file = paste0(path_pack, "/extdata/stutter/", fn_est_sr, ".RData"))

          removeModal()
        })
      }
    }
  )
}
