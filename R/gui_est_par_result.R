#' est_par_result_ui
#'
#' @description The function to create the ui module to create the result tab for estimating parameters
est_par_result_ui <- function(id){
  ns <- NS(id)

  tabPanel("Result",
           br(),
           fluidRow(
             column(2,
                    h4("Files and settings"),
                    br(),
                    actionButton(ns("act_saveas_est_par_1"), label = "Save as", class = "btn btn-primary btn-lg")
             ),
             column(5,
                    wellPanel(
                      h4("Locus-specific amplification efficiency")
                    ),
                    wellPanel(
                      h4("Gamma distribution")
                    )
             ),
             column(5,
                    wellPanel(
                      h4("Diagnostics")
                    )
             )
           )
  )
}

#' est_par_result_server
#'
#' @description The function to create the server module to create the result tab for estimating parameters
est_par_result_server <- function(id, path_pack, rv_est_par){
  moduleServer(
    id,
    function(input, output, session){
      data_est_par <- rv_est_par$data_est_par

      if(!is.null(data_est_par)){
        median_alpha <- data_est_par$median_alpha
        rhat_alpha <- data_est_par$rhat_alpha
        median_lp_ae <- data_est_par$median_lp_ae
        rhat_lp_ae <- data_est_par$rhat_lp_ae
        traceplot_ae <- data_est_par$traceplot_ae
        median_lp_gamma <- data_est_par$median_lp_gamma
        median_cv <- data_est_par$median_cv
        rhat_lp_gamma <- data_est_par$rhat_lp_gamma
        rhat_cv <- data_est_par$rhat_cv
        traceplot_cv <- data_est_par$traceplot_cv
        list_dt_average_sr <- rv_est_par$list_dt_average_sr
        fn_stutter_file <- rv_est_par$fn_stutter_file

        rv_par_file_names <- reactiveVal(NULL)

        observeEvent(input$act_saveas_est_par_1, {
          par_file_names <- list.files(paste0(path_pack, "/extdata/param"))
          rv_par_file_names(gsub(".{6}$", "", par_file_names))

          showModal(modalDialog(
            title = "Save parameter data",
            textInput(session$ns("fn_est_par"), label = "Enter the file name"),
            footer = tagList(
              disabled(actionButton(session$ns("act_saveas_est_par_2"), "Save")),
              modalButton("Close")
            )
          ))
        })

        observeEvent(input$fn_est_par, {
          fn_est_par <- input$fn_est_par
          if(nchar(fn_est_par) == 0){
            showFeedbackDanger(inputId = "fn_est_par", text = "Required")
            disable("act_saveas_est_par_2")
          }else if(is.element(fn_est_par, rv_par_file_names())){
            showFeedbackDanger(inputId = "fn_est_par", text = "This name is already registered.")
            disable("act_saveas_est_par_2")
          }else{
            hideFeedback("fn_est_par")
            enable("act_saveas_est_par_2")
          }
        })

        observeEvent(input$act_saveas_est_par_2, {
          fn_est_par <- input$fn_est_par

          save(list = c("median_alpha", "rhat_alpha", "median_lp_ae", "rhat_lp_ae", "traceplot_ae",
                        "median_lp_gamma", "median_cv", "rhat_lp_gamma", "rhat_cv", "traceplot_cv",
                        "list_dt_average_sr", "fn_stutter_file"),
               file = paste0(path_pack, "/extdata/param/", fn_est_par, ".RData"))

          removeModal()
        })
      }
    }
  )
}
