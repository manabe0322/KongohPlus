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
        fn_exp <- rv_est_par$fn_exp
        fn_gt <- rv_est_par$fn_gt
        dt_setting_est_par <- rv_est_par$dt_setting_est_par
        list_dt_average_sr <- rv_est_par$list_dt_average_sr
        name_kit <- rv_est_par$name_kit
        name_sr <- rv_est_par$name_sr
        warmup <- rv_est_par$warmup
        sampling <- rv_est_par$sampling
        chains <- rv_est_par$chains

        lp_ae <- data_est_par$lp_ae
        alpha <- data_est_par$alpha
        rhat_lp_ae <- data_est_par$rhat_lp_ae
        rhat_alpha <- data_est_par$rhat_alpha
        traceplot_ae <- data_est_par$traceplot_ae
        lp_gamma <- data_est_par$lp_gamma
        cv_allele <- data_est_par$cv_allele
        cv_stutter <- data_est_par$cv_stutter
        rhat_lp_gamma <- data_est_par$rhat_lp_gamma
        rhat_cv_allele <- data_est_par$rhat_cv_allele
        rhat_cv_stutter <- data_est_par$rhat_cv_stutter
        traceplot_cv_allele <- data_est_par$traceplot_cv_allele
        traceplot_cv_stutter <- data_est_par$traceplot_cv_stutter

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

          save(list = c("fn_exp", "fn_gt", "dt_setting_est_par", "list_dt_average_sr", "name_kit", "name_sr", "warmup", "sampling", "chains",
                        "lp_ae", "alpha", "rhat_lp_ae", "rhat_alpha", "traceplot_ae",
                        "lp_gamma", "cv_allele", "cv_stutter", "rhat_lp_gamma", "rhat_cv_allele", "rhat_cv_stutter", "traceplot_cv_allele", "traceplot_cv_stutter"),
               file = paste0(path_pack, "/extdata/param/", fn_est_par, ".RData"))

          removeModal()
        })
      }
    }
  )
}
