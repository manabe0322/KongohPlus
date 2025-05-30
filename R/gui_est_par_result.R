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
        fn_exp_param <- rv_est_par$fn_exp_param
        fn_gt_param <- rv_est_par$fn_gt_param
        dt_setting_est_par <- rv_est_par$dt_setting_est_par
        list_dt_average_sr <- rv_est_par$list_dt_average_sr
        name_kit_param <- rv_est_par$name_kit_param
        name_sr <- rv_est_par$name_sr
        min_alpha <- rv_est_par$min_alpha
        max_alpha <- rv_est_par$max_alpha
        min_cv_allele <- rv_est_par$min_cv_allele
        max_cv_allele <- rv_est_par$max_cv_allele
        min_cv_stutter <- rv_est_par$min_cv_stutter
        max_cv_stutter <- rv_est_par$max_cv_stutter
        initial_point_factor <- rv_est_par$initial_point_factor
        alpha <- data_est_par$alpha
        expect_ae <- data_est_par$expect_ae
        cv_allele <- data_est_par$cv_allele
        cv_stutter <- data_est_par$cv_stutter

        observeEvent(input$act_saveas_est_par_1, {
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
          par_file_names <- list.files(paste0(path_pack, "/extdata/param"))
          if(is.element(paste0(fn_est_par, ".RData"), par_file_names)){
            showFeedbackDanger(inputId = "fn_est_par", text = "A file with the same name already exists.")
            disable("act_saveas_est_par_2")
          }else if(nchar(fn_est_par) == 0){
            hideFeedback("fn_est_par")
            disable("act_saveas_est_par_2")
          }else{
            hideFeedback("fn_est_par")
            enable("act_saveas_est_par_2")
          }
        })

        observeEvent(input$act_saveas_est_par_2, {
          fn_est_par <- input$fn_est_par

          save(list = c("fn_exp_param", "fn_gt_param", "dt_setting_est_par", "list_dt_average_sr", "name_kit_param", "name_sr", "alpha", "expect_ae", "cv_allele", "cv_stutter",
                        "min_alpha", "max_alpha", "min_cv_allele", "max_cv_allele", "min_cv_stutter", "max_cv_stutter", "initial_point_factor"),
               file = paste0(path_pack, "/extdata/param/", fn_est_par, ".RData"))

          removeModal()
        })
      }
    }
  )
}
