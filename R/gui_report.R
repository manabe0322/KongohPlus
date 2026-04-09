#' report_ui
#'
#' @description The function to create the ui module for generating a report
#' @param id The ui ID
report_ui <- function(id){
  ns <- NS(id)

  tabPanel("Report",
           useWaiter(),
           textInput(ns("case_name"), label = h4("Case name")),
           textInput(ns("user_name"), label = h4("User name")),
           numericInput(ns("min_gt_prob_display"), label = h4("Minimum genotype probability to display"), value = 0.001, step = 0.001),
           numericInput(ns("max_gtc_display"), label = h4("Maximum number of genotypes to display"), value = 20),
           br(),
           disabled(downloadButton(ns("download_report"), "Report", class = "btn btn-primary btn-lg"))
  )
}

#' report_server
#'
#' @description The function to create the server module for generating a report
#' @param id The server ID
#' @param path_pack The path of the package KongohPlus
#' @param rv_deconvo The reactive values derived from the setting server for deconvolution
#' @param rv_lr The reactive values derived from the LR server
report_server <- function(id, path_pack, rv_deconvo, rv_lr){
  moduleServer(
    id,
    function(input, output, session){
      rv_report <- reactiveValues()
      rv_report$ver_soft <- NULL
      rv_report$dt_csp <- NULL
      rv_report$fn_csp <- NULL
      rv_report$fn_af <- NULL
      rv_report$dt_anameth <- NULL
      rv_report$fn_anameth <- NULL
      rv_report$markers_use <- NULL
      rv_report$data_hnc_list <- NULL
      rv_report$dt_ignore_mk <- NULL
      rv_report$hnc_all <- NULL
      rv_report$message_sr <- NULL
      rv_report$min_size <- NULL
      rv_report$max_size <- NULL
      rv_report$min_temp <- 1
      rv_report$max_temp <- 32000
      rv_report$min_rep_ae <- 0.2
      rv_report$max_rep_ae <- 5
      rv_report$min_deg <- 0.0001
      rv_report$max_iter <- 500
      rv_report$num_consistent <- 5
      rv_report$max_mle_runs <- 10
      rv_report$max_diff_like <- 0.005
      rv_report$theta <- NULL
      rv_report$dt_ref <- NULL
      rv_report$fn_ref <- NULL
      rv_report$dt_lr_calc_record <- NULL
      rv_report$lr_list <- NULL
      rv_report$bool_correct_min_gt_prob_display <- TRUE
      rv_report$bool_correct_max_gtc_display <- TRUE

      observe({
        req(rv_deconvo)
        data_deconvo <- rv_deconvo$data_deconvo

        if(!is.null(data_deconvo)){
          rv_report$ver_soft <- data_deconvo$ver_soft
          rv_report$dt_csp <- data_deconvo$dt_csp
          rv_report$fn_csp <- data_deconvo$fn_csp
          rv_report$fn_af <- data_deconvo$fn_af
          rv_report$dt_anameth <- data_deconvo$dt_anameth
          rv_report$fn_anameth <- data_deconvo$fn_anameth
          rv_report$markers_use <- data_deconvo$markers_use
          rv_report$data_hnc_list <- data_deconvo$data_hnc_list
          rv_report$dt_ignore_mk <- data_deconvo$dt_ignore_mk
          rv_report$hnc_all <- data_deconvo$hnc_all
          rv_report$message_sr <- data_deconvo$message_sr
          rv_report$min_size <- data_deconvo$min_size
          rv_report$max_size <- data_deconvo$max_size
          rv_report$min_temp <- data_deconvo$min_temp
          rv_report$max_temp <- data_deconvo$max_temp
          rv_report$min_rep_ae <- data_deconvo$min_rep_ae
          rv_report$max_rep_ae <- data_deconvo$max_rep_ae
          rv_report$min_deg <- data_deconvo$min_deg
          rv_report$max_iter <- data_deconvo$max_iter
          rv_report$num_consistent <- data_deconvo$num_consistent
          rv_report$max_mle_runs <- data_deconvo$max_mle_runs
          rv_report$max_diff_like <- data_deconvo$max_diff_like
          rv_report$theta <- data_deconvo$theta

          enable("download_report")
        }else{
          disable("download_report")
        }
      })

      observe({
        req(rv_lr)
        rv_report$dt_ref <- rv_lr$dt_ref
        rv_report$fn_ref <- rv_lr$fn_ref
        rv_report$dt_lr_calc_record <- rv_lr$dt_lr_calc_record
        rv_report$lr_list <- rv_lr$lr_list
      })

      observeEvent(input$min_gt_prob_display, {
        min_gt_prob_display <- input$min_gt_prob_display

        hideFeedback("min_gt_prob_display")
        bool_correct_min_gt_prob_display <- TRUE

        if(!is.numeric(min_gt_prob_display) || min_gt_prob_display < 0 || min_gt_prob_display > 1){
          showFeedbackDanger(inputId = "min_gt_prob_display", text = "Value must be between 0 and 1.")
          bool_correct_min_gt_prob_display <- FALSE
        }

        rv_report$bool_correct_min_gt_prob_display <- bool_correct_min_gt_prob_display
      })

      observeEvent(input$max_gtc_display, {
        max_gtc_display <- input$max_gtc_display

        hideFeedback("max_gtc_display")
        bool_correct_max_gtc_display <- TRUE

        if(!is.integer(max_gtc_display) || max_gtc_display < 1){
          showFeedbackDanger(inputId = "max_gtc_display", text = "Value must be a positive integer.")
          bool_correct_max_gtc_display <- FALSE
        }

        rv_report$bool_correct_max_gtc_display <- bool_correct_max_gtc_display
      })

      observe({
        req(rv_report)
        req(rv_deconvo)
        data_deconvo <- rv_deconvo$data_deconvo

        if(!is.null(data_deconvo)){
          if(all(rv_report$bool_correct_min_gt_prob_display, rv_report$bool_correct_max_gtc_display)){
            enable("download_report")
          }else{
            disable("download_report")
          }
        }else{
          disable("download_report")
        }
      })

      output$download_report <- downloadHandler(
        filename = function() {
          paste0("report_", gsub(" ", "_", format(as.POSIXct(Sys.time()), "%Y-%m-%d %H%M%S")), "_", input$case_name, "_", input$user_name, ".html")
        },
        content = function(file){
          req(rv_deconvo$data_deconvo)

          message_pandoc <- check_pandoc()
          if(nzchar(message_pandoc)){
            showModal(modalDialog(title = "Error", HTML(message_pandoc), easyClose = TRUE, footer = NULL))
          }else{
            waiter_ui <- Waiter$new(html = tagList(spin_3k(),
                                                   h3("Generating the report...")),
                                    color = transparent(.5))

            waiter_ui$show()

            # Temporary directory for work
            tmpdir <- tempfile("report_")
            dir.create(tmpdir, recursive = TRUE)

            # Copy Rmd to a temporary directory (render requires the working directory)
            rmd_tmp <- file.path(tmpdir, "report.Rmd")
            file.copy(paste0(path_pack, "/rmd/report.Rmd"), rmd_tmp, overwrite = TRUE)

            # Define params
            params <- list(
              case_name = input$case_name,
              user_name = input$user_name,
              min_gt_prob_display = input$min_gt_prob_display,
              max_gtc_display = input$max_gtc_display,
              ver_soft = rv_report$ver_soft,
              dt_csp = rv_report$dt_csp,
              fn_csp = rv_report$fn_csp,
              fn_af = rv_report$fn_af,
              dt_anameth = rv_report$dt_anameth,
              fn_anameth = rv_report$fn_anameth,
              markers_use = rv_report$markers_use,
              data_hnc_list = rv_report$data_hnc_list,
              dt_ignore_mk = rv_report$dt_ignore_mk,
              hnc_all = rv_report$hnc_all,
              message_sr = rv_report$message_sr,
              min_size = rv_report$min_size,
              max_size = rv_report$max_size,
              min_temp = rv_report$min_temp,
              max_temp = rv_report$max_temp,
              min_rep_ae = rv_report$min_rep_ae,
              max_rep_ae = rv_report$max_rep_ae,
              min_deg = rv_report$min_deg,
              max_iter = rv_report$max_iter,
              num_consistent = rv_report$num_consistent,
              max_mle_runs = rv_report$max_mle_runs,
              max_diff_like = rv_report$max_diff_like,
              theta = rv_report$theta,
              dt_ref = rv_report$dt_ref,
              fn_ref = rv_report$fn_ref,
              dt_lr_calc_record = rv_report$dt_lr_calc_record,
              lr_list = rv_report$lr_list
            )

            # Generate HTML
            css_path <- system.file("rmd/report.css", package = "KongohPlus")
            out_html <- render(
              input = rmd_tmp,
              output_format = html_document(self_contained = TRUE, css = css_path),
              params = params,
              envir = new.env(parent = globalenv()),
              quiet = TRUE
            )

            ok <- file.copy(out_html, file, overwrite = TRUE)
            if(!ok){
              stop("Failed to copy rendered HTML to download target.")
            }

            waiter_ui$hide()
          }
        }
      )
    }
  )
}
