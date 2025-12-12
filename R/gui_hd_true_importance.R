#' hd_true_importance_ui
#'
#' @description The function to create the ui module for the Hd true test
hd_true_importance_ui <- function(id){
  ns <- NS(id)

  tabPanel("Hd true test (importance sampling)",
           useWaiter(),
           titlePanel("Hd true test (importance sampling)"),
           fluidRow(
             column(3,
                    wellPanel(
                      uiOutput(ns("title_csp")),
                      textOutput(ns("fn_csp")),
                      br(),
                      uiOutput(ns("title_af")),
                      textOutput(ns("fn_af"))
                    ),
                    wellPanel(
                      uiOutput(ns("noc")),
                      br(),
                      uiOutput(ns("ui_n_simu")),
                      br(),
                      br(),
                      uiOutput(ns("act_analyze_hd_true")),
                      br(),
                      textOutput(ns("act_reason"))
                    )
             ),
             column(2,
                    dataTableOutput(ns("dt_summary_importance")),
                    br(),
                    br(),
                    uiOutput(ns("ui_download_dt_lr_importance")),
                    br(),
                    uiOutput(ns("ui_download_dt_gt_non_contrib"))
             ),
             column(7,
                    plotOutput(ns("hist_lr_importance"), height = "600px")
             )
           )
  )
}

#' hd_true_importance_server
#'
#' @description The function to create the server module for the Hd true test
hd_true_importance_server <- function(id, path_pack, rv_deconvo){
  moduleServer(
    id,
    function(input, output, session){
      rv_hd_true_importance <- reactiveValues()
      rv_hd_true_importance$dt_lr_importance <- NULL
      rv_hd_true_importance$dt_gt_non_contrib <- NULL
      rv_hd_true_importance$lr_importance <- NULL
      rv_hd_true_importance$bias <- NULL
      rv_hd_true_importance$expect_lr_hd <- NULL
      rv_hd_true_importance$data_deconvo <- NULL
      rv_hd_true_importance$hnc_all <- NULL
      rv_hd_true_importance$theta <- NULL
      rv_hd_true_importance$fn_csp <- NULL
      rv_hd_true_importance$fn_af <- NULL
      rv_hd_true_importance$data_hnc_list <- NULL

      observe({
        req(rv_deconvo)
        req(rv_hd_true_importance)

        rv_hd_true_importance$data_deconvo <- data_deconvo <- rv_deconvo$data_deconvo
        rv_hd_true_importance$hnc_all <- data_deconvo$hnc_all
        rv_hd_true_importance$theta <- data_deconvo$theta
        rv_hd_true_importance$fn_csp <- data_deconvo$fn_csp
        rv_hd_true_importance$fn_af <- data_deconvo$fn_af
        rv_hd_true_importance$data_hnc_list <- data_deconvo$data_hnc_list

        bool_perform_deconvo <- rv_deconvo$bool_perform_deconvo
        if(bool_perform_deconvo){
          rv_hd_true_importance$dt_lr_importance <- NULL
          rv_hd_true_importance$dt_gt_non_contrib <- NULL
          rv_hd_true_importance$lr_importance <- NULL
          rv_hd_true_importance$bias <- NULL
          rv_hd_true_importance$expect_lr_hd <- NULL
        }
      })

      output$title_csp <- renderUI({
        req(rv_hd_true_importance$fn_csp)
        h4("Crime stain profile")
      })

      output$fn_csp <- renderText({
        req(rv_hd_true_importance$fn_csp)
        paste0(rv_hd_true_importance$fn_csp)
      })

      output$title_af <- renderUI({
        req(rv_hd_true_importance$fn_af)
        h4("Allele frequency")
      })

      output$fn_af <- renderText({
        req(rv_hd_true_importance$fn_af)
        paste0(rv_hd_true_importance$fn_af)
      })

      output$noc <- renderUI({
        req(rv_hd_true_importance$hnc_all)
        hnc_all <- rv_hd_true_importance$hnc_all
        selectInput(session$ns("noc"), label = "Number of contributors", choices = hnc_all, selected = hnc_all[1], selectize = FALSE)
      })

      output$ui_n_simu <- renderUI({
        req(rv_hd_true_importance$hnc_all)
        numericInput(session$ns("n_simu"), label = "Number of simulations", value = 100)
      })

      observeEvent(input$n_simu, {
        req(rv_hd_true_importance$theta)
        theta <- rv_hd_true_importance$theta
        n_simu <- input$n_simu

        if(!is.integer(n_simu) || n_simu <= 0){
          showFeedbackDanger(inputId = "n_simu", text = "The value must be a positive integer.")
          bool_n_simu <- FALSE
        }else{
          hideFeedback("n_simu")
          bool_n_simu <- TRUE
        }

        if(bool_n_simu && theta == 0){
          enable("act_analyze_hd_true")
        }else{
          disable("act_analyze_hd_true")
        }
      })

      output$act_analyze_hd_true <- renderUI({
        req(rv_hd_true_importance$theta)
        theta <- rv_hd_true_importance$theta

        if(theta == 0){
          actionButton(session$ns("act_analyze_hd_true"), label = "Analyze", class = "btn btn-primary btn-lg")
        }else{
          disabled(actionButton(session$ns("act_analyze_hd_true"), label = "Analyze", class = "btn btn-primary btn-lg"))
        }
      })

      observeEvent(input$act_analyze_hd_true, {
        req(rv_hd_true_importance$hnc_all)
        hnc_all <- rv_hd_true_importance$hnc_all

        waiter_ui <- Waiter$new(html = tagList(spin_3k(),
                                               h3("Checking the settings...")),
                                color = transparent(.5))
        waiter_ui$show()

        hnc <- as.numeric(input$noc)
        n_simu <- input$n_simu
        tmp <- try(perform_hd_true(hnc, rv_hd_true_importance$data_deconvo, n_simu, dt_gt_non_contrib = NULL, importance_sampling = TRUE, waiter_ui))
        if(inherits(tmp, "try-error")){
          waiter_ui$hide()
          showModal(modalDialog(title = "Error", "An unexpected error has occurred.", easyClose = TRUE, footer = NULL))
        }else{
          rv_hd_true_importance$dt_lr_importance <- dt_lr_importance <- tmp$dt_lr_rand

          if(is.null(dt_lr_importance)){
            data_hnc_list <- rv_hd_true_importance$data_hnc_list
            data_h1 <- data_hnc_list[[which(hnc_all == hnc)]]
            error_message_h1 <- data_h1$error_message
            waiter_ui$hide()
            showModal(modalDialog(title = "Error", HTML(error_message_h1), easyClose = TRUE, footer = NULL))
          }else{
            rv_hd_true_importance$lr_importance <- lr_importance <- as.numeric(dt_lr_importance[nrow(dt_lr_importance), .SD, .SDcols = grep("Sample", names(dt_lr_importance))])
            rv_hd_true_importance$bias <- bias <- tmp$bias
            rv_hd_true_importance$dt_gt_non_contrib <- tmp$dt_gt_non_contrib
            rv_hd_true_importance$expect_lr_hd <- sum(lr_importance * bias) / n_simu
            waiter_ui$hide()
          }
        }
      })

      output$act_reason <- renderText({
        req(rv_hd_true_importance$theta)
        theta <- rv_hd_true_importance$theta

        if(theta == 0){
          ""
        }else{
          "This feature is available only when theta is set to 0."
        }
      })

      output$dt_summary_importance <- renderDataTable({
        req(rv_hd_true_importance$expect_lr_hd)
        expect_lr_hd <- rv_hd_true_importance$expect_lr_hd

        dt_summary_importance <- data.table("Average Hd true LR" = expect_lr_hd)

        datatable(
          dt_summary_importance,
          selection = "none",
          options = list(dom = "t", iDisplayLength = nrow(dt_summary_importance), ordering = FALSE, autowidth = TRUE),
          rownames = FALSE
        )
      })

      output$ui_download_dt_lr_importance <- renderUI({
        req(rv_hd_true_importance$dt_lr_importance)
        downloadButton(session$ns("download_dt_lr_importance"), "Download LRs")
      })

      output$download_dt_lr_importance <- downloadHandler(
        filename = paste0(gsub(" ", "_", format(as.POSIXct(Sys.time()), "%Y-%m-%d %H%M%S")), "_Hd-true_LR_importance.csv"),
        content = function(file){
          write.csv(rv_hd_true_importance$dt_lr_importance, file, row.names = FALSE)
        }
      )

      output$ui_download_dt_gt_non_contrib <- renderUI({
        req(rv_hd_true_importance$dt_gt_non_contrib)
        downloadButton(session$ns("download_dt_gt_non_contrib"), "Download genotypes")
      })

      output$download_dt_gt_non_contrib <- downloadHandler(
        filename = paste0(gsub(" ", "_", format(as.POSIXct(Sys.time()), "%Y-%m-%d %H%M%S")), "_Hd-true_genotypes.csv"),
        content = function(file){
          write.csv(rv_hd_true_importance$dt_gt_non_contrib, file, row.names = FALSE)
        }
      )

      output$hist_lr_importance <- renderPlot({
        req(rv_hd_true_importance$lr_importance)
        lr_importance <- rv_hd_true_importance$lr_importance
        log10_lr_importance <- log10(lr_importance)
        log10_lr_importance[log10_lr_importance < -300] <- -300
        hist(log10_lr_importance)
      })

      return(rv_hd_true_importance)
    }
  )
}
