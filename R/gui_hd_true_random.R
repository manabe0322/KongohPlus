#' hd_true_random_ui
#'
#' @description The function to create the ui module for the Hd true test
hd_true_random_ui <- function(id){
  ns <- NS(id)
  hd_true_random_ui_sub(ns("hd_true_random_sub"))
}

#' hd_true_random_ui_sub
#'
#' @description The function to create the sub-ui module for the Hd true test
hd_true_random_ui_sub <- function(id){
  ns <- NS(id)

  tabPanel("Hd true test (random genotypes)",
           useWaiter(),
           titlePanel("Hd true test (random genotypes)"),
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
                      uiOutput(ns("title_gt_non_contrib")),
                      uiOutput(ns("check_random")),
                      uiOutput(ns("file_gt_non_contrib")),
                      uiOutput(ns("ui_n_simu")),
                      br(),
                      br(),
                      uiOutput(ns("act_analyze_hd_true"))
                    )
             ),
             column(2,
                    dataTableOutput(ns("dt_summary_lr_rand")),
                    br(),
                    br(),
                    uiOutput(ns("ui_download_dt_lr_rand")),
                    br(),
                    uiOutput(ns("ui_download_dt_gt_non_contrib"))
             ),
             column(7,
                    plotOutput(ns("hist_lr_rand"), height = "600px")
             )
           )
  )
}

#' hd_true_random_server
#'
#' @description The function to create the server module for the Hd true test
hd_true_random_server <- function(id, path_pack, rv_deconvo){
  moduleServer(
    id,
    function(input, output, session){
      rv_hd_true_random <- reactiveValues()
      rv_hd_true_random$bool_n_simu <- TRUE
      rv_hd_true_random$lr_rand <- NULL
      rv_hd_true_random$dt_lr_rand <- NULL
      rv_hd_true_random$dt_gt_non_contrib <- NULL

      observe({
        req(rv_deconvo)
        req(rv_hd_true_random)
        rv_hd_true_random <- hd_true_random_server_sub("hd_true_random_sub", path_pack, rv_deconvo, rv_hd_true_random)
      })

      return(rv_hd_true_random)
    }
  )
}

#' hd_true_random_server_sub
#'
#' @description The function to create the sub-server module for the Hd true test
hd_true_random_server_sub <- function(id, path_pack, rv_deconvo, rv_hd_true_random){
  moduleServer(
    id,
    function(input, output, session){
      data_deconvo <- rv_deconvo$data_deconvo

      if(!is.null(data_deconvo)){
        hnc_all <- data_deconvo$hnc_all

        output$title_csp <- renderUI({h4("Crime stain profile")})
        output$fn_csp <- renderText({paste0(data_deconvo$fn_csp)})
        output$title_af <- renderUI({h4("Allele frequency")})
        output$fn_af <- renderText({paste0(data_deconvo$fn_af)})
        output$title_gt_non_contrib <- renderUI({h4("Genotypes of non-contributor")})

        output$noc <- renderUI({
          selectInput(session$ns("noc"), label = h4("Number of contributors"), choices = hnc_all, selected = hnc_all[1], selectize = FALSE)
        })

        output$check_random <- renderUI({
          checkboxInput(session$ns("check_random"), label = "Generate randomly", value = FALSE)
        })

        output$file_gt_non_contrib <- renderUI({
          fileInput(session$ns("file_gt_non_contrib"), label = "", accept = ".csv", buttonLabel = "Load")
        })

        output$ui_n_simu <- renderUI({
          disabled(numericInput(session$ns("n_simu"), label = "Number of random genotypes", value = 100))
        })

        output$act_analyze_hd_true <- renderUI({
          disabled(actionButton(session$ns("act_analyze_hd_true"), label = "Analyze", class = "btn btn-primary btn-lg"))
        })

        observeEvent(input$check_random, {
          if(input$check_random){
            disable("file_gt_non_contrib")
            enable("n_simu")

            if(rv_hd_true_random$bool_n_simu){
              enable("act_analyze_hd_true")
            }else{
              disable("act_analyze_hd_true")
            }
          }else{
            enable("file_gt_non_contrib")
            disable("n_simu")

            if(!is.null(rv_hd_true_random$dt_gt_non_contrib)){
              enable("act_analyze_hd_true")
            }else{
              disable("act_analyze_hd_true")
            }
          }
        })

        observeEvent(input$file_gt_non_contrib, {
          file_input <- input$file_gt_non_contrib
          if(!is.null(file_input)){
            dt_gt_non_contrib <- fread(file_input$datapath)
            dt_gt_non_contrib <- try(change_data_type_ref(dt_gt_non_contrib))
            if(!inherits(dt_gt_non_contrib, "try-error")){
              rv_hd_true_random$dt_gt_non_contrib <- dt_gt_non_contrib
            }else{
              rv_hd_true_random$dt_gt_non_contrib <- NULL
              showModal(modalDialog(title = "Error", "The selected file format is invalid. Please upload a file in the correct format.", easyClose = TRUE, footer = NULL))
            }
          }else{
            rv_hd_true_random$dt_gt_non_contrib <- NULL
          }

          if(!is.null(rv_hd_true_random$dt_gt_non_contrib) && !input$check_random){
            enable("act_analyze_hd_true")
          }else{
            disable("act_analyze_hd_true")
          }
        })

        observeEvent(input$n_simu, {
          n_simu <- input$n_simu

          if(!is.integer(n_simu) || n_simu <= 0){
            showFeedbackDanger(inputId = "n_simu", text = "The value must be a positive integer.")
            rv_hd_true_random$bool_n_simu <- FALSE
          }else{
            hideFeedback("n_simu")
            rv_hd_true_random$bool_n_simu <- TRUE
          }

          if(rv_hd_true_random$bool_n_simu && input$check_random){
            enable("act_analyze_hd_true")
          }else{
            disable("act_analyze_hd_true")
          }
        })

        observeEvent(input$act_analyze_hd_true, {

          waiter_show(html = spin_3k(), color = "white")

          # --- Check input genotypes for non-contributors --- #

          if(!input$check_random){
            dt_gt_non_contrib <- rv_hd_true_random$dt_gt_non_contrib
            error_message <- check_gt_non_contrib(dt_gt_non_contrib, data_deconvo$markers_use)
            n_simu <- length(unique(dt_gt_non_contrib[, SampleName]))
          }else{
            dt_gt_non_contrib <- NULL
            error_message <- ""
            n_simu <- input$n_simu
          }

          if(error_message != ""){
            waiter_hide()
            showModal(modalDialog(title = "Error", HTML(error_message), easyClose = TRUE, footer = NULL))
          }else{
            hnc <- as.numeric(input$noc)
            tmp <- try(perform_hd_true(hnc, data_deconvo, n_simu, dt_gt_non_contrib, importance_sampling = FALSE))
            if(inherits(tmp, "try-error")){
              waiter_hide()
              showModal(modalDialog(title = "Error", "An unexpected error has occurred.", easyClose = TRUE, footer = NULL))
            }else{
              rv_hd_true_random$dt_lr_rand <- dt_lr_rand <- tmp$dt_lr_rand

              if(is.null(dt_lr_rand)){
                data_hnc_list <- data_deconvo$data_hnc_list
                data_h1 <- data_hnc_list[[which(hnc_all == hnc)]]
                error_message_h1 <- data_h1$error_message
                waiter_hide()
                showModal(modalDialog(title = "Error", HTML(error_message_h1), easyClose = TRUE, footer = NULL))
              }else{
                rv_hd_true_random$lr_rand <- as.numeric(dt_lr_rand[nrow(dt_lr_rand), .SD, .SDcols = !c("Marker")])
                rv_hd_true_random$dt_gt_non_contrib <- tmp$dt_gt_non_contrib
                waiter_hide()
              }
            }
          }
        })

        output$dt_summary_lr_rand <- renderDataTable({
          lr_rand <- rv_hd_true_random$lr_rand

          if(!is.null(lr_rand)){
            dt_summary_lr_rand <- data.table(Quantile = c("Maximum", "75%", "Median", "25%", "Minimum"),
                                             Value = c(max(lr_rand), quantile(lr_rand, 0.75), median(lr_rand), quantile(lr_rand, 0.25), min(lr_rand)))

            datatable(
              dt_summary_lr_rand,
              selection = "none",
              options = list(dom = "t", iDisplayLength = nrow(dt_summary_lr_rand), ordering = FALSE, autowidth = TRUE),
              rownames = FALSE
            ) %>%
              formatSignif(columns = c("Value"), digits = 3)
          }
        })

        output$ui_download_dt_lr_rand <- renderUI({
          if(!is.null(rv_hd_true_random$dt_lr_rand)){
            downloadButton(session$ns("download_dt_lr_rand"), "Download LRs")
          }
        })

        output$download_dt_lr_rand <- downloadHandler(
          filename = paste0(gsub(" ", "_", format(as.POSIXct(Sys.time()), "%Y-%m-%d %H%M%S")), "_Hd-true_LR.csv"),
          content = function(file){
            write.csv(rv_hd_true_random$dt_lr_rand, file, row.names = FALSE)
          }
        )

        output$ui_download_dt_gt_non_contrib <- renderUI({
          if(!is.null(rv_hd_true_random$dt_lr_rand)){
            downloadButton(session$ns("download_dt_gt_non_contrib"), "Download genotypes")
          }
        })

        output$download_dt_gt_non_contrib <- downloadHandler(
          filename = paste0(gsub(" ", "_", format(as.POSIXct(Sys.time()), "%Y-%m-%d %H%M%S")), "_Hd-true_genotypes.csv"),
          content = function(file){
            write.csv(rv_hd_true_random$dt_gt_non_contrib, file, row.names = FALSE)
          }
        )

        output$hist_lr_rand <- renderPlot({
          lr_rand <- rv_hd_true_random$lr_rand
          if(!is.null(lr_rand)){
            log10_lr_rand <- log10(lr_rand)
            log10_lr_rand[log10_lr_rand < -300] <- -300
            hist(log10_lr_rand)
          }
        })
      }

      return(rv_hd_true_random)
    }
  )
}
