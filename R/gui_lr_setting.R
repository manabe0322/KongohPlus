#' lr_setting_ui
#'
#' @description The function to create the ui module for setting conditions of the likelihood ratio
lr_setting_ui <- function(id){
  ns <- NS(id)

  tabPanel("Setting",
           useWaiter(),
           tags$style(HTML("
                         .handsontable .htCore td {
                           vertical-align: middle !important;
                           text-align: center !important;
                         }"
           )),
           br(),
           fluidRow(
             column(3,
                    wellPanel(
                      uiOutput(ns("use_deconvo")),
                      uiOutput(ns("hp")),
                      numericInput(ns("noc_lr_hp"), label = "Number of contributors", value = 1),
                      uiOutput(ns("hd")),
                      numericInput(ns("noc_lr_hd"), label = "Number of contributors", value = 1),
                      br(),
                      h4("Analysis method"),
                      shinyFilesButton(ns("load_anameth_lr"), label = "Load", title = "Select a file", multiple = FALSE),
                      uiOutput(ns("anameth_lr")),
                      br(),
                      checkboxInput(ns("sub_source"), label = "Sub-source LR", value = TRUE),
                      br(),
                      br(),
                      actionButton(ns("act_analyze_lr"), label = "Analyze", class = "btn btn-primary btn-lg")
                    )
             ),
             column(4,
                    h4(textOutput(ns("title_ignore_mk"))),
                    rHandsontableOutput(ns("dt_ignore_mk"))
             ),
             column(5,
                    fluidRow(
                      column(4,
                             h4(textOutput(ns("title_at"))),
                             rHandsontableOutput(ns("anameth_at"))
                      ),
                      column(8,
                             h4(textOutput(ns("title_al_gt_prob"))),
                             rHandsontableOutput(ns("anameth_algt_prob")),
                             br(),
                             h4(textOutput(ns("title_gt_fltr"))),
                             rHandsontableOutput(ns("anameth_gt_fltr")),
                             br(),
                             h4(textOutput(ns("title_dropin"))),
                             rHandsontableOutput(ns("anameth_dropin")),
                             br(),
                             h4(textOutput(ns("title_mcmc"))),
                             rHandsontableOutput(ns("anameth_mcmc"))
                      )
                    )
             )
           )
  )
}

#' lr_setting_server
#'
#' @description The function to create the server module for setting conditions of the likelihood ratio
lr_setting_server <- function(id, path_pack, rv_file, rv_deconvo){
  moduleServer(
    id,
    function(input, output, session){
      rv_lr_setting <- reactiveValues()
      rv_lr_setting$data_deconvo <- NULL
      rv_lr_setting$dt_csp <- NULL
      rv_lr_setting$dt_ref <- NULL
      rv_lr_setting$dt_af <- NULL
      rv_lr_setting$fn_csp <- NULL
      rv_lr_setting$fn_ref <- NULL
      rv_lr_setting$fn_af <- NULL
      rv_lr_setting$refs <- NULL
      rv_lr_setting$fn_anameth_lr <- ""
      rv_lr_setting$dt_anameth <- NULL
      rv_lr_setting$dt_ignore_mk <- NULL
      rv_lr_setting$readonly_ignore_mk <- FALSE
      rv_lr_setting$warmup <- 1200
      rv_lr_setting$sampling <- 300
      rv_lr_setting$chains <- 4

      rv_lr <- reactiveValues()
      rv_lr$data_lr <- NULL
      rv_lr$data_deconvo_hp <- NULL
      rv_lr$data_deconvo_hd <- NULL

      observe({
        req(rv_file)
        rv_lr_setting$dt_csp <- rv_file$dt_csp
        rv_lr_setting$dt_ref <- dt_ref <- rv_file$dt_ref
        rv_lr_setting$dt_af <- rv_file$dt_af
        rv_lr_setting$fn_csp <- rv_file$fn_csp
        rv_lr_setting$fn_ref <- rv_file$fn_ref
        rv_lr_setting$fn_af <- rv_file$fn_af

        if(!is.null(dt_ref)){
          rv_lr_setting$refs <- unique(dt_ref[, SampleName])
        }
      })

      observe({
        req(rv_deconvo)
        rv_lr_setting$data_deconvo <- rv_deconvo$data_deconvo
        fn_anameth_lr <- rv_deconvo$data_deconvo$anameth
        if(!is.null(fn_anameth_lr)){
          dt_anameth <- fread(paste0(path_pack, "/extdata/anameth/", fn_anameth_lr))
          markers_auto <- dt_anameth$Subcategory[dt_anameth$Category == "AT"]
          markers_use <- rv_deconvo$data_deconvo$markers_use
          rv_lr_setting$dt_ignore_mk <- data.table(Ignore = !is.element(markers_auto, markers_use),
                                                   Marker = markers_auto)
          rv_lr_setting$readonly_ignore_mk <- TRUE
        }
      })

      output$use_deconvo <- renderUI({
        if(is.null(rv_lr_setting$data_deconvo)){
          disabled(checkboxInput(session$ns("use_deconvo"), label = "Use the result of deconvolution", value = FALSE))
        }else{
          checkboxInput(session$ns("use_deconvo"), label = "Use the result of deconvolution", value = TRUE)
        }
      })

      observeEvent(input$use_deconvo, {
        use_deconvo <- input$use_deconvo
        data_deconvo <- rv_lr_setting$data_deconvo
        if(use_deconvo){
          rv_lr_setting$fn_anameth_lr <- data_deconvo$anameth

          disable("noc_lr_hp")
          disable("noc_lr_hd")
          disable("load_anameth_lr")
          disable("act_advance_lr")
          updateNumericInput(session, "noc_lr_hp", value = data_deconvo$hnc)
          updateNumericInput(session, "noc_lr_hd", value = data_deconvo$hnc)
          updateTextInput(session, "anameth_lr", value = data_deconvo$anameth)

          dt_anameth <- fread(paste0(path_pack, "/extdata/anameth/", rv_lr_setting$fn_anameth_lr))
          markers_auto <- dt_anameth$Subcategory[dt_anameth$Category == "AT"]
          markers_use <- data_deconvo$markers_use

          rv_lr_setting$dt_ignore_mk <- data.table(Ignore = !is.element(markers_auto, markers_use),
                                                   Marker = markers_auto)
          rv_lr_setting$readonly_ignore_mk <- TRUE
        }else{
          enable("noc_lr_hp")
          enable("noc_lr_hd")
          enable("load_anameth_lr")
          enable("act_advance_lr")
          rv_lr_setting$readonly_ignore_mk <- FALSE
        }
      })

      observeEvent(input$noc_lr_hp, {
        noc_lr_hp <- input$noc_lr_hp
        if(!is.integer(noc_lr_hp) || noc_lr_hp < 1 || noc_lr_hp > 4){
          showFeedbackDanger(inputId = "noc_lr_hp", text = "Only 1, 2, 3, and 4 is allowed.")
          disable("act_analyze_lr")
        }else{
          hideFeedback("noc_lr_hp")
          enable("act_analyze_lr")
        }
      })

      observeEvent(input$noc_lr_hd, {
        noc_lr_hd <- input$noc_lr_hd
        if(!is.integer(noc_lr_hd) || noc_lr_hd < 1 || noc_lr_hd > 4){
          showFeedbackDanger(inputId = "noc_lr_hd", text = "Only 1, 2, 3, and 4 is allowed.")
          disable("act_analyze_lr")
        }else{
          hideFeedback("noc_lr_hd")
          enable("act_analyze_lr")
        }
      })

      output$hp <- renderUI({
        checkboxGroupInput(session$ns("hp"), label = h4("Prosecutor hypothesis"), choices = rv_lr_setting$refs)
      })

      output$hd <- renderUI({
        checkboxGroupInput(session$ns("hd"), label = h4("Defense hypothesis"), choices = rv_lr_setting$refs)
      })

      shinyFileChoose(input, "load_anameth_lr", roots = c(documents = paste0(path_pack, "/extdata/anameth")), session = session, filetypes = c("csv"))

      observeEvent(input$load_anameth_lr, {
        req(input$load_anameth_lr)
        path_anameth <- parseFilePaths(roots = c(documents = paste0(path_pack, "/extdata/anameth")), input$load_anameth_lr)$datapath
        if(length(path_anameth) > 0){
          rv_lr_setting$fn_anameth_lr <- basename(path_anameth)
        }
      })

      output$anameth_lr <- renderUI({
        disabled(textInput(session$ns("anameth_lr"), label = "", value = rv_lr_setting$fn_anameth_lr))
      })

      observeEvent(input$anameth_lr, {
        fn_anameth_lr <- input$anameth_lr
        if(fn_anameth_lr != ""){
          rv_lr_setting$dt_anameth <- dt_anameth <- fread(paste0(path_pack, "/extdata/anameth/", fn_anameth_lr))
          if(!input$use_deconvo){
            markers_auto <- dt_anameth$Subcategory[dt_anameth$Category == "AT"]
            rv_lr_setting$dt_ignore_mk <- dt_ignore_mk <- data.table(Ignore = rep(FALSE, length(markers_auto)),
                                                                     Marker = markers_auto)
          }

          output$title_ignore_mk <- renderText({"Ignored marker"})
          output$title_at <- renderText({"Analytical threshold"})
          output$title_al_gt_prob <- renderText({"Genotype probability"})
          output$title_gt_fltr <- renderText({"Genotype filter"})
          output$title_dropin <- renderText({"Drop-in"})
          output$title_mcmc <- renderText({"MCMC"})

          output$anameth_at <- renderRHandsontable({
            dt_at <- data.table(Marker = dt_anameth$Subcategory[dt_anameth$Category == "AT"],
                                AT = dt_anameth$Value[dt_anameth$Category == "AT"])
            rhandsontable(dt_at, rowHeaders = FALSE) %>%
              hot_table(fillHandle = FALSE, rowHeights = 30) %>%
              hot_col("Marker", type = "text", readOnly = TRUE, width = 100) %>%
              hot_col("AT", type = "numeric", format = "0", readOnly = TRUE, width = 100)
          })

          output$anameth_algt_prob <- renderRHandsontable({
            dt_algt_prob <- data.table(Subcategory = dt_anameth$Subcategory[dt_anameth$Category == "Genotype probability"],
                                       Value = dt_anameth$Value[dt_anameth$Category == "Genotype probability"])
            rhandsontable(dt_algt_prob, rowHeaders = FALSE, colHeaders = FALSE) %>%
              hot_table(fillHandle = FALSE, rowHeights = 30) %>%
              hot_col(col = 1, type = "text", readOnly = TRUE, width = 300) %>%
              hot_col(col = 2, type = "numeric", readOnly = TRUE, width = 100)
          })

          output$anameth_gt_fltr <- renderRHandsontable({
            dt_gt_fltr <- data.table(Subcategory = dt_anameth$Subcategory[dt_anameth$Category == "Genotype filter"],
                                     Value = dt_anameth$Value[dt_anameth$Category == "Genotype filter"])
            rhandsontable(dt_gt_fltr , rowHeaders = FALSE, colHeaders = FALSE) %>%
              hot_table(fillHandle = FALSE, rowHeights = 30) %>%
              hot_col(col = 1, type = "text", readOnly = TRUE, width = 300) %>%
              hot_col(col = 2, type = "numeric", readOnly = TRUE, width = 100)
          })

          output$anameth_dropin <- renderRHandsontable({
            dt_dropin <- data.table(Subcategory = dt_anameth$Subcategory[dt_anameth$Category == "Drop-in"],
                                    Value = dt_anameth$Value[dt_anameth$Category == "Drop-in"])
            rhandsontable(dt_dropin, rowHeaders = FALSE, colHeaders = FALSE) %>%
              hot_table(fillHandle = FALSE, rowHeights = 30) %>%
              hot_col(col = 1, type = "text", readOnly = TRUE, width = 300) %>%
              hot_col(col = 2, type = "numeric", readOnly = TRUE, width = 100)
          })

          output$anameth_mcmc <- renderRHandsontable({
            dt_mcmc <- data.table(Subcategory = dt_anameth$Subcategory[dt_anameth$Category == "MCMC"],
                                  Value = dt_anameth$Value[dt_anameth$Category == "MCMC"])
            rhandsontable(dt_mcmc, rowHeaders = FALSE, colHeaders = FALSE) %>%
              hot_table(fillHandle = FALSE, rowHeights = 30) %>%
              hot_col(col = 1, type = "text", readOnly = TRUE, width = 300) %>%
              hot_col(col = 2, type = "numeric", readOnly = TRUE, width = 100)
          })
        }
      })

      output$dt_ignore_mk <- renderRHandsontable({
        dt_ignore_mk <- rv_lr_setting$dt_ignore_mk
        if(!is.null(dt_ignore_mk)){
          rhandsontable(dt_ignore_mk, rowHeaders = FALSE) %>%
            hot_table(fillHandle = FALSE, rowHeights = 30) %>%
            hot_col("Ignore", type = "checkbox", readOnly = rv_lr_setting$readonly_ignore_mk, width = 100) %>%
            hot_col("Marker", type = "text", readOnly = TRUE, width = 100)
        }
      })

      observeEvent(input$act_analyze_lr, {
        waiter_show(html = spin_3k(), color = "white")

        # Fix setting
        hp <- input$hp
        hd <- input$hd

        refs <- rv_lr_setting$refs
        dt_csp <- rv_lr_setting$dt_csp
        dt_ref <- rv_lr_setting$dt_ref
        dt_af <- rv_lr_setting$dt_af
        fn_csp <- rv_lr_setting$fn_csp
        fn_ref <- rv_lr_setting$fn_ref
        fn_af <- rv_lr_setting$fn_af
        anameth <- rv_lr_setting$fn_anameth_lr
        dt_anameth <- rv_lr_setting$dt_anameth
        use_deconvo <- input$use_deconvo
        sub_source <- input$sub_source

        # Extract information on the kit
        dt_kit <- fread(paste0(path_pack, "/extdata/kit/", dt_anameth$Name[dt_anameth$Category == "Kit"], ".csv"))
        dt_kit <- change_data_type_kit(dt_kit)
        dt_kit_auto <- extract_kit_auto(dt_kit)
        markers_auto <- unique(dt_kit_auto[, Marker])

        if(use_deconvo){
          error_message <- check_main_file(dt_csp, dt_ref, dt_af, markers_auto, method = "lr")
          if(error_message != ""){
            waiter_hide()
            showModal(modalDialog(title = "Error", HTML(error_message), easyClose = TRUE, footer = NULL))
          }else{
            # Extract data from data_deconvo
            data_deconvo <- rv_lr_setting$data_deconvo
            hnc <- data_deconvo$hnc
            markers_use <- data_deconvo$markers_use
            peaks_list <- data_deconvo$peaks_list
            pop_al_list <- data_deconvo$pop_al_list
            pop_prob_list <- data_deconvo$pop_prob_list
            theta <- as.numeric(dt_anameth$Value[dt_anameth$Subcategory == "Theta"])
            gtc_use <- data_deconvo$gtc_use
            lpd_gtc_use_1st <- data_deconvo$lpd_gtc_use_1st
            lpd_gtc_use_2nd <- data_deconvo$lpd_gtc_use_2nd

            log10like_hp <- assign_log10like(hnc, markers_use, peaks_list, pop_al_list, pop_prob_list, theta, gtc_use, lpd_gtc_use_1st, lpd_gtc_use_2nd, dt_ref, hp, refs, sub_source)
            log10like_hd <- assign_log10like(hnc, markers_use, peaks_list, pop_al_list, pop_prob_list, theta, gtc_use, lpd_gtc_use_1st, lpd_gtc_use_2nd, dt_ref, hd, refs, sub_source)
            dt_lr <- data.table(Marker = c(markers_use, "Total"),
                                Log10_Likelihood_Hp = log10like_hp,
                                Log10_Likelihood_Hd = log10like_hd,
                                LR = 10^(log10like_hp - log10like_hd))
            names(dt_lr) <- c("Marker", "log10 (f(E | Hp))", "log10 (f(E | Hd))", "LR")

            data_lr <- list("dt_ref" = dt_ref,
                            "fn_ref" = rv_lr_setting$fn_ref,
                            "dt_lr" = dt_lr,
                            "refs" = refs,
                            "hp" = hp,
                            "hd" = hd)
            rv_lr$data_lr <- data_lr
            rv_lr$data_deconvo_hp <- data_deconvo
            rv_lr$data_deconvo_hd <- data_deconvo

            waiter_hide()
          }
        }else{
          # Fix setting
          hnc_lr_hp <- input$noc_lr_hp
          hnc_lr_hd <- input$noc_lr_hd

          load(paste0(path_pack, "/extdata/param/", dt_anameth$Name[dt_anameth$Category == "Parameter"], ".RData"))
          load(paste0(path_pack, "/extdata/stutter_file/", dt_anameth$Name[dt_anameth$Category == "Stutter file"], ".RData"))

          dt_at <- data.table(Marker = dt_anameth$Subcategory[dt_anameth$Category == "AT"],
                              AT = dt_anameth$Value[dt_anameth$Category == "AT"])

          dt_al_prob <- data.table(Subcategory = dt_anameth$Subcategory[dt_anameth$Category == "Genotype probability"],
                                   Value = dt_anameth$Value[dt_anameth$Category == "Genotype probability"])
          dt_gt_fltr <- data.table(Subcategory = dt_anameth$Subcategory[dt_anameth$Category == "Genotype filter"],
                                   Value = dt_anameth$Value[dt_anameth$Category == "Genotype filter"])
          dt_dropin <- data.table(Subcategory = dt_anameth$Subcategory[dt_anameth$Category == "Drop-in"],
                                  Value = dt_anameth$Value[dt_anameth$Category == "Drop-in"])
          dt_mcmc <- data.table(Subcategory = dt_anameth$Subcategory[dt_anameth$Category == "MCMC"],
                                Value = dt_anameth$Value[dt_anameth$Category == "MCMC"])

          dt_ignore_mk <- hot_to_r(input$dt_ignore_mk)

          # Check input files and setting
          error_message <- check_main_file(dt_csp, dt_ref, dt_af, markers_auto, method = "lr")
          if(error_message != ""){
            waiter_hide()
            showModal(modalDialog(title = "Error", HTML(error_message), easyClose = TRUE, footer = NULL))
          }else{
            # Check registration of the stutter ratio
            message_sr <- check_csp_sr(dt_csp, dt_ignore_mk, dt_at, list_dt_average_sr)
            if(message_sr == "error"){
              waiter_hide()
              showModal(modalDialog(title = "Error", "All markers are ignored!", easyClose = TRUE, footer = NULL))
            }else{
              data_deconvo_hp <- perform_deconvo(dt_csp, dt_af, dt_kit, dt_ignore_mk, dt_at, dt_al_prob, dt_gt_fltr, dt_dropin, dt_mcmc,
                                                 median_alpha, median_cv, list_dt_average_sr,
                                                 hnc_lr_hp)
              if(!is.null(data_deconvo_hp$error_message)){
                waiter_hide()
                showModal(modalDialog(title = "Error", HTML(data_deconvo_hp$error_message), easyClose = TRUE, footer = NULL))
              }else{
                if(hnc_lr_hp == hnc_lr_hd){
                  data_deconvo_hd <- data_deconvo_hp
                }else{
                  data_deconvo_hd <- perform_deconvo(dt_csp, dt_af, dt_kit, dt_ignore_mk, dt_at, dt_al_prob, dt_gt_fltr, dt_dropin, dt_mcmc,
                                                     median_alpha, median_cv, list_dt_average_sr,
                                                     hnc_lr_hd)
                }
                if(!is.null(data_deconvo_hd$error_message)){
                  waiter_hide()
                  showModal(modalDialog(title = "Error", HTML(data_deconvo_hd$error_message), easyClose = TRUE, footer = NULL))
                }else{
                  data_deconvo_hp$fn_csp <- fn_csp
                  data_deconvo_hp$dt_csp <- dt_csp
                  data_deconvo_hp$fn_af <- fn_af
                  data_deconvo_hp$dt_af <- dt_af
                  data_deconvo_hp$anameth <- anameth
                  data_deconvo_hp$dt_anameth <- dt_anameth
                  data_deconvo_hp$dt_kit <- dt_kit
                  data_deconvo_hp$dt_ignore_mk <- dt_ignore_mk
                  data_deconvo_hp$alpha_auto <- median_alpha
                  data_deconvo_hp$cv <- median_cv
                  data_deconvo_hp$list_dt_average_sr <- list_dt_average_sr
                  data_deconvo_hp$hnc <- hnc_lr_hp
                  data_deconvo_hp$message_sr <- message_sr

                  data_deconvo_hd$fn_csp <- fn_csp
                  data_deconvo_hd$dt_csp <- dt_csp
                  data_deconvo_hd$fn_af <- fn_af
                  data_deconvo_hd$dt_af <- dt_af
                  data_deconvo_hd$anameth <- anameth
                  data_deconvo_hd$dt_anameth <- dt_anameth
                  data_deconvo_hd$dt_kit <- dt_kit
                  data_deconvo_hd$dt_ignore_mk <- dt_ignore_mk
                  data_deconvo_hd$alpha_auto <- median_alpha
                  data_deconvo_hd$cv <- median_cv
                  data_deconvo_hd$list_dt_average_sr <- list_dt_average_sr
                  data_deconvo_hd$hnc <- hnc_lr_hd
                  data_deconvo_hd$message_sr <- message_sr

                  # Extract data from data_deconvo for Hp
                  hnc <- data_deconvo_hp$hnc
                  markers_use <- data_deconvo_hp$markers_use
                  peaks_list <- data_deconvo_hp$peaks_list
                  pop_al_list <- data_deconvo_hp$pop_al_list
                  pop_prob_list <- data_deconvo_hp$pop_prob_list
                  gtc_use <- data_deconvo_hp$gtc_use
                  lpd_gtc_use_1st <- data_deconvo_hp$lpd_gtc_use_1st
                  lpd_gtc_use_2nd <- data_deconvo_hp$lpd_gtc_use_2nd

                  theta <- as.numeric(dt_al_prob$Value[dt_al_prob$Subcategory == "Theta"])

                  log10like_hp <- assign_log10like(hnc_lr_hp, markers_use, peaks_list, pop_al_list, pop_prob_list, theta, gtc_use, lpd_gtc_use_1st, lpd_gtc_use_2nd, dt_ref, hp, refs, sub_source)

                  # Extract data from data_deconvo for Hd
                  hnc <- data_deconvo_hd$hnc
                  markers_use <- data_deconvo_hd$markers_use
                  peaks_list <- data_deconvo_hd$peaks_list
                  pop_al_list <- data_deconvo_hd$pop_al_list
                  pop_prob_list <- data_deconvo_hd$pop_prob_list
                  gtc_use <- data_deconvo_hd$gtc_use
                  lpd_gtc_use_1st <- data_deconvo_hd$lpd_gtc_use_1st
                  lpd_gtc_use_2nd <- data_deconvo_hd$lpd_gtc_use_2nd

                  log10like_hd <- assign_log10like(hnc_lr_hd, markers_use, peaks_list, pop_al_list, pop_prob_list, theta, gtc_use, lpd_gtc_use_1st, lpd_gtc_use_2nd, dt_ref, hd, refs, sub_source)

                  dt_lr <- data.table(Marker = c(markers_use, "Total"),
                                      Log10_Likelihood_Hp = log10like_hp,
                                      Log10_Likelihood_Hd = log10like_hd,
                                      LR = 10^(log10like_hp - log10like_hd))
                  names(dt_lr) <- c("Marker", "log10 (f(E | Hp))", "log10 (f(E | Hd))", "LR")

                  data_lr <- list("dt_ref" = dt_ref,
                                  "fn_ref" = fn_ref,
                                  "dt_lr" = dt_lr,
                                  "refs" = refs,
                                  "hp" = hp,
                                  "hd" = hd)
                  rv_lr$data_lr <- data_lr
                  rv_lr$data_deconvo_hp <- data_deconvo_hp
                  rv_lr$data_deconvo_hd <- data_deconvo_hd

                  waiter_hide()
                }
              }
            }
          }
        }
      })

      return(rv_lr)
    }
  )
}
