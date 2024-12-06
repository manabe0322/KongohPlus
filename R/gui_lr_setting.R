#' lr_setting_ui
#'
#' @description The function to create the ui module for setting conditions of the likelihood ratio
lr_setting_ui <- function(id){
  ns <- NS(id)

  tabPanel("Setting",
           useWaiter(),
           br(),
           sidebarLayout(
             sidebarPanel(width = 3,
                          uiOutput(ns("use_deconvo")),
                          uiOutput(ns("hp")),
                          numericInput(ns("noc_lr_hp"), label = "Number of contributors", value = 1),
                          uiOutput(ns("hd")),
                          numericInput(ns("noc_lr_hd"), label = "Number of contributors", value = 1),
                          br(),
                          uiOutput(ns("anameth_lr")),
                          br(),
                          actionButton(ns("act_advance_lr"), label = "Advanced setting"),
                          br(),
                          br(),
                          actionButton(ns("act_analyze_lr"), label = "Analyze", class = "btn btn-primary btn-lg")
             ),
             mainPanel(width = 9,
                       tags$style(HTML("
                         .handsontable .htCore td {
                           vertical-align: middle !important;
                           text-align: center !important;
                         }"
                       )),
                       rHandsontableOutput(ns("dt_info_mk"))
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
      rv_lr_setting$dt_anameth <- NULL
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
          disable("noc_lr_hp")
          disable("noc_lr_hd")
          disable("anameth_lr")
          disable("act_advance_lr")
          updateNumericInput(session, "noc_lr_hp", value = data_deconvo$hnc)
          updateNumericInput(session, "noc_lr_hd", value = data_deconvo$hnc)
          updateSelectInput(session, "anameth_lr", selected = data_deconvo$anameth)
        }else{
          enable("noc_lr_hp")
          enable("noc_lr_hd")
          enable("anameth_lr")
          enable("act_advance_lr")
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

      observeEvent(input$act_advance_lr, {
        showModal(modalDialog(
          title = "Advanced setting",
          fluidRow(
            column(6,
                   h4("MCMC conditions"),
                   numericInput(session$ns("warmup"), label = "Warmup iterations", value = rv_lr_setting$warmup),
                   numericInput(session$ns("sampling"), label = "Sampling iterations", value = rv_lr_setting$sampling),
                   numericInput(session$ns("chains"), label = "Number of chains", value = rv_lr_setting$chains)
            ),
            column(6)
          ),
          footer = tagList(
            actionButton(session$ns("act_advance_save"), "Save"),
            modalButton("Cancel")
          ),
          size = "l"
        ))
      })

      observeEvent(input$warmup, {
        warmup <- input$warmup
        if(!is.integer(warmup)){
          showFeedbackDanger(inputId = "warmup", text = "A positive integer is required!")
        }else if(warmup <= 0){
          showFeedbackDanger(inputId = "warmup", text = "A positive integer is required!")
        }else{
          hideFeedback("warmup")
        }
      })

      observeEvent(input$sampling, {
        sampling <- input$sampling
        if(!is.integer(sampling)){
          showFeedbackDanger(inputId = "sampling", text = "A positive integer is required!")
        }else if(sampling <= 0){
          showFeedbackDanger(inputId = "sampling", text = "A positive integer is required!")
        }else{
          hideFeedback("sampling")
        }
      })

      observeEvent(input$chains, {
        chains <- input$chains
        if(!is.integer(chains)){
          showFeedbackDanger(inputId = "chains", text = "A positive integer is required!")
        }else if(chains <= 0){
          showFeedbackDanger(inputId = "chains", text = "A positive integer is required!")
        }else{
          hideFeedback("chains")
        }
      })

      observeEvent(input$act_advance_save, {
        req(input$warmup)
        req(input$sampling)
        req(input$chains)

        warmup <- input$warmup
        sampling <- input$sampling
        chains <- input$chains

        if(all(is.integer(warmup), is.integer(sampling), is.integer(chains))){
          rv_lr_setting$warmup <- warmup
          rv_lr_setting$sampling <- sampling
          rv_lr_setting$chains <- chains
          removeModal()
        }
      })

      output$anameth_lr <- renderUI({
        anameth_all <- list.files(paste0(path_pack, "/extdata/anameth"))
        anameth_all <- gsub(".{4}$", "", anameth_all)
        selectInput(session$ns("anameth_lr"), label = h4("Analysis method"), choices = anameth_all, selected = anameth_all[1], selectize = FALSE)
      })

      observeEvent(input$anameth_lr, {
        # Load the selected analysis method
        rv_lr_setting$dt_anameth <- dt_anameth <- fread(paste0(path_pack, "/extdata/anameth/", input$anameth_lr, ".csv"))

        # Create a data.table for information on markers
        dt_at <- data.table(Marker = dt_anameth$Subcategory[dt_anameth$Category == "AT"],
                            AT = as.numeric(dt_anameth$Value[dt_anameth$Category == "AT"]))
        dt_info_mk <- copy(dt_at)
        dt_info_mk <- cbind(
          Ignore = rep(FALSE, nrow(dt_info_mk)),
          dt_info_mk
        )

        output$dt_info_mk <- renderRHandsontable({
          rhandsontable(dt_info_mk, rowHeaders = FALSE) %>%
            hot_table(fillHandle = FALSE, rowHeights = 30) %>%
            hot_col("Ignore", type = "checkbox") %>%
            hot_col("Marker", type = "text", readOnly = TRUE) %>%
            hot_col("AT", type = "numeric", format = "0", readOnly = TRUE)
        })
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
        dt_anameth <- rv_lr_setting$dt_anameth
        use_deconvo <- input$use_deconvo

        # Extract information on the kit
        dt_kit <- fread(paste0(path_pack, "/extdata/kit/", dt_anameth$Value[dt_anameth$Category == "Kit"], ".csv"))
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
            theta <- data_deconvo$theta
            gtc_use <- data_deconvo$gtc_use
            lpd_gtc_use_1st <- data_deconvo$lpd_gtc_use_1st
            lpd_gtc_use_2nd <- data_deconvo$lpd_gtc_use_2nd

            log10like_hp <- assign_log10like(hnc, markers_use, peaks_list, pop_al_list, pop_prob_list, theta, gtc_use, lpd_gtc_use_1st, lpd_gtc_use_2nd, dt_ref, hp, refs)
            log10like_hd <- assign_log10like(hnc, markers_use, peaks_list, pop_al_list, pop_prob_list, theta, gtc_use, lpd_gtc_use_1st, lpd_gtc_use_2nd, dt_ref, hd, refs)
            dt_lr <- data.table(Marker = c(markers_use, "Total"),
                                Log10_Likelihood_Hp = log10like_hp,
                                Log10_Likelihood_Hd = log10like_hd,
                                LR = 10^(log10like_hp - log10like_hd))
            names(dt_lr) <- c("Marker", "log10 (f(E | Hp))", "log10 (f(E | Hd))", "LR")

            data_lr <- list("dt_ref" = dt_ref,
                            "fn_ref" = rv_lr_setting$fn_ref,
                            "dt_lr" = dt_lr)
            rv_lr$data_lr <- data_lr
            rv_lr$data_deconvo_hp <- data_deconvo
            rv_lr$data_deconvo_hd <- data_deconvo

            waiter_hide()
          }
        }else{
          # Fix setting
          hnc_lr_hp <- input$noc_lr_hp
          hnc_lr_hd <- input$noc_lr_hd
          warmup <- rv_lr_setting$warmup
          sampling <- rv_lr_setting$sampling
          chains <- rv_lr_setting$chains
          anameth <- input$anameth_lr

          load(paste0(path_pack, "/extdata/param/", dt_anameth$Value[dt_anameth$Category == "Parameter"], ".RData"))
          load(paste0(path_pack, "/extdata/stutter_file/", dt_anameth$Value[dt_anameth$Category == "Stutter file"], ".RData"))

          dt_al_prob <- data.table(Subcategory = dt_anameth$Subcategory[dt_anameth$Category == "Allele_probability"],
                                   Value = dt_anameth$Value[dt_anameth$Category == "Allele_probability"])
          dt_gt_fltr <- data.table(Subcategory = dt_anameth$Subcategory[dt_anameth$Category == "Genotype_filter"],
                                   Value = dt_anameth$Value[dt_anameth$Category == "Genotype_filter"])
          dt_dropin <- data.table(Subcategory = dt_anameth$Subcategory[dt_anameth$Category == "Drop-in"],
                                  Value = dt_anameth$Value[dt_anameth$Category == "Drop-in"])

          dt_info_mk <- hot_to_r(input$dt_info_mk)

          # Check input files and setting
          error_message <- check_main_file(dt_csp, dt_ref, dt_af, markers_auto, method = "lr")
          if(error_message != ""){
            waiter_hide()
            showModal(modalDialog(title = "Error", HTML(error_message), easyClose = TRUE, footer = NULL))
          }else{
            # Check registration of the stutter ratio
            message_sr <- check_csp_sr(dt_csp, dt_info_mk, list_dt_average_sr)
            if(message_sr == "error"){
              waiter_hide()
              showModal(modalDialog(title = "Error", "All markers are ignored!", easyClose = TRUE, footer = NULL))
            }else{
              data_deconvo_hp <- perform_deconvo(dt_csp, dt_af, dt_kit, dt_info_mk, dt_al_prob, dt_gt_fltr, dt_dropin,
                                                 alpha_auto, cv, list_dt_average_sr,
                                                 hnc_lr_hp, warmup, sampling, chains)
              if(!is.null(data_deconvo_hp$error_message)){
                waiter_hide()
                showModal(modalDialog(title = "Error", HTML(data_deconvo_hp$error_message), easyClose = TRUE, footer = NULL))
              }else{
                if(hnc_lr_hp == hnc_lr_hd){
                  data_deconvo_hd <- data_deconvo_hp
                }else{
                  data_deconvo_hd <- perform_deconvo(dt_csp, dt_af, dt_kit, dt_info_mk, dt_al_prob, dt_gt_fltr, dt_dropin,
                                                     alpha_auto, cv, list_dt_average_sr,
                                                     hnc_lr_hd, warmup, sampling, chains)
                }
                if(!is.null(data_deconvo_hd$error_message)){
                  waiter_hide()
                  showModal(modalDialog(title = "Error", HTML(data_deconvo_hd$error_message), easyClose = TRUE, footer = NULL))
                }else{
                  data_deconvo_hp$message_sr <- message_sr
                  data_deconvo_hp$fn_csp <- rv_lr_setting$fn_csp
                  data_deconvo_hp$fn_af <- rv_lr_setting$fn_af
                  data_deconvo_hp$anameth <- anameth
                  data_deconvo_hd$message_sr <- message_sr
                  data_deconvo_hd$fn_csp <- rv_lr_setting$fn_csp
                  data_deconvo_hd$fn_af <- rv_lr_setting$fn_af
                  data_deconvo_hd$anameth <- anameth

                  # Extract data from data_deconvo for Hp
                  hnc <- data_deconvo_hp$hnc
                  markers_use <- data_deconvo_hp$markers_use
                  peaks_list <- data_deconvo_hp$peaks_list
                  pop_al_list <- data_deconvo_hp$pop_al_list
                  pop_prob_list <- data_deconvo_hp$pop_prob_list
                  theta <- data_deconvo_hp$theta
                  gtc_use <- data_deconvo_hp$gtc_use
                  lpd_gtc_use_1st <- data_deconvo_hp$lpd_gtc_use_1st
                  lpd_gtc_use_2nd <- data_deconvo_hp$lpd_gtc_use_2nd

                  log10like_hp <- assign_log10like(hnc, markers_use, peaks_list, pop_al_list, pop_prob_list, theta, gtc_use, lpd_gtc_use_1st, lpd_gtc_use_2nd, dt_ref, hp, refs)

                  # Extract data from data_deconvo for Hd
                  hnc <- data_deconvo_hd$hnc
                  markers_use <- data_deconvo_hd$markers_use
                  peaks_list <- data_deconvo_hd$peaks_list
                  pop_al_list <- data_deconvo_hd$pop_al_list
                  pop_prob_list <- data_deconvo_hd$pop_prob_list
                  theta <- data_deconvo_hd$theta
                  gtc_use <- data_deconvo_hd$gtc_use
                  lpd_gtc_use_1st <- data_deconvo_hd$lpd_gtc_use_1st
                  lpd_gtc_use_2nd <- data_deconvo_hd$lpd_gtc_use_2nd

                  log10like_hd <- assign_log10like(hnc, markers_use, peaks_list, pop_al_list, pop_prob_list, theta, gtc_use, lpd_gtc_use_1st, lpd_gtc_use_2nd, dt_ref, hd, refs)

                  dt_lr <- data.table(Marker = c(markers_use, "Total"),
                                      Log10_Likelihood_Hp = log10like_hp,
                                      Log10_Likelihood_Hd = log10like_hd,
                                      LR = 10^(log10like_hp - log10like_hd))
                  names(dt_lr) <- c("Marker", "log10 (f(E | Hp))", "log10 (f(E | Hd))", "LR")

                  data_lr <- list("dt_ref" = dt_ref,
                                  "fn_ref" = rv_lr_setting$fn_ref,
                                  "dt_lr" = dt_lr)
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
