#' lr_setting_ui
#'
#' @description The function to create the ui module for setting conditions of the likelihood ratio
lr_setting_ui <- function(id){
  ns <- NS(id)

  tabPanel("Setting",
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
                       tags$head( #CSS added to shrink the column with
                         tags$style('td.small .shiny-input-container{width:auto;}
                                     td.small{width:30px;}
                                     th { font-size: 16px !important; }'
                         )
                       ),
                       dataTableOutput(ns("dt_info_mk_display_lr"))
             )
           )
  )
}

#' lr_setting_server
#'
#' @description The function to create the server module for setting conditions of the likelihood ratio
lr_setting_server <- function(id, path_pack, rv_file, rv_deconvo, session_lr_top){
  moduleServer(
    id,
    function(input, output, session){
      shinyInput <- function(FUN, len, id, ...) {
        inputs <- character(len)
        for (i in seq_len(len)) {
          inputs[i] <- as.character(FUN(paste0(id, i), label = NULL, ...))
        }
        inputs
      }

      js <- function(dtid, ns){
        c(
          "$('[id^=checkb]').on('click', function(){",
          "  var id = this.getAttribute('id');",
          "  var i = parseInt(/checkb(\\d+)/.exec(id)[1]);",
          "  var value = $(this).prop('checked');",
          "  var info = [{row: i, col: 1, value: value}];",
          sprintf(
            "Shiny.setInputValue('%s', info);",
            ns(sprintf("%s_cell_edit:DT.cellInfo", dtid))
          ),
          "})"
        )
      }

      rv_lr_setting <- reactiveValues()
      rv_lr_setting$data_deconvo <- NULL
      rv_lr_setting$dt_csp <- NULL
      rv_lr_setting$dt_ref <- NULL
      rv_lr_setting$dt_af <- NULL
      rv_lr_setting$refs <- NULL
      rv_lr_setting$dt_anameth <- NULL
      rv_lr_setting$dt_info_mk <- NULL
      rv_lr_setting$warmup <- 1200
      rv_lr_setting$sampling <- 300
      rv_lr_setting$chains <- 4

      rv_lr <- reactiveValues()
      rv_lr$dt_lr <- NULL
      rv_lr$dt_csp <- NULL
      rv_lr$dt_ref <- NULL
      rv_lr$dt_af <- NULL
      rv_lr$anameth <- NULL
      rv_lr$dt_kit <- NULL
      rv_lr$dt_at <- NULL
      rv_lr$dt_param <- NULL
      rv_lr$dt_alcor <- NULL
      rv_lr$dt_gtc_fltr <- NULL
      rv_lr$dt_gtc_prob <- NULL
      rv_lr$warmup <- NULL
      rv_lr$sampling <- NULL
      rv_lr$chains <- NULL

      observe({
        req(rv_file)
        rv_lr_setting$dt_csp <- rv_file$dt_csp
        rv_lr_setting$dt_ref <- dt_ref <- rv_file$dt_ref
        rv_lr_setting$dt_af <- rv_file$dt_af

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
          updateSelectInput(session, "anameth_lr", selected = rv_deconvo$anameth)
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
        rv_lr_setting$dt_info_mk <- dt_info_mk

        # Create a displayed data.table for information on markers
        dt_info_mk_display_lr <- copy(dt_at)
        dt_info_mk_display_lr <- cbind(
          Ignore = shinyInput(checkboxInput, nrow(dt_info_mk_display_lr), "checkb"),
          dt_info_mk_display_lr
        )

        output$dt_info_mk_display_lr <- renderDataTable(server = FALSE, {
          datatable(
            dt_info_mk_display_lr,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_info_mk_display_lr), ordering = FALSE, autoWidth = TRUE, scrollX = TRUE, scrollY = "calc(100vh - 200px)",
                           columnDefs = list(list(targets = 0, className = "small"), list(targets = 1:2, width = "100px"))),
            rownames = FALSE,
            escape = FALSE,
            callback = JS(js("dt_info_mk_display_lr", session$ns))
          )
        })
      })

      observeEvent(input[["dt_info_mk_display_lr_cell_edit"]], {
        info <- input[["dt_info_mk_display_lr_cell_edit"]]
        dt_info_mk <- rv_lr_setting$dt_info_mk
        dt_info_mk[[info[1, "col"]]][info[1, "row"]] <- info[1, "value"]
        rv_lr_setting$dt_info_mk <- dt_info_mk
      })

      observeEvent(input$act_analyze_lr, {
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
            showModal(modalDialog(title = "Error", HTML(error_message), easyClose = TRUE, footer = NULL))
          }else{
            # Extract data from data_deconvo
            data_deconvo <- rv_lr_setting$data_deconvo
            hnc <- data_deconvo$hnc
            markers_use <- data_deconvo$markers_use
            pop_al_list <- data_deconvo$pop_al_list
            pop_prob_list <- data_deconvo$pop_prob_list
            theta <- data_deconvo$theta
            gtc_use <- data_deconvo$gtc_use
            lpd_gtc_use_1st <- data_deconvo$lpd_gtc_use_1st
            lpd_gtc_use_2nd <- data_deconvo$lpd_gtc_use_2nd
            dt_mpdeg <- data_deconvo$dt_mpdeg

            data_hp <- assign_log10like(hnc, markers_use, pop_al_list, pop_prob_list, theta, gtc_use, lpd_gtc_use_1st, lpd_gtc_use_2nd, dt_mpdeg, dt_ref, hp, refs)
            data_hd <- assign_log10like(hnc, markers_use, pop_al_list, pop_prob_list, theta, gtc_use, lpd_gtc_use_1st, lpd_gtc_use_2nd, dt_mpdeg, dt_ref, hd, refs)
            dt_lr <- make_dt_lr(markers_use,
                                data_hp$log10like_adopt, data_hp$combined_log10like_adopt,
                                data_hd$log10like_adopt, data_hd$combined_log10like_adopt)

            data_lr <- list("dt_mpdeg_hp" = data_hp$dt_mpdeg_hyp,
                            "dt_mpdeg_hd" = data_hd$dt_mpdeg_hyp,
                            "dt_pg_list_hp" = data_hp$dt_pg_list_hyp,
                            "dt_pg_list_hd" = data_hd$dt_pg_list_hyp,
                            "dt_lr" = dt_lr)
            rv_lr$data_lr <- data_lr

            updateTabsetPanel(session = session_lr_top, inputId = "tabset_lr", selected = "Result")
          }
        }else{
          # Fix setting
          hnc_lr_hp <- input$noc_lr_hp
          hnc_lr_hd <- input$noc_lr_hd
          warmup <- rv_lr_setting$warmup
          sampling <- rv_lr_setting$sampling
          chains <- rv_lr_setting$chains
          anameth <- input$anameth_lr

          load(paste0(path_pack, "/extdata/param/", dt_anameth$Value[dt_anameth$Category == "Parameter"], ".Rdata"))
          dt_al_prob <- data.table(Subcategory = dt_anameth$Subcategory[dt_anameth$Category == "Allele_probability"],
                                   Value = dt_anameth$Value[dt_anameth$Category == "Allele_probability"])
          dt_gt_fltr <- data.table(Subcategory = dt_anameth$Subcategory[dt_anameth$Category == "Genotype_filter"],
                                   Value = dt_anameth$Value[dt_anameth$Category == "Genotype_filter"])
          dt_dropin <- data.table(Subcategory = dt_anameth$Subcategory[dt_anameth$Category == "Drop-in"],
                                  Value = dt_anameth$Value[dt_anameth$Category == "Drop-in"])

          dt_info_mk <- rv_lr_setting$dt_info_mk

          # Check input files and setting
          error_message <- check_main_file(dt_csp, dt_ref, dt_af, markers_auto, method = "lr")
          if(error_message != ""){
            showModal(modalDialog(title = "Error", HTML(error_message), easyClose = TRUE, footer = NULL))
          }else{
            data_deconvo_hp <- perform_deconvo(dt_csp, dt_af, dt_kit_auto, dt_info_mk, dt_par, dt_sr, dt_al_prob, dt_gt_fltr,
                                               hnc_lr_hp, warmup, sampling, chains)
            if(hnc_lr_hp == hnc_lr_hd){
              data_deconvo_hd <- data_deconvo_hp
            }else{
              data_deconvo_hd <- perform_deconvo(dt_csp, dt_af, dt_kit_auto, dt_info_mk, dt_par, dt_sr, dt_al_prob, dt_gt_fltr,
                                                 hnc_lr_hd, warmup, sampling, chains)
            }

            # Extract data from data_deconvo for Hp
            data_deconvo <- rv_lr_setting$data_deconvo_hp
            hnc <- data_deconvo$hnc
            markers_use <- data_deconvo$markers_use
            pop_al_list <- data_deconvo$pop_al_list
            pop_prob_list <- data_deconvo$pop_prob_list
            theta <- data_deconvo$theta
            gtc_use <- data_deconvo$gtc_use
            lpd_gtc_use_1st <- data_deconvo$lpd_gtc_use_1st
            lpd_gtc_use_2nd <- data_deconvo$lpd_gtc_use_2nd
            dt_mpdeg <- data_deconvo$dt_mpdeg

            data_hp <- assign_log10like(hnc, markers_use, pop_al_list, pop_prob_list, theta, gtc_use, lpd_gtc_use_1st, lpd_gtc_use_2nd, dt_mpdeg, dt_ref, hp, refs)

            # Extract data from data_deconvo for Hd
            data_deconvo <- rv_lr_setting$data_deconvo_hd
            hnc <- data_deconvo$hnc
            markers_use <- data_deconvo$markers_use
            pop_al_list <- data_deconvo$pop_al_list
            pop_prob_list <- data_deconvo$pop_prob_list
            theta <- data_deconvo$theta
            gtc_use <- data_deconvo$gtc_use
            lpd_gtc_use_1st <- data_deconvo$lpd_gtc_use_1st
            lpd_gtc_use_2nd <- data_deconvo$lpd_gtc_use_2nd
            dt_mpdeg <- data_deconvo$dt_mpdeg

            data_hd <- assign_log10like(hnc, markers_use, pop_al_list, pop_prob_list, theta, gtc_use, lpd_gtc_use_1st, lpd_gtc_use_2nd, dt_mpdeg, dt_ref, hd, refs)

            dt_lr <- make_dt_lr(markers_use,
                                data_hp$log10like_adopt, data_hp$combined_log10like_adopt,
                                data_hd$log10like_adopt, data_hd$combined_log10like_adopt)
            rv_lr$data_lr <- data_lr

            updateTabsetPanel(session = session_lr_top, inputId = "tabset_lr", selected = "Result")
          }
        }
      })

      return(rv_lr)
    }
  )
}
