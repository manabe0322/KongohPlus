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
                       fluidRow(
                         column(4,
                                dataTableOutput(ns("dt_at")),
                                br()
                         ),
                         column(8,
                                h4("Kit"),
                                textOutput(ns("kit")),
                                br(),
                                h4("Allele probability"),
                                textOutput(ns("alprob")),
                                br(),
                                h4("Theta"),
                                textOutput(ns("theta")),
                                br(),
                                h4("Parameters"),
                                actionButton(ns("act_ae"), label = "AE"),
                                actionButton(ns("act_bsr"), label = "BSR"),
                                actionButton(ns("act_fsr"), label = "FSR"),
                                actionButton(ns("act_dsr"), label = "DSR"),
                                actionButton(ns("act_m2sr"), label = "M2SR"),
                                actionButton(ns("act_gamma"), label = "Gamma"),
                                br(),
                                br(),
                                h4("Conditions to remove unrealistic genotype combinations"),
                                textOutput(ns("bsr_fltr")),
                                textOutput(ns("fsr_fltr")),
                                textOutput(ns("dsr_fltr")),
                                textOutput(ns("m2sr_fltr"))
                         )
                       )
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
      rv_lr_setting <- reactiveValues()
      rv_lr_setting$data_deconvo <- NULL
      rv_lr_setting$dt_csp <- NULL
      rv_lr_setting$dt_ref <- NULL
      rv_lr_setting$dt_af <- NULL
      rv_lr_setting$refs <- NULL
      rv_lr_setting$anameth <- NULL
      rv_lr_setting$dt_anameth <- NULL
      rv_lr_setting$dt_kit <- NULL
      rv_lr_setting$dt_at <- NULL
      rv_lr_setting$dt_param <- NULL
      rv_lr_setting$dt_alcor <- NULL
      rv_lr_setting$gtc_fltr <- NULL
      rv_lr_setting$gtc_prob <- NULL
      rv_lr_setting$dt_ae <- NULL
      rv_lr_setting$dt_bsr <- NULL
      rv_lr_setting$dt_fsr <- NULL
      rv_lr_setting$dt_dsr <- NULL
      rv_lr_setting$dt_m2sr <- NULL
      rv_lr_setting$dt_gamma <- NULL
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

      anameth <- list.files(paste0(path_pack, "/extdata/anameth"))

      output$anameth_lr <- renderUI({
        anameth_all <- list.files(paste0(path_pack, "/extdata/anameth"))
        anameth_all <- gsub(".{4}$", "", anameth_all)
        selectInput(session$ns("anameth_lr"), label = h4("Analysis method"), choices = anameth_all, selected = anameth_all[1], selectize = FALSE)
      })

      observeEvent(input$anameth_lr, {
        # Load the selected analysis method
        anameth <- input$anameth_lr
        rv_lr_setting$anameth <- anameth
        rv_lr_setting$dt_anameth <- dt_anameth <- fread(paste0(path_pack, "/extdata/anameth/", anameth, ".csv"))

        # Extract data per column
        folders <- dt_anameth[, Folder]
        names <- dt_anameth[, Name]

        # Load each data.table
        rv_lr_setting$dt_kit <- fread(paste0(path_pack, "/extdata/kit/", names[folders == "kit"], ".csv"))
        rv_lr_setting$dt_at <- dt_at <- fread(paste0(path_pack, "/extdata/at/", names[folders == "at"], ".csv"))
        rv_lr_setting$dt_param <- fread(paste0(path_pack, "/extdata/param/", names[folders == "param"], ".csv"))
        rv_lr_setting$dt_alcor <- fread(paste0(path_pack, "/extdata/alcor/", names[folders == "alcor"], ".csv"))
        rv_lr_setting$dt_gtc_fltr <- dt_gtc_fltr <- fread(paste0(path_pack, "/extdata/gtc_fltr/", names[folders == "gtc_fltr"], ".csv"))
        rv_lr_setting$dt_gtc_prob <- dt_gtc_prob <- fread(paste0(path_pack, "/extdata/gtc_prob/", names[folders == "gtc_prob"], ".csv"))

        # Display analytical thresholds
        output$dt_at <- renderDataTable(server = FALSE, {
          datatable(
            dt_at,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_at), ordering = FALSE, autowidth = TRUE),
            rownames = FALSE
          )
        })

        gtc_prob_param <- dt_gtc_prob[, Parameter]
        gtc_prob_value <- dt_gtc_prob[, Value]

        gtc_fltr_type <- dt_gtc_fltr[, Type]
        gtc_fltr_value <- dt_gtc_fltr[, Value]

        output$kit <- renderText({names[folders == "kit"]})
        output$alprob <- renderText({gtc_prob_value[gtc_prob_param == "Method"]})
        output$theta <- renderText({gtc_prob_value[gtc_prob_param == "Theta"]})
        output$bsr_fltr <- renderText({paste0("Filter of back stutter ratio: ", gtc_fltr_value[gtc_fltr_type == "BSR"])})
        output$fsr_fltr <- renderText({paste0("Filter of forward stutter ratio: ", gtc_fltr_value[gtc_fltr_type == "FSR"])})
        output$dsr_fltr <- renderText({paste0("Filter of double-back stutter ratio: ", gtc_fltr_value[gtc_fltr_type == "DSR"])})
        output$m2sr_fltr <- renderText({paste0("Filter of minus 2-nt stutter ratio: ", gtc_fltr_value[gtc_fltr_type == "M2SR"])})
      })

      observeEvent(input$act_ae, {
        dt_param <- rv_lr_setting$dt_param
        setkey(dt_param, Parameter)
        rv_lr_setting$dt_ae <- dt_param["AE"]

        showModal(modalDialog(
          title = "Locus-specific amplification efficiency",
          dataTableOutput(session$ns("dt_ae")),
          easyClose = TRUE,
          footer = NULL
        ))
      })

      output$dt_ae <- renderDataTable(server = FALSE, {
        dt_ae <- rv_lr_setting$dt_ae
        datatable(
          dt_ae,
          selection = "none",
          options = list(dom = "t", iDisplayLength = nrow(dt_ae), ordering = FALSE, autoWidth = TRUE,
                         columnDefs = list(list(targets = c(1, 2, 4:16), visible = FALSE))),
          rownames = FALSE
        ) %>%
          formatRound(columns = c("Alpha"), digits = 3)
      })

      observeEvent(input$act_bsr, {
        print(input$act_bsr)

        # Extract parameters for BSR
        dt_param <- rv_lr_setting$dt_param
        setkey(dt_param, Parameter)
        rv_lr_setting$dt_bsr <- dt_bsr <- dt_param["BSR"]
        markers <- dt_bsr[, Marker]

        showModal(modalDialog(
          title = "Back stutter ratio",
          selectInput(session$ns("mk_bsr"), label = h4("Marker"), choices = markers, selected = markers[1], selectize = FALSE),
          h4("Parameter"),
          dataTableOutput(session$ns("dt_param_bsr")),
          br(),
          h4("Corrected alleles and average back stutter ratios"),
          dataTableOutput(session$ns("dt_alcor_bsr")),
          easyClose = TRUE,
          footer = NULL
        ))
      })

      observeEvent(input$mk_bsr, {
        # Extract parameters
        dt_bsr <- rv_lr_setting$dt_bsr
        mk_bsr <- input$mk_bsr
        setkey(dt_bsr, Marker)
        dt_bsr_select <- dt_bsr[mk_bsr]
        model_bsr_select <- dt_bsr_select[, Model]

        if(model_bsr_select != "Not consider"){
          # Extract corrected alleles
          dt_alcor <- rv_lr_setting$dt_alcor
          setkey(dt_alcor, Marker)
          dt_alcor_select <- dt_alcor[mk_bsr]
        }else{
          dt_alcor_select <- NULL
        }

        output$dt_param_bsr <- renderDataTable(server = FALSE, {
          datatable(
            dt_bsr_select,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_bsr_select), ordering = FALSE, autoWidth = TRUE,
                           columnDefs = list(list(targets = c(0, 1, 3, 9:16), visible = FALSE))),
            rownames = FALSE
          ) %>%
            formatRound(columns = c("Myu", "Slope", "Intercept"), digits = 3) %>%
            formatRound(columns = c("Sigma", "MinUS"), digits = 0)
        })

        output$dt_alcor_bsr <- renderDataTable(server = FALSE, {
          col_names <- names(dt_alcor_select)
          col_names[4] <- "Corrected allele"
          col_names[5] <- "Average back stutter ratio"
          datatable(
            dt_alcor_select,
            colnames = col_names,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_alcor_select), ordering = FALSE, autoWidth = TRUE,
                           columnDefs = list(list(targets = c(0, 1, 5:10), visible = FALSE))),
            rownames = FALSE
          ) %>%
            formatRound(columns = c("Corrected_BSR", "Average_BSR"), digits = 3)
        })
      })

      observeEvent(input$act_fsr, {
        # Extract parameters for FSR
        dt_param <- rv_lr_setting$dt_param
        setkey(dt_param, Parameter)
        rv_lr_setting$dt_fsr <- dt_fsr <- dt_param["FSR"]
        markers <- dt_fsr[, Marker]

        showModal(modalDialog(
          title = "Forward stutter ratio",
          selectInput(session$ns("mk_fsr"), label = h4("Marker"), choices = markers, selected = markers[1], selectize = FALSE),
          h4("Parameter"),
          dataTableOutput(session$ns("dt_param_fsr")),
          br(),
          h4("Corrected alleles and average forward stutter ratios"),
          dataTableOutput(session$ns("dt_alcor_fsr")),
          easyClose = TRUE,
          footer = NULL
        ))
      })

      observeEvent(input$mk_fsr, {
        # Extract parameters
        dt_fsr <- rv_lr_setting$dt_fsr
        mk_fsr <- input$mk_fsr
        setkey(dt_fsr, Marker)
        dt_fsr_select <- dt_fsr[mk_fsr]
        model_fsr_select <- dt_fsr_select[, Model]

        if(model_fsr_select != "Not consider"){
          # Extract corrected alleles
          dt_alcor <- rv_lr_setting$dt_alcor
          setkey(dt_alcor, Marker)
          dt_alcor_select <- dt_alcor[mk_fsr]
        }else{
          dt_alcor_select <- NULL
        }

        output$dt_param_fsr <- renderDataTable(server = FALSE, {
          datatable(
            dt_fsr_select,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_fsr_select), ordering = FALSE, autoWidth = TRUE,
                           columnDefs = list(list(targets = c(0, 1, 3, 9:16), visible = FALSE))),
            rownames = FALSE
          ) %>%
            formatRound(columns = c("Myu", "Slope", "Intercept"), digits = 3) %>%
            formatRound(columns = c("Sigma", "MinUS"), digits = 0)
        })

        output$dt_alcor_fsr <- renderDataTable(server = FALSE, {
          col_names <- names(dt_alcor_select)
          col_names[7] <- "Corrected allele"
          col_names[8] <- "Average forward stutter ratio"
          datatable(
            dt_alcor_select,
            colnames = col_names,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_alcor_select), ordering = FALSE, autoWidth = TRUE,
                           columnDefs = list(list(targets = c(0, 1, 3:5, 8:10), visible = FALSE))),
            rownames = FALSE
          ) %>%
            formatRound(columns = c("Corrected_FSR", "Average_FSR"), digits = 3)
        })
      })

      observeEvent(input$act_dsr, {
        # Extract parameters for DSR
        dt_param <- rv_lr_setting$dt_param
        setkey(dt_param, Parameter)
        rv_lr_setting$dt_dsr <- dt_dsr <- dt_param["DSR"]
        markers <- dt_dsr[, Marker]

        showModal(modalDialog(
          title = "Double-back stutter ratio",
          selectInput(session$ns("mk_dsr"), label = h4("Marker"), choices = markers, selected = markers[1], selectize = FALSE),
          h4("Parameter"),
          dataTableOutput(session$ns("dt_param_dsr")),
          br(),
          h4("Corrected alleles and average double-back stutter ratios"),
          dataTableOutput(session$ns("dt_alcor_dsr")),
          easyClose = TRUE,
          footer = NULL
        ))
      })

      observeEvent(input$mk_dsr, {
        # Extract parameters
        dt_dsr <- rv_lr_setting$dt_dsr
        mk_dsr <- input$mk_dsr
        setkey(dt_dsr, Marker)
        dt_dsr_select <- dt_dsr[mk_dsr]
        model_dsr_select <- dt_dsr_select[, Model]

        if(model_dsr_select != "Not consider"){
          # Extract corrected alleles
          dt_alcor <- rv_lr_setting$dt_alcor
          setkey(dt_alcor, Marker)
          dt_alcor_select <- dt_alcor[mk_dsr]
        }else{
          dt_alcor_select <- NULL
        }

        output$dt_param_dsr <- renderDataTable(server = FALSE, {
          datatable(
            dt_dsr_select,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_dsr_select), ordering = FALSE, autoWidth = TRUE,
                           columnDefs = list(list(targets = c(0, 1, 3, 9:16), visible = FALSE))),
            rownames = FALSE
          ) %>%
            formatRound(columns = c("Myu", "Slope", "Intercept"), digits = 3) %>%
            formatRound(columns = c("Sigma", "MinUS"), digits = 0)
        })

        output$dt_alcor_dsr <- renderDataTable(server = FALSE, {
          col_names <- names(dt_alcor_select)
          col_names[10] <- "Corrected allele"
          col_names[11] <- "Average double-back stutter ratio"
          datatable(
            dt_alcor_select,
            colnames = col_names,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_alcor_select), ordering = FALSE, autoWidth = TRUE,
                           columnDefs = list(list(targets = c(0, 1, 3:8), visible = FALSE))),
            rownames = FALSE
          ) %>%
            formatRound(columns = c("Corrected_DSR", "Average_DSR"), digits = 3)
        })
      })

      observeEvent(input$act_m2sr, {
        # Extract parameters for M2SR
        dt_param <- rv_lr_setting$dt_param
        setkey(dt_param, Parameter)
        rv_lr_setting$dt_m2sr <- dt_m2sr <- dt_param["M2SR"]
        markers <- dt_m2sr[, Marker]

        showModal(modalDialog(
          title = "Minus 2-nt stutter ratio",
          selectInput(session$ns("mk_m2sr"), label = h4("Marker"), choices = markers, selected = markers[1], selectize = FALSE),
          dataTableOutput(session$ns("dt_param_m2sr")),
          easyClose = TRUE,
          footer = NULL
        ))
      })

      observeEvent(input$mk_m2sr, {
        # Extract parameters
        dt_m2sr <- rv_lr_setting$dt_m2sr
        mk_m2sr <- input$mk_m2sr
        setkey(dt_m2sr, Marker)
        dt_m2sr_select <- dt_m2sr[mk_m2sr]
        model_m2sr_select <- dt_m2sr_select[, Model]

        output$dt_param_m2sr <- renderDataTable(server = FALSE, {
          datatable(
            dt_m2sr_select,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_m2sr_select), ordering = FALSE, autoWidth = TRUE,
                           columnDefs = list(list(targets = c(0, 1, 3, 5, 6, 8:16), visible = FALSE))),
            rownames = FALSE
          ) %>%
            formatRound(columns = c("Myu"), digits = 3) %>%
            formatRound(columns = c("Sigma"), digits = 0)
        })
      })

      observeEvent(input$act_gamma, {
        dt_param <- rv_lr_setting$dt_param
        setkey(dt_param, Parameter)
        rv_lr_setting$dt_gamma <- dt_param["Gamma"]

        showModal(modalDialog(
          title = "Gamma distribution",
          dataTableOutput(session$ns("dt_gamma")),
          easyClose = TRUE,
          footer = NULL
        ))
      })

      output$dt_gamma <- renderDataTable(server = FALSE, {
        dt_gamma <- rv_lr_setting$dt_gamma
        datatable(
          dt_gamma,
          selection = "none",
          options = list(dom = "t", iDisplayLength = nrow(dt_gamma), ordering = FALSE, autoWidth = TRUE,
                         columnDefs = list(list(targets = c(0:8, 10:16), visible = FALSE))),
          rownames = FALSE
        ) %>%
          formatRound(columns = c("Cv"), digits = 3)
      })

      observeEvent(input$act_analyze_lr, {
        error_message <- ""
        hp <- input$hp
        hd <- input$hd
        refs <- rv_lr_setting$refs
        dt_ref <- rv_lr_setting$dt_ref
        use_deconvo <- input$use_deconvo

        if(use_deconvo){
          data_deconvo_hp <- data_deconvo_hd <- rv_lr_setting$data_deconvo
          data_hp <- assign_log10like(data_deconvo_hp, dt_ref, hp, refs)
          data_hd <- assign_log10like(data_deconvo_hd, dt_ref, hd, refs)
          dt_lr <- make_dt_lr(data_deconvo_hd$markers_auto,
                              data_hp$log10like_adopt, data_hp$combined_log10like_adopt,
                              data_hd$log10like_adopt, data_hd$combined_log10like_adopt)
        }else{
          # Fix setting
          hnc_lr_hp <- input$noc_lr_hp
          hnc_lr_hd <- input$noc_lr_hd
          dt_csp <- rv_lr_setting$dt_csp
          dt_ref <- rv_lr_setting$dt_ref
          dt_af <- rv_lr_setting$dt_af
          anameth <- rv_lr_setting$anameth
          dt_kit <- rv_lr_setting$dt_kit
          dt_at <- rv_lr_setting$dt_at
          dt_param <- rv_lr_setting$dt_param
          dt_alcor <- rv_lr_setting$dt_alcor
          dt_gtc_fltr <- rv_lr_setting$dt_gtc_fltr
          dt_gtc_prob <- rv_lr_setting$dt_gtc_prob
          warmup <- rv_lr_setting$warmup
          sampling <- rv_lr_setting$sampling
          chains <- rv_lr_setting$chains

          # Check input files and setting
          error_message <- check_main_file(dt_csp, dt_ref, dt_af, dt_kit, dt_at, dt_alcor, method = "deconvo")
          if(error_message != ""){
            showModal(modalDialog(title = "Error", HTML(error_message), easyClose = TRUE, footer = NULL))
          }else{
            data_deconvo_hp <- perform_deconvo(path_pack, dt_csp, dt_af, dt_kit, dt_at, dt_param, dt_alcor, dt_gtc_fltr, dt_gtc_prob, hnc_lr_hp,
                                               warmup, sampling, chains)
            if(hnc_lr_hp == hnc_lr_hd){
              data_deconvo_hd <- data_deconvo_hp
            }else{
              data_deconvo_hd <- perform_deconvo(path_pack, dt_csp, dt_af, dt_kit, dt_at, dt_param, dt_alcor, dt_gtc_fltr, dt_gtc_prob, hnc_lr_hd,
                                                 warmup, sampling, chains)
            }

            data_hp <- assign_log10like(data_deconvo_hp, dt_ref, hp, refs)
            data_hd <- assign_log10like(data_deconvo_hd, dt_ref, hd, refs)
            dt_lr <- make_dt_lr(data_deconvo_hd$markers_auto,
                                data_hp$log10like_adopt, data_hp$combined_log10like_adopt,
                                data_hd$log10like_adopt, data_hd$combined_log10like_adopt)
          }
        }

        if(error_message == ""){
          data_lr <- list("dt_mpdeg_hp" = data_hp$dt_mpdeg_hyp,
                          "dt_mpdeg_hd" = data_hd$dt_mpdeg_hyp,
                          "dt_pg_list_hp" = data_hp$dt_pg_list_hyp,
                          "dt_pg_list_hd" = data_hd$dt_pg_list_hyp,
                          "dt_lr" = dt_lr)
          rv_lr$data_lr <- data_lr

          updateTabsetPanel(session = session_lr_top, inputId = "tabset_lr", selected = "Result")
        }
      })

      return(rv_lr)
    }
  )
}
