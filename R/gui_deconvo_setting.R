#' deconvo_setting_ui
#'
#' @description The function to create the ui module for setting conditions of deconvolution
deconvo_setting_ui <- function(id){
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
                      h4("Number of contributors"),
                      numericInput(ns("noc_from"), label = "From", value = 1),
                      numericInput(ns("noc_to"), label = "To", value = 1),
                      br(),
                      h4("Analysis method"),
                      shinyFilesButton(ns("load_anameth_deconvo"), label = "Load", title = "Select a file", multiple = FALSE),
                      uiOutput(ns("anameth_deconvo")),
                      br(),
                      br(),
                      actionButton(ns("act_analyze_deconvo"), label = "Analyze", class = "btn btn-primary btn-lg")
                    )
             ),
             column(3,
                    h4(textOutput(ns("title_ignore_mk"))),
                    rHandsontableOutput(ns("dt_ignore_mk"))
             ),
             column(3,
                    h4(textOutput(ns("title_at"))),
                    dataTableOutput(ns("anameth_at")),
                    br()
             ),
             column(3,
                    h4(textOutput(ns("title_al_gt_prob"))),
                    dataTableOutput(ns("anameth_algt_prob")),
                    br(),
                    h4(textOutput(ns("title_gt_fltr"))),
                    dataTableOutput(ns("anameth_gt_fltr")),
                    br(),
                    h4(textOutput(ns("title_dropin"))),
                    dataTableOutput(ns("anameth_dropin")),
                    br(),
                    h4(textOutput(ns("title_mcmc"))),
                    dataTableOutput(ns("anameth_mcmc")),
                    br()
             )
           )
  )
}

#' deconvo_setting_server
#'
#' @description The function to create the server module for setting conditions of deconvolution
deconvo_setting_server <- function(id, path_pack, rv_file){
  moduleServer(
    id,
    function(input, output, session){
      rv_deconvo_setting <- reactiveValues()
      rv_deconvo_setting$dt_csp <- NULL
      rv_deconvo_setting$dt_af <- NULL
      rv_deconvo_setting$fn_csp <- NULL
      rv_deconvo_setting$fn_af <- NULL
      rv_deconvo_setting$fn_anameth_deconvo <- ""
      rv_deconvo_setting$dt_anameth <- NULL
      rv_deconvo_setting$bool_correct_noc <- TRUE

      rv_deconvo <- reactiveValues()
      rv_deconvo$data_deconvo <- NULL

      observe({
        req(rv_file)
        rv_deconvo_setting$dt_csp <- rv_file$dt_csp
        rv_deconvo_setting$dt_af <- rv_file$dt_af
        rv_deconvo_setting$fn_csp <- rv_file$fn_csp
        rv_deconvo_setting$fn_af <- rv_file$fn_af
      })

      observeEvent(input$noc_from | input$noc_to, {
        noc_from <- input$noc_from
        noc_to <- input$noc_to

        hideFeedback("noc_from")
        hideFeedback("noc_to")
        bool_correct_noc <- TRUE

        if(!is.integer(noc_from) || noc_from < 1 || noc_from > 4){
          showFeedbackDanger(inputId = "noc_from", text = "Only 1, 2, 3, and 4 is allowed.")
          bool_correct_noc <- FALSE
        }

        if(!is.integer(noc_to) || noc_to < 1 || noc_to > 4){
          showFeedbackDanger(inputId = "noc_to", text = "Only 1, 2, 3, and 4 is allowed.")
          bool_correct_noc <- FALSE
        }

        if(is.integer(noc_from) && is.integer(noc_to) && noc_from > noc_to){
          showFeedbackDanger(inputId = "noc_from", text = "The From number must not exceed the To number.")
          showFeedbackDanger(inputId = "noc_to", text = "The From number must not exceed the To number.")
          bool_correct_noc <- FALSE
        }

        rv_deconvo_setting$bool_correct_noc <- bool_correct_noc

        if(bool_correct_noc){
          if(rv_deconvo_setting$fn_anameth_deconvo == ""){
            disable("act_analyze_deconvo")
          }else{
            enable("act_analyze_deconvo")
          }
        }else{
          disable("act_analyze_deconvo")
        }
      })

      shinyFileChoose(input, "load_anameth_deconvo", roots = c(documents = paste0(path_pack, "/extdata/anameth")), session = session, filetypes = c("csv"))

      observeEvent(input$load_anameth_deconvo, {
        req(input$load_anameth_deconvo)
        path_anameth <- parseFilePaths(roots = c(documents = paste0(path_pack, "/extdata/anameth")), input$load_anameth_deconvo)$datapath
        if(length(path_anameth) > 0){
          rv_deconvo_setting$fn_anameth_deconvo <- basename(path_anameth)
          if(rv_deconvo_setting$bool_correct_noc){
            enable("act_analyze_deconvo")
          }else{
            disable("act_analyze_deconvo")
          }
        }else{
          disable("act_analyze_deconvo")
        }
      })

      output$anameth_deconvo <- renderUI({
        disabled(textInput(session$ns("anameth_deconvo"), label = "", value = rv_deconvo_setting$fn_anameth_deconvo))
      })

      observeEvent(input$anameth_deconvo, {
        fn_anameth_deconvo <- input$anameth_deconvo
        if(fn_anameth_deconvo != ""){
          rv_deconvo_setting$dt_anameth <- dt_anameth <- fread(paste0(path_pack, "/extdata/anameth/", fn_anameth_deconvo))
          markers_auto <- dt_anameth$Subcategory[dt_anameth$Category == "AT"]

          dt_ignore_mk <- data.table(Ignore = rep(FALSE, length(markers_auto)),
                                     Marker = markers_auto)

          output$dt_ignore_mk <- renderRHandsontable({
            rhandsontable(dt_ignore_mk, rowHeaders = FALSE) %>%
              hot_table(fillHandle = FALSE, rowHeights = 30) %>%
              hot_col("Ignore", type = "checkbox", width = 100) %>%
              hot_col("Marker", type = "text", readOnly = TRUE, width = 100)
          })

          output$title_ignore_mk <- renderText({"Ignored marker"})
          output$title_at <- renderText({"Analytical threshold"})
          output$title_al_gt_prob <- renderText({"Genotype probability"})
          output$title_gt_fltr <- renderText({"Genotype filter"})
          output$title_dropin <- renderText({"Drop-in"})
          output$title_mcmc <- renderText({"MCMC"})

          output$anameth_at <- renderRHandsontable({
            dt_at <- data.table(Marker = dt_anameth$Subcategory[dt_anameth$Category == "AT"],
                                AT = dt_anameth$Value[dt_anameth$Category == "AT"])

            datatable(
              dt_at,
              selection = "none",
              options = list(dom = "t", iDisplayLength = nrow(dt_at), ordering = FALSE, autowidth = TRUE),
              rownames = FALSE,
              escape = FALSE
            ) %>%
              formatStyle(names(dt_at), lineHeight = "80%")
          })

          output$anameth_algt_prob <- renderRHandsontable({
            dt_algt_prob <- data.table(Subcategory = dt_anameth$Subcategory[dt_anameth$Category == "Genotype probability"],
                                       Value = dt_anameth$Value[dt_anameth$Category == "Genotype probability"])

            datatable(
              dt_algt_prob,
              selection = "none",
              options = list(dom = "t", iDisplayLength = nrow(dt_algt_prob), ordering = FALSE, autowidth = TRUE),
              rownames = FALSE,
              escape = FALSE
            ) %>%
              formatStyle(names(dt_algt_prob), lineHeight = "80%")
          })

          output$anameth_gt_fltr <- renderRHandsontable({
            dt_gt_fltr <- data.table(Subcategory = dt_anameth$Subcategory[dt_anameth$Category == "Genotype filter"],
                                     Value = dt_anameth$Value[dt_anameth$Category == "Genotype filter"])

            datatable(
              dt_gt_fltr,
              selection = "none",
              options = list(dom = "t", iDisplayLength = nrow(dt_gt_fltr), ordering = FALSE, autowidth = TRUE),
              rownames = FALSE,
              escape = FALSE
            ) %>%
              formatStyle(names(dt_gt_fltr), lineHeight = "80%")
          })

          output$anameth_dropin <- renderRHandsontable({
            dt_dropin <- data.table(Subcategory = dt_anameth$Subcategory[dt_anameth$Category == "Drop-in"],
                                    Value = dt_anameth$Value[dt_anameth$Category == "Drop-in"])

            datatable(
              dt_dropin,
              selection = "none",
              options = list(dom = "t", iDisplayLength = nrow(dt_dropin), ordering = FALSE, autowidth = TRUE),
              rownames = FALSE,
              escape = FALSE
            ) %>%
              formatStyle(names(dt_dropin), lineHeight = "80%")
          })

          output$anameth_mcmc <- renderRHandsontable({
            dt_mcmc <- data.table(Subcategory = dt_anameth$Subcategory[dt_anameth$Category == "MCMC"],
                                  Value = dt_anameth$Value[dt_anameth$Category == "MCMC"])

            datatable(
              dt_mcmc,
              selection = "none",
              options = list(dom = "t", iDisplayLength = nrow(dt_mcmc), ordering = FALSE, autowidth = TRUE),
              rownames = FALSE,
              escape = FALSE
            ) %>%
              formatStyle(names(dt_mcmc), lineHeight = "80%")
          })
        }
      })

      observeEvent(input$act_analyze_deconvo, {
        waiter_show(html = spin_3k(), color = "white")

        # Fix setting
        hnc_from <- input$noc_from
        hnc_to <- input$noc_to
        fn_csp <- rv_deconvo_setting$fn_csp
        dt_csp <- rv_deconvo_setting$dt_csp
        fn_af <- rv_deconvo_setting$fn_af
        dt_af <- rv_deconvo_setting$dt_af
        anameth <- rv_deconvo_setting$fn_anameth_deconvo
        dt_anameth <- rv_deconvo_setting$dt_anameth

        # Extract information on the kit
        dt_kit <- fread(paste0(path_pack, "/extdata/kit/", dt_anameth$Name[dt_anameth$Category == "Kit"], ".csv"))
        dt_kit <- change_data_type_kit(dt_kit)
        dt_kit_auto <- extract_kit_auto(dt_kit)
        markers_auto <- unique(dt_kit_auto[, Marker])

        load(paste0(path_pack, "/extdata/param/", dt_anameth$Name[dt_anameth$Category == "Parameter"], ".RData"))
        load(paste0(path_pack, "/extdata/stutter_file/", dt_anameth$Name[dt_anameth$Category == "Stutter file"], ".RData"))

        dt_ignore_mk <- hot_to_r(input$dt_ignore_mk)

        # Check input files
        error_message <- check_main_file(dt_csp, NULL, dt_af, markers_auto, method = "deconvo")
        if(error_message != ""){
          waiter_hide()
          showModal(modalDialog(title = "Error", HTML(error_message), easyClose = TRUE, footer = NULL))
        }else{
          # Check registration of the stutter ratio
          message_sr <- check_csp_sr(dt_csp, dt_ignore_mk, dt_anameth, list_dt_average_sr)

          if(message_sr == "error"){
            waiter_hide()
            showModal(modalDialog(title = "Error", "All markers are ignored!", easyClose = TRUE, footer = NULL))
          }else{
            hnc_all <- hnc_from:hnc_to
            data_deconvo <- perform_deconvo_all_noc(dt_csp, dt_af, dt_kit_auto, dt_ignore_mk, dt_anameth, median_alpha, median_cv, list_dt_average_sr, hnc_all)
            data_deconvo$hnc_all <- hnc_all
            data_deconvo$dt_csp <- dt_csp
            data_deconvo$dt_af <- dt_af
            data_deconvo$fn_csp <- fn_csp
            data_deconvo$fn_af <- fn_af
            data_deconvo$anameth <- anameth
            data_deconvo$dt_anameth <- dt_anameth
            data_deconvo$dt_kit <- dt_kit
            data_deconvo$dt_ignore_mk <- dt_ignore_mk
            data_deconvo$alpha_auto <- median_alpha
            data_deconvo$cv <- median_cv
            data_deconvo$list_dt_average_sr <- list_dt_average_sr
            data_deconvo$message_sr <- message_sr
            rv_deconvo$data_deconvo <- data_deconvo

            waiter_hide()
          }
        }
      })

      return(rv_deconvo)
    }
  )
}
