#' hd_true_setting_ui
#'
#' @description The function to create the ui module for setting conditions of the Hd true test
hd_true_setting_ui <- function(id){
  ns <- NS(id)

  tabPanel("Setting",
           useWaiter(),
           br(),
           sidebarLayout(
             sidebarPanel(width = 3,
                          uiOutput(ns("use_deconvo")),
                          br(),
                          numericInput(ns("noc_hd_true"), label = "Number of contributors", value = 1),
                          br(),
                          uiOutput(ns("anameth_hd_true")),
                          br(),
                          actionButton(ns("act_advance_hd_true"), label = "Advanced setting"),
                          br(),
                          br(),
                          numericInput(ns("n_simu"), label = h4("Number of simulations"), value = 100),
                          br(),
                          actionButton(ns("act_analyze_hd_true"), label = "Analyze", class = "btn btn-primary btn-lg")
             ),
             mainPanel(width = 9,
                       tags$head( #CSS added to shrink the column with
                         tags$style('td.small .shiny-input-container{width:auto;}
                                     td.small{width:30px;}
                                     th.small{font-size: 14px !important;}'
                         )
                       ),
                       dataTableOutput(ns("dt_info_mk_display_hd_true")))
           )
  )
}

#' hd_true_setting_server
#'
#' @description The function to create the server module for setting conditions of the Hd true test
hd_true_setting_server <- function(id, path_pack, rv_file, rv_deconvo){
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

      rv_hd_true_setting <- reactiveValues()
      rv_hd_true_setting$data_deconvo <- NULL
      rv_hd_true_setting$dt_csp <- NULL
      rv_hd_true_setting$dt_af <- NULL
      rv_hd_true_setting$fn_csp <- NULL
      rv_hd_true_setting$fn_af <- NULL
      rv_hd_true_setting$dt_anameth <- NULL
      rv_hd_true_setting$dt_info_mk <- NULL
      rv_hd_true_setting$warmup <- 1200
      rv_hd_true_setting$sampling <- 300
      rv_hd_true_setting$chains <- 4

      rv_hd_true <- reactiveValues()
      rv_hd_true$lr_rand <- NULL
      rv_hd_true$data_deconvo <- NULL

      observe({
        req(rv_file)
        rv_hd_true_setting$dt_csp <- rv_file$dt_csp
        rv_hd_true_setting$dt_af <- rv_file$dt_af
        rv_hd_true_setting$fn_csp <- rv_file$fn_csp
        rv_hd_true_setting$fn_af <- rv_file$fn_af
      })

      observe({
        req(rv_deconvo)
        rv_hd_true_setting$data_deconvo <- rv_deconvo$data_deconvo
      })

      output$use_deconvo <- renderUI({
        if(is.null(rv_hd_true_setting$data_deconvo)){
          disabled(checkboxInput(session$ns("use_deconvo"), label = "Use the result of deconvolution", value = FALSE))
        }else{
          checkboxInput(session$ns("use_deconvo"), label = "Use the result of deconvolution", value = TRUE)
        }
      })

      observeEvent(input$use_deconvo, {
        use_deconvo <- input$use_deconvo
        data_deconvo <- rv_hd_true_setting$data_deconvo
        if(use_deconvo){
          disable("noc_hd_true")
          disable("anameth_hd_true")
          disable("act_advance_hd_true")
          updateNumericInput(session, "noc_hd_true", value = data_deconvo$hnc)
          updateSelectInput(session, "anameth_hd_true", selected = data_deconvo$anameth)
        }else{
          enable("noc_hd_true")
          enable("anameth_hd_true")
          enable("act_advance_hd_true")
        }
      })

      observeEvent(input$noc_hd_true, {
        noc_hd_true <- input$noc_hd_true
        if(!is.integer(noc_hd_true) || noc_hd_true < 1 || noc_hd_true > 4){
          showFeedbackDanger(inputId = "noc_hd_true", text = "Only 1, 2, 3, and 4 is allowed.")
          disable("act_analyze_hd_true")
        }else{
          hideFeedback("noc_hd_true")
          enable("act_analyze_hd_true")
        }
      })

      observeEvent(input$act_advance_hd_true, {
        showModal(modalDialog(
          title = "Advanced setting",
          fluidRow(
            column(6,
                   h4("MCMC conditions"),
                   numericInput(session$ns("warmup"), label = "Warmup iterations", value = rv_hd_true_setting$warmup),
                   numericInput(session$ns("sampling"), label = "Sampling iterations", value = rv_hd_true_setting$sampling),
                   numericInput(session$ns("chains"), label = "Number of chains", value = rv_hd_true_setting$chains)
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
          rv_hd_true_setting$warmup <- warmup
          rv_hd_true_setting$sampling <- sampling
          rv_hd_true_setting$chains <- chains
          removeModal()
        }
      })

      output$anameth_hd_true <- renderUI({
        anameth_all <- list.files(paste0(path_pack, "/extdata/anameth"))
        anameth_all <- gsub(".{4}$", "", anameth_all)
        selectInput(session$ns("anameth_hd_true"), label = h4("Analysis method"), choices = anameth_all, selected = anameth_all[1], selectize = FALSE)
      })

      observeEvent(input$anameth_hd_true, {
        # Load the selected analysis method
        rv_hd_true_setting$dt_anameth <- dt_anameth <- fread(paste0(path_pack, "/extdata/anameth/", input$anameth_hd_true, ".csv"))

        # Create a data.table for information on markers
        dt_at <- data.table(Marker = dt_anameth$Subcategory[dt_anameth$Category == "AT"],
                            AT = as.numeric(dt_anameth$Value[dt_anameth$Category == "AT"]))
        dt_info_mk <- copy(dt_at)
        dt_info_mk <- cbind(
          Ignore = rep(FALSE, nrow(dt_info_mk)),
          dt_info_mk
        )
        rv_hd_true_setting$dt_info_mk <- dt_info_mk

        # Create a displayed data.table for information on markers
        dt_info_mk_display_hd_true <- copy(dt_at)
        dt_info_mk_display_hd_true <- cbind(
          Ignore = shinyInput(checkboxInput, nrow(dt_info_mk_display_hd_true), "checkb"),
          dt_info_mk_display_hd_true
        )

        output$dt_info_mk_display_hd_true <- renderDataTable(server = FALSE, {
          datatable(
            dt_info_mk_display_hd_true,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_info_mk_display_hd_true), ordering = FALSE, autoWidth = TRUE, scrollX = TRUE, scrollY = "calc(100vh - 200px)",
                           columnDefs = list(list(targets = 0, className = "small"), list(targets = 1:2, width = "100px"))),
            rownames = FALSE,
            escape = FALSE,
            callback = JS(js("dt_info_mk_display_hd_true", session$ns))
          )
        })
      })

      observeEvent(input[["dt_info_mk_display_hd_true_cell_edit"]], {
        info <- input[["dt_info_mk_display_hd_true_cell_edit"]]
        dt_info_mk <- rv_hd_true_setting$dt_info_mk
        dt_info_mk[[info[1, "col"]]][info[1, "row"]] <- info[1, "value"]
        rv_hd_true_setting$dt_info_mk <- dt_info_mk
      })

      observeEvent(input$act_analyze_hd_true, {
        waiter_show(html = spin_3k(), color = "white")

        # Fix setting
        dt_csp <- rv_hd_true_setting$dt_csp
        dt_af <- rv_hd_true_setting$dt_af
        dt_anameth <- rv_hd_true_setting$dt_anameth
        use_deconvo <- input$use_deconvo
        n_simu <- input$n_simu

        # Extract information on the kit
        dt_kit <- fread(paste0(path_pack, "/extdata/kit/", dt_anameth$Value[dt_anameth$Category == "Kit"], ".csv"))
        dt_kit <- change_data_type_kit(dt_kit)
        dt_kit_auto <- extract_kit_auto(dt_kit)
        markers_auto <- unique(dt_kit_auto[, Marker])

        error_message <- check_main_file(dt_csp, NULL, dt_af, markers_auto, method = "deconvo")
        if(error_message != ""){
          waiter_hide()
          showModal(modalDialog(title = "Error", HTML(error_message), easyClose = TRUE, footer = NULL))
        }else{
          if(use_deconvo){
            data_deconvo <- rv_hd_true_setting$data_deconvo
          }else{
            # Fix setting
            hnc <- input$noc_hd_true
            warmup <- rv_hd_true_setting$warmup
            sampling <- rv_hd_true_setting$sampling
            chains <- rv_hd_true_setting$chains
            anameth <- input$anameth_hd_true

            load(paste0(path_pack, "/extdata/param/", dt_anameth$Value[dt_anameth$Category == "Parameter"], ".Rdata"))
            dt_al_prob <- data.table(Subcategory = dt_anameth$Subcategory[dt_anameth$Category == "Allele_probability"],
                                     Value = dt_anameth$Value[dt_anameth$Category == "Allele_probability"])
            dt_gt_fltr <- data.table(Subcategory = dt_anameth$Subcategory[dt_anameth$Category == "Genotype_filter"],
                                     Value = dt_anameth$Value[dt_anameth$Category == "Genotype_filter"])
            dt_dropin <- data.table(Subcategory = dt_anameth$Subcategory[dt_anameth$Category == "Drop-in"],
                                    Value = dt_anameth$Value[dt_anameth$Category == "Drop-in"])

            dt_info_mk <- rv_hd_true_setting$dt_info_mk

            # Check registration of the stutter ratio
            message_sr <- check_csp_sr(dt_csp, dt_info_mk, dt_sr)
            if(message_sr == "error"){
              data_deconvo <- NULL
              waiter_hide()
              showModal(modalDialog(title = "Error", "All markers are ignored!", easyClose = TRUE, footer = NULL))
            }else{
              data_deconvo <- perform_deconvo(dt_csp, dt_af, dt_kit, dt_info_mk, dt_par, dt_sr, dt_al_prob, dt_gt_fltr, dt_dropin,
                                              hnc, warmup, sampling, chains)
            }
          }
          if(!is.null(data_deconvo)){
            if(!is.null(data_deconvo$error_message)){
              waiter_hide()
              showModal(modalDialog(title = "Error", HTML(data_deconvo$error_message), easyClose = TRUE, footer = NULL))
            }else{
              # Extract data from data_deconvo
              hnc <- data_deconvo$hnc
              markers_use <- data_deconvo$markers_use
              peaks_list <- data_deconvo$peaks_list
              pop_al_list <- data_deconvo$pop_al_list
              pop_prob_list <- data_deconvo$pop_prob_list
              theta <- data_deconvo$theta
              gtc_use <- data_deconvo$gtc_use
              lpd_gtc_use_1st <- data_deconvo$lpd_gtc_use_1st
              lpd_gtc_use_2nd <- data_deconvo$lpd_gtc_use_2nd

              method_al_prob <- dt_anameth$Value[dt_anameth$Category == "Allele_probability" & dt_anameth$Subcategory == "Method"]
              tmp <- set_al_prob_randgt(markers_use, dt_af, method_al_prob)
              pop_al_list_randgt <- tmp$pop_al_list
              pop_prob_list_randgt <- tmp$pop_prob_list
              randgt <- t(generate_randgt(n_simu, pop_al_list_randgt, pop_prob_list_randgt))

              lr_rand <- rep(NA, n_simu)
              for(i in 1:n_simu){
                ref1 <- randgt[, c(2 * i - 1, 2 * i)]
                ref1 <- t(apply(ref1, 1, sort))
                refs <- paste0("Ref", i)
                dt_ref_rand <- data.table(SampleName = rep(refs, length(markers_use)),
                                          Marker = markers_use,
                                          Allele1 = ref1[, 1],
                                          Allele2 = ref1[, 2])
                data_hp <- assign_log10like_randgt(hnc, markers_use, peaks_list, pop_al_list, pop_prob_list, theta, gtc_use, lpd_gtc_use_1st, lpd_gtc_use_2nd, dt_ref_rand, refs, refs)
                data_hd <- assign_log10like_randgt(hnc, markers_use, peaks_list, pop_al_list, pop_prob_list, theta, gtc_use, lpd_gtc_use_1st, lpd_gtc_use_2nd, dt_ref_rand, NULL, refs)
                lr_rand[i] <- 10^(data_hp$combined_log10like_adopt - data_hd$combined_log10like_adopt)
              }
              rv_hd_true$lr_rand <- lr_rand
              rv_hd_true$data_deconvo <- data_deconvo

              waiter_hide()
            }
          }
        }
      })

      return(rv_hd_true)
    }
  )
}
