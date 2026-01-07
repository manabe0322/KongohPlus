#' anameth_ui
#'
#' @description The function to create the ui module for the analysis method
anameth_ui <- function(id){
  ns <- NS(id)

  tabPanel("7. Set the analysis method",

           tags$style("
             .edit-modal .modal-dialog {
               width: 90% !important;   /* Expand up to 90% of the screen width */
               max-width: 2000px;       /* Set maximum width */
             }
             .edit-modal .modal-body {
               max-height: 80vh;        /* Limit height to 80% of the viewport */
               overflow-y: auto;        /* Enable vertical scrolling if needed */
             }
           "),
           titlePanel("7. Set the analysis method"),
           br(),
           actionButton(ns("act_anameth_new"), "New"),
           actionButton(ns("act_anameth_copy"), "Copy"),
           actionButton(ns("act_anameth_edit"), "Edit"),
           actionButton(ns("act_anameth_delete"), "Delete"),
           br(),
           br(),
           dataTableOutput(ns("display_info_amameth"))
  )
}

#' anameth_server
#'
#' @description The function to create the server module for the analysis method
anameth_server <- function(id, path_pack){
  moduleServer(
    id,
    function(input, output, session){

      # --- Define the initial reactive values --- #

      rv_anameth <- reactiveValues()
      rv_anameth$dt_info_anameth <- create_dt_info_anameth(path_pack)
      rv_anameth$dt_anameth_edit <- NULL
      rv_anameth$anameth_edit <- NULL
      rv_anameth$list_dt_average_sr_origin <- NULL
      rv_anameth$list_dt_average_sr <- NULL
      rv_anameth$dt_average_sr <- NULL
      rv_anameth$current_mk <- NULL
      rv_anameth$alart_anameth_edit <- character(0)

      # --- Define a new analysis method --- #

      observeEvent(input$act_anameth_new, {
        kit_all <- list.files(paste0(path_pack, "/extdata/kit"))
        kit_all <- gsub(".{4}$", "", kit_all)

        stutter_all <- list.files(paste0(path_pack, "/extdata/stutter"))
        stutter_all <- gsub(".{6}$", "", stutter_all)

        lsae_all <- list.files(paste0(path_pack, "/extdata/lsae"))
        lsae_all <- gsub(".{6}$", "", lsae_all)

        gamma_all <- list.files(paste0(path_pack, "/extdata/gamma"))
        gamma_all <- gsub(".{6}$", "", gamma_all)

        dropin_all <- list.files(paste0(path_pack, "/extdata/dropin"))
        dropin_all <- gsub(".{6}$", "", dropin_all)

        showModal(modalDialog(
          title = "Define a new analysis method",
          textInput(session$ns("anameth_name"), label = "Enter a name of the analysis method"),
          br(),
          selectInput(session$ns("kit_name"), label = "Select kit", choices = kit_all, selected = kit_all[1], selectize = FALSE),
          br(),
          selectInput(session$ns("stutter_name"), label = "Select data for the stutter ratio", choices = stutter_all, selected = stutter_all[1], selectize = FALSE),
          br(),
          selectInput(session$ns("lsae_name"), label = "Select data for the locus-specific amplification efficiencies", choices = lsae_all, selected = lsae_all[1], selectize = FALSE),
          br(),
          selectInput(session$ns("gamma_name"), label = "Select data for the gamma scale parameter", choices = gamma_all, selected = gamma_all[1], selectize = FALSE),
          br(),
          selectInput(session$ns("dropin_name"), label = "Select data for the drop-in", choices = dropin_all, selected = dropin_all[1], selectize = FALSE),
          footer = tagList(
            actionButton(session$ns("act_add_anameth_new"), "Add the analysis method"),
            modalButton("Close")
          )
        ))
      })

      observeEvent(input$anameth_name, {
        anameth_name <- input$anameth_name
        anameth_all <- list.files(paste0(path_pack, "/extdata/anameth"))
        anameth_all <- gsub(".{4}$", "", anameth_all)

        if(is.element(anameth_name, anameth_all)){
          showFeedbackDanger(inputId = "anameth_name", text = "A file with the same name already exists.")
          disable("act_add_anameth_new")
        }else if(anameth_name == ""){
          hideFeedback("anameth_name")
          disable("act_add_anameth_new")
        }else{
          hideFeedback("anameth_name")
          enable("act_add_anameth_new")
        }
      })

      observeEvent(input$act_add_anameth_new, {
        anameth_name <- input$anameth_name
        kit_name <- input$kit_name
        stutter_name <- input$stutter_name
        lsae_name <- input$lsae_name
        gamma_name <- input$gamma_name
        dropin_name <- input$dropin_name

        if(all(!is.null(kit_name), !is.null(stutter_name), !is.null(lsae_name), !is.null(gamma_name), !is.null(dropin_name))){
          load(paste0(path_pack, "/extdata/stutter/", stutter_name, ".RData"))
          load(paste0(path_pack, "/extdata/lsae/", lsae_name, ".RData"))
          load(paste0(path_pack, "/extdata/gamma/", gamma_name, ".RData"))
          load(paste0(path_pack, "/extdata/dropin/", dropin_name, ".RData"))

          if(kit_name == data_est_sr$name_kit_stutter && kit_name == data_est_lsae$name_kit_lsae && kit_name == data_est_gamma$name_kit_gamma && kit_name == data_est_dropin$name_kit_dropin){
            dt_tmp <- data.table("Analysis method" = anameth_name,
                                 "Kit" = kit_name,
                                 "Stutter data" = stutter_name,
                                 "LSAE data" = lsae_name,
                                 "Gamma data" = gamma_name,
                                 "Drop-in data" = dropin_name)
            rv_anameth$dt_info_anameth <- rbind(rv_anameth$dt_info_anameth, dt_tmp)

            default_dt_anameth <- create_default_dt_anameth(path_pack, kit_name, stutter_name, lsae_name, gamma_name, dropin_name)
            write.csv(default_dt_anameth, paste0(path_pack, "/extdata/anameth/", anameth_name, ".csv"), row.names = FALSE)

            removeModal()
          }else{
            removeModal()
            showModal(modalDialog(title = "Error", "Kit names are inconsistent across datasets!", easyClose = TRUE, footer = NULL))
          }
        }
      })

      # --- Copy an analysis method --- #

      observeEvent(input$act_anameth_copy, {
        index_select <- input$display_info_amameth_rows_selected
        if(length(index_select) == 1){
          dt_info_anameth <- rv_anameth$dt_info_anameth
          anameth_copy <- paste0(dt_info_anameth[["Analysis method"]][index_select], "_copy")

          showModal(modalDialog(
            title = "Copy the analysis method",
            textInput(session$ns("name_copy_anameth"), label = "Enter the name of the analysis method", value = anameth_copy),
            footer = tagList(
              actionButton(session$ns("act_anameth_copy_save"), "Save"),
              modalButton("Cancel")
            )
          ))
        }else{
          showModal(modalDialog(title = "Error", "Select an analysis method!", easyClose = TRUE, footer = NULL))
        }
      })

      # Check the name of the analysis method
      observeEvent(input$name_copy_anameth, {
        name_copy_anameth <- input$name_copy_anameth
        dt_info_anameth <- rv_anameth$dt_info_anameth

        if(name_copy_anameth == ""){
          showFeedbackDanger(inputId = "name_copy_anameth", text = "This field cannot be empty.")
          disable("act_anameth_copy_save")
        }else if(is.element(name_copy_anameth, dt_info_anameth[["Analysis method"]])){
          showFeedbackDanger(inputId = "name_copy_anameth", text = "This name is already registered.")
          disable("act_anameth_copy_save")
        }else{
          hideFeedback("name_copy_anameth")
          enable("act_anameth_copy_save")
        }
      })

      observeEvent(input$act_anameth_copy_save, {
        index_select <- input$display_info_amameth_rows_selected
        name_copy_anameth <- input$name_copy_anameth
        dt_info_anameth <- rv_anameth$dt_info_anameth

        # Update the data.table for information on the analysis method
        dt_info_anameth_copy <- copy(dt_info_anameth)
        dt_info_anameth_copy <- dt_info_anameth_copy[index_select, ]
        anameth_select <- dt_info_anameth_copy[["Analysis method"]]
        dt_info_anameth_copy[["Analysis method"]] <- name_copy_anameth
        dt_info_anameth <- rbind(dt_info_anameth, dt_info_anameth_copy)
        rv_anameth$dt_info_anameth <- dt_info_anameth

        # Make a copy of the selected analysis method
        dt_anameth_copy <- fread(paste0(path_pack, "/extdata/anameth/", anameth_select, ".csv"))
        write.csv(dt_anameth_copy, paste0(path_pack, "/extdata/anameth/", name_copy_anameth, ".csv"), row.names = FALSE)

        removeModal()
      })

      # --- Edit an analysis method --- #

      observeEvent(input$act_anameth_edit, {
        index_select <- input$display_info_amameth_rows_selected
        if(length(index_select) == 1){

          # Load the selected analysis method
          dt_info_anameth <- rv_anameth$dt_info_anameth
          rv_anameth$anameth_edit <- anameth_edit <- dt_info_anameth[["Analysis method"]][index_select]
          rv_anameth$dt_anameth_edit <- dt_anameth_edit <- fread(paste0(path_pack, "/extdata/anameth/", anameth_edit, ".csv"))

          markers_auto <- unique(dt_anameth_edit$Subcategory[dt_anameth_edit$Category == "SR"])
          rv_anameth$current_mk <- markers_auto[1]
          list_dt_average_sr <- create_list_dt_average_sr(dt_anameth_edit, markers_auto)
          rv_anameth$list_dt_average_sr_origin <- rv_anameth$list_dt_average_sr <- list_dt_average_sr
          rv_anameth$dt_average_sr <- list_dt_average_sr[[1]]

          # Create UI
          showModal(modalDialog(
            title = "Edit the analysis method",
            uiOutput(session$ns("alart_edit_anameth")),
            tabsetPanel(
              id = session$ns("tabset_edit_anameth"),
              type = "tabs",
              tabPanel(
                "General",
                br(),
                dataTableOutput(session$ns("dt_edit_general"))
              ),
              tabPanel(
                "AT",
                br(),
                rHandsontableOutput(session$ns("dt_edit_at"))
              ),
              tabPanel(
                "Genotype filter",
                br(),
                rHandsontableOutput(session$ns("dt_edit_gtfltr"))
              ),
              tabPanel(
                "Degradation",
                br(),
                rHandsontableOutput(session$ns("dt_edit_deg"))
              ),
              tabPanel(
                "LSAE",
                br(),
                rHandsontableOutput(session$ns("dt_edit_lsae"))
              ),
              tabPanel(
                "Peak height",
                br(),
                rHandsontableOutput(session$ns("dt_edit_gamma"))
              ),
              tabPanel(
                "Drop-in",
                br(),
                rHandsontableOutput(session$ns("dt_edit_dropin"))
              ),
              tabPanel(
                "SR",
                br(),
                splitLayout(
                  cellWidths = c("300px", "400px", "400px"),
                  wellPanel(
                    selectInput(session$ns("select_mk"), label = "Marker", choices = markers_auto, selected = markers_auto[1], selectize = FALSE)
                  ),
                  wellPanel(
                    fluidRow(
                      column(
                        9,
                        numericInput(session$ns("add_al"), label = "Add an allele", value = NULL)
                      ),
                      column(
                        3,
                        br(),
                        actionButton(session$ns("act_add_al"), label = "Add")
                      )
                    )
                  ),
                  wellPanel(
                    fluidRow(
                      column(
                        9,
                        uiOutput(session$ns("select_del_al"))
                      ),
                      column(
                        3,
                        br(),
                        actionButton(session$ns("act_del_al"), label = "Delete")
                      )
                    )
                  )
                ),
                br(),
                tags$style(HTML("
                         .handsontable .htCore td {
                           vertical-align: middle !important;
                           text-align: center !important;
                         }"
                )),
                rHandsontableOutput(session$ns("dt_edit_sr"))
              )
            ),
            footer = tagList(
              actionButton(session$ns("act_anameth_edit_save"), "Save"),
              actionButton(session$ns("act_anameth_edit_cancel"), "Cancel")
            )
          ) %>% tagAppendAttributes(class = "edit-modal")
          )
        }else{
          showModal(modalDialog(title = "Error", "Select an analysis method!", easyClose = TRUE, footer = NULL))
        }
      })

      observeEvent(input$act_anameth_edit_cancel, {
        # Reset reactive values for editing the analysis method
        rv_anameth$dt_anameth_edit <- NULL
        rv_anameth$anameth_edit <- NULL
        rv_anameth$list_dt_average_sr_origin <- NULL
        rv_anameth$list_dt_average_sr <- NULL
        rv_anameth$dt_average_sr <- NULL
        rv_anameth$current_mk <- NULL
        rv_anameth$alart_anameth_edit <- character(0)

        removeModal()
      })

      # UI for alart message
      output$alart_edit_anameth <- renderUI({
        HTML(paste0('<span style="color: red;">', rv_anameth$alart_anameth_edit, '</span>'))
      })

      # UI for general information
      output$dt_edit_general <- renderDataTable({
        dt_anameth_edit <- rv_anameth$dt_anameth_edit
        anameth_edit <- rv_anameth$anameth_edit

        if(!is.null(dt_anameth_edit)){
          dt_edit_general <- data.table(Category = c("Analysis method", "Kit", "Stutter data", "LSAE data", "Gamma data", "Drop-in data"),
                                        Name = c(anameth_edit,
                                                 dt_anameth_edit$Name[dt_anameth_edit$Category == "Kit"],
                                                 dt_anameth_edit$Name[dt_anameth_edit$Category == "Stutter data"],
                                                 dt_anameth_edit$Name[dt_anameth_edit$Category == "LSAE data"],
                                                 dt_anameth_edit$Name[dt_anameth_edit$Category == "Gamma data"],
                                                 dt_anameth_edit$Name[dt_anameth_edit$Category == "Drop-in data"]
                                        ))

          datatable(
            dt_edit_general,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_edit_general), ordering = FALSE, autowidth = TRUE),
            rownames = FALSE,
            colnames = NULL,
            escape = FALSE
          )
        }
      })

      # UI for editing AT
      output$dt_edit_at <- renderRHandsontable({
        dt_anameth_edit <- rv_anameth$dt_anameth_edit

        if(!is.null(dt_anameth_edit)){
          dt_edit_at <- data.table(Marker = dt_anameth_edit$Subcategory[dt_anameth_edit$Category == "AT"],
                                   AT = dt_anameth_edit$Value[dt_anameth_edit$Category == "AT"])

          rhandsontable(dt_edit_at, rowHeaders = FALSE, height = 550) %>%
            hot_table(fillHandle = FALSE, rowHeights = 30) %>%
            hot_col("Marker", type = "text", readOnly = TRUE, width = 150) %>%
            hot_col("AT", type = "numeric", format = "0", width = 150)
        }
      })

      # UI for editing genotype filter
      output$dt_edit_gtfltr <- renderRHandsontable({
        dt_anameth_edit <- rv_anameth$dt_anameth_edit

        if(!is.null(dt_anameth_edit)){
          dt_edit_gtfltr <- data.table(Filter = dt_anameth_edit$Subcategory[dt_anameth_edit$Category == "Genotype filter"],
                                       Value = dt_anameth_edit$Value[dt_anameth_edit$Category == "Genotype filter"])

          rhandsontable(dt_edit_gtfltr, rowHeaders = FALSE) %>%
            hot_table(fillHandle = FALSE, rowHeights = 30) %>%
            hot_col("Filter", type = "text", readOnly = TRUE, width = 250) %>%
            hot_col("Value", type = "numeric", format = "0.00", width = 150)
        }
      })

      output$dt_edit_deg <- renderRHandsontable({
        dt_anameth_edit <- rv_anameth$dt_anameth_edit

        if(!is.null(dt_anameth_edit)){
          dt_edit_deg <- data.table(Filter = dt_anameth_edit$Subcategory[dt_anameth_edit$Category == "Degradation"],
                                    Value = dt_anameth_edit$Value[dt_anameth_edit$Category == "Degradation"])

          rhandsontable(dt_edit_deg, rowHeaders = FALSE) %>%
            hot_table(fillHandle = FALSE, rowHeights = 30) %>%
            hot_col("Filter", type = "text", readOnly = TRUE, width = 250) %>%
            hot_col("Value", type = "numeric", format = "0", width = 150)
        }
      })

      # UI for editing LSAE
      output$dt_edit_lsae <- renderRHandsontable({
        dt_anameth_edit <- rv_anameth$dt_anameth_edit

        if(!is.null(dt_anameth_edit)){
          dt_edit_lsae <- data.table(Marker = dt_anameth_edit$Subcategory[dt_anameth_edit$Category == "LSAE"],
                                     LSAE = as.character(dt_anameth_edit$Value[dt_anameth_edit$Category == "LSAE"]))

          rhandsontable(dt_edit_lsae, rowHeaders = FALSE, height = 550) %>%
            hot_table(fillHandle = FALSE, rowHeights = 30) %>%
            hot_col("Marker", type = "text", readOnly = TRUE, width = 150) %>%
            hot_col("LSAE", type = "text", width = 150,
                    renderer = "
                    function (instance, td, row, col, prop, value, cellProperties) {
                      Handsontable.renderers.TextRenderer.apply(this, arguments);
                      if(isNaN(parseFloat(value))) {
                        td.style.background = '#fdd';
                      } else {
                        td.style.background = '';
                      }
                    }
                  ")
        }
      })

      # UI for editing gamma parameters
      output$dt_edit_gamma <- renderRHandsontable({
        dt_anameth_edit <- rv_anameth$dt_anameth_edit

        if(!is.null(dt_anameth_edit)){
          dt_edit_gamma <- data.table(Parameter = dt_anameth_edit$Subcategory[dt_anameth_edit$Category == "Peak height"],
                                      Value = as.character(dt_anameth_edit$Value[dt_anameth_edit$Category == "Peak height"]))

          rhandsontable(dt_edit_gamma, rowHeaders = FALSE) %>%
            hot_table(fillHandle = FALSE, rowHeights = 30) %>%
            hot_col("Parameter", type = "text", readOnly = TRUE, width = 150) %>%
            hot_col("Value", type = "text", width = 150,
                    renderer = "
                    function (instance, td, row, col, prop, value, cellProperties) {
                      Handsontable.renderers.TextRenderer.apply(this, arguments);
                      if(isNaN(parseFloat(value))) {
                        td.style.background = '#fdd';
                      } else {
                        td.style.background = '';
                      }
                    }
                  ")
        }
      })

      # UI for editing drop-in parameters
      output$dt_edit_dropin <- renderRHandsontable({
        dt_anameth_edit <- rv_anameth$dt_anameth_edit

        if(!is.null(dt_anameth_edit)){
          dt_edit_dropin <- data.table(Parameter = dt_anameth_edit$Subcategory[dt_anameth_edit$Category == "Drop-in"],
                                       Value = as.character(dt_anameth_edit$Value[dt_anameth_edit$Category == "Drop-in"]))

          rhandsontable(dt_edit_dropin, rowHeaders = FALSE) %>%
            hot_table(fillHandle = FALSE, rowHeights = 30) %>%
            hot_col("Parameter", type = "text", readOnly = TRUE, width = 150) %>%
            hot_col("Value", type = "text", width = 150,
                    renderer = "
                    function (instance, td, row, col, prop, value, cellProperties) {
                      Handsontable.renderers.TextRenderer.apply(this, arguments);
                      if(isNaN(parseFloat(value))) {
                        td.style.background = '#fdd';
                      } else {
                        td.style.background = '';
                      }
                    }
                  ")
        }
      })

      # UI for editing SRs
      output$select_del_al <- renderUI({
        current_alleles <- rv_anameth$dt_average_sr$Allele
        selectInput(session$ns("del_al"), label = "Delete an allele", choices = current_alleles, selected = current_alleles[1], selectize = FALSE)
      })

      output$dt_edit_sr <- renderRHandsontable({
        dt_average_sr <- rv_anameth$dt_average_sr

        if(!is.null(dt_average_sr)){
          dt_average_sr[["Average BSR"]] <- as.character(dt_average_sr[["Average BSR"]])
          dt_average_sr[["Average FSR"]] <- as.character(dt_average_sr[["Average FSR"]])
          dt_average_sr[["Average DSR"]] <- as.character(dt_average_sr[["Average DSR"]])
          dt_average_sr[["Average M2SR"]] <- as.character(dt_average_sr[["Average M2SR"]])

          rhandsontable(dt_average_sr, rowHeaders = FALSE, height = 500) %>%
            hot_table(fillHandle = FALSE, rowHeights = 30) %>%
            hot_col(c("Allele"), type = "text", readOnly = TRUE, width = 50) %>%
            hot_col(c("Corrected allele (BSR)", "Corrected allele (FSR)", "Corrected allele (DSR)", "Corrected allele (M2SR)"), type = "numeric", format = "0.000", readOnly = TRUE, width = 175) %>%
            hot_col(c("Average BSR", "Average FSR", "Average DSR", "Average M2SR"), type = "text", width = 120,
                    renderer = "
                    function (instance, td, row, col, prop, value, cellProperties) {
                      Handsontable.renderers.TextRenderer.apply(this, arguments);
                      if(isNaN(parseFloat(value))) {
                        td.style.background = '#fdd';
                      } else {
                        td.style.background = '';
                      }
                    }
                  ")
        }
      })

      observeEvent(input$select_mk, {
        select_mk <- input$select_mk
        current_mk <- rv_anameth$current_mk

        # Update edited information
        if(!is.null(input$dt_edit_sr) && !is.null(select_mk) && !is.null(current_mk) && select_mk != rv_anameth$current_mk){
          dt_edit_sr <- hot_to_r(input$dt_edit_sr)
          options(warn = -1)
          dt_edit_sr[["Average BSR"]] <- as.numeric(dt_edit_sr[["Average BSR"]])
          dt_edit_sr[["Average FSR"]] <- as.numeric(dt_edit_sr[["Average FSR"]])
          dt_edit_sr[["Average DSR"]] <- as.numeric(dt_edit_sr[["Average DSR"]])
          dt_edit_sr[["Average M2SR"]] <- as.numeric(dt_edit_sr[["Average M2SR"]])
          options(warn = 0)

          # Change BSR = NA to BSR = 0
          tmp_bsr <- dt_edit_sr[["Average BSR"]]
          tmp_bsr[is.na(tmp_bsr)] <- 0
          dt_edit_sr[["Average BSR"]] <- tmp_bsr

          # Change FSR = NA to FSR = 0
          tmp_fsr <- dt_edit_sr[["Average FSR"]]
          tmp_fsr[is.na(tmp_fsr)] <- 0
          dt_edit_sr[["Average FSR"]] <- tmp_fsr

          # Change DSR = NA to DSR = 0
          tmp_dsr <- dt_edit_sr[["Average DSR"]]
          tmp_dsr[is.na(tmp_dsr)] <- 0
          dt_edit_sr[["Average DSR"]] <- tmp_dsr

          # Change M2SR = NA to M2SR = 0
          tmp_m2sr <- dt_edit_sr[["Average M2SR"]]
          tmp_m2sr[is.na(tmp_m2sr)] <- 0
          dt_edit_sr[["Average M2SR"]] <- tmp_m2sr

          rv_anameth$list_dt_average_sr[[current_mk]] <- dt_edit_sr

          rv_anameth$current_mk <- select_mk
          rv_anameth$dt_average_sr <- rv_anameth$list_dt_average_sr[[select_mk]]
        }
      })

      observeEvent(input$act_add_al, {
        add_al <- input$add_al
        if(!is.na(add_al)){
          dt_average_sr <- rv_anameth$dt_average_sr
          dt_average_sr <- rbind(dt_average_sr, list(as.character(add_al), NA, 0, NA, 0, NA, 0, NA, 0))

          # Sort allele
          tmp_al <- as.numeric(dt_average_sr$Allele)
          tmp_al_sort <- sort(tmp_al)
          dt_average_sr <- dt_average_sr[match(tmp_al_sort, tmp_al)]

          rv_anameth$dt_average_sr <- dt_average_sr
        }
      })

      observeEvent(input$act_del_al, {
        del_al <- input$del_al
        dt_average_sr <- rv_anameth$dt_average_sr
        current_alleles <- dt_average_sr$Allele
        rv_anameth$dt_average_sr <- dt_average_sr[- which(current_alleles == del_al)]
      })

      observeEvent(input$act_anameth_edit_save, {
        dt_anameth_edit <- rv_anameth$dt_anameth_edit
        anameth_edit <- rv_anameth$anameth_edit
        index_select <- input$display_info_amameth_rows_selected
        dt_info_anameth <- rv_anameth$dt_info_anameth
        bool_alart_anameth <- FALSE
        alart_anameth_edit <- rv_anameth$alart_anameth_edit

        if(!is.null(input$dt_edit_at)){
          dt_edit_at <- hot_to_r(input$dt_edit_at)
          options(warn = -1)
          dt_edit_at$AT <- as.numeric(dt_edit_at$AT)
          options(warn = 0)

          if(any(is.na(dt_edit_at$AT))){
            bool_alart_anameth <- TRUE
          }else{
            dt_anameth_edit$Value[dt_anameth_edit$Category == "AT"] <- dt_edit_at$AT
          }
        }

        if(!is.null(input$dt_edit_gtfltr)){
          dt_edit_gtfltr <- hot_to_r(input$dt_edit_gtfltr)
          options(warn = -1)
          dt_edit_gtfltr$Value <- as.numeric(dt_edit_gtfltr$Value)
          options(warn = 0)

          if(any(is.na(dt_edit_gtfltr$Value))){
            bool_alart_anameth <- TRUE
          }else{
            dt_anameth_edit$Value[dt_anameth_edit$Category == "Genotype filter"] <- dt_edit_gtfltr$Value
          }
        }

        if(!is.null(input$dt_edit_deg)){
          dt_edit_deg <- hot_to_r(input$dt_edit_deg)
          options(warn = -1)
          dt_edit_deg$Value <- as.numeric(dt_edit_deg$Value)
          options(warn = 0)

          if(any(is.na(dt_edit_deg$Value))){
            bool_alart_anameth <- TRUE
          }else{
            dt_anameth_edit$Value[dt_anameth_edit$Category == "Degradation"] <- dt_edit_deg$Value
          }
        }

        if(!is.null(input$dt_edit_lsae)){
          dt_edit_lsae <- hot_to_r(input$dt_edit_lsae)
          options(warn = -1)
          dt_edit_lsae$LSAE <- as.numeric(dt_edit_lsae$LSAE)
          options(warn = 0)

          if(any(is.na(dt_edit_lsae$LSAE))){
            bool_alart_anameth <- TRUE
          }else{
            # Mark as edited when the LSAE parameters or gamma parameters are changed.
            if(any(round(dt_anameth_edit$Value[dt_anameth_edit$Category == "LSAE"], 6) != round(dt_edit_lsae$LSAE, 6))){
              name_lsae <- dt_anameth_edit$Name[dt_anameth_edit$Category == "LSAE data"]

              if(substring(name_lsae, nchar(name_lsae) - 8) != " (edited)"){
                dt_anameth_edit$Name[dt_anameth_edit$Category == "LSAE data"] <- paste0(name_lsae, " (edited)")
                dt_info_anameth[["LSAE data"]][index_select] <- paste0(name_lsae, " (edited)")
              }
            }

            dt_anameth_edit$Value[dt_anameth_edit$Category == "LSAE"] <- dt_edit_lsae$LSAE
          }
        }

        if(!is.null(input$dt_edit_gamma)){
          dt_edit_gamma <- hot_to_r(input$dt_edit_gamma)
          options(warn = -1)
          dt_edit_gamma$Value <- as.numeric(dt_edit_gamma$Value)
          options(warn = 0)

          if(any(is.na(dt_edit_gamma$Value))){
            bool_alart_anameth <- TRUE
          }else{
            # Mark as edited when the LSAE parameters or gamma parameters are changed.
            if(any(round(dt_anameth_edit$Value[dt_anameth_edit$Category == "Peak height"], 6) != round(dt_edit_gamma$Value, 6))){
              name_gamma <- dt_anameth_edit$Name[dt_anameth_edit$Category == "Gamma data"]

              if(substring(name_gamma, nchar(name_gamma) - 8) != " (edited)"){
                dt_anameth_edit$Name[dt_anameth_edit$Category == "Gamma data"] <- paste0(name_gamma, " (edited)")
                dt_info_anameth[["Gamma data"]][index_select] <- paste0(name_gamma, " (edited)")
              }
            }

            dt_anameth_edit$Value[dt_anameth_edit$Category == "Peak height"] <- dt_edit_gamma$Value
          }
        }

        if(!is.null(input$dt_edit_dropin)){
          dt_edit_dropin <- hot_to_r(input$dt_edit_dropin)
          options(warn = -1)
          dt_edit_dropin$Value <- as.numeric(dt_edit_dropin$Value)
          options(warn = 0)

          if(any(is.na(dt_edit_dropin$Value))){
            bool_alart_anameth <- TRUE
          }else{
            # Mark as edited when the drop-in parameters are changed.
            if(any(round(dt_anameth_edit$Value[dt_anameth_edit$Category == "Drop-in"], 6) != round(dt_edit_dropin$Value, 6))){
              name_dropin <- dt_anameth_edit$Name[dt_anameth_edit$Category == "Drop-in data"]

              if(substring(name_dropin, nchar(name_dropin) - 8) != " (edited)"){
                dt_anameth_edit$Name[dt_anameth_edit$Category == "Drop-in data"] <- paste0(name_dropin, " (edited)")
                dt_info_anameth[["Drop-in data"]][index_select] <- paste0(name_dropin, " (edited)")
              }
            }

            dt_anameth_edit$Value[dt_anameth_edit$Category == "Drop-in"] <- dt_edit_dropin$Value
          }
        }

        if(!is.null(input$dt_edit_sr)){
          dt_edit_sr <- hot_to_r(input$dt_edit_sr)
          options(warn = -1)
          dt_edit_sr[["Average BSR"]] <- as.numeric(dt_edit_sr[["Average BSR"]])
          dt_edit_sr[["Average FSR"]] <- as.numeric(dt_edit_sr[["Average FSR"]])
          dt_edit_sr[["Average DSR"]] <- as.numeric(dt_edit_sr[["Average DSR"]])
          dt_edit_sr[["Average M2SR"]] <- as.numeric(dt_edit_sr[["Average M2SR"]])
          options(warn = 0)

          # Change BSR = NA to BSR = 0
          tmp_bsr <- dt_edit_sr[["Average BSR"]]
          tmp_bsr[is.na(tmp_bsr)] <- 0
          dt_edit_sr[["Average BSR"]] <- tmp_bsr

          # Change FSR = NA to FSR = 0
          tmp_fsr <- dt_edit_sr[["Average FSR"]]
          tmp_fsr[is.na(tmp_fsr)] <- 0
          dt_edit_sr[["Average FSR"]] <- tmp_fsr

          # Change DSR = NA to DSR = 0
          tmp_dsr <- dt_edit_sr[["Average DSR"]]
          tmp_dsr[is.na(tmp_dsr)] <- 0
          dt_edit_sr[["Average DSR"]] <- tmp_dsr

          # Change M2SR = NA to M2SR = 0
          tmp_m2sr <- dt_edit_sr[["Average M2SR"]]
          tmp_m2sr[is.na(tmp_m2sr)] <- 0
          dt_edit_sr[["Average M2SR"]] <- tmp_m2sr

          rv_anameth$list_dt_average_sr[[rv_anameth$current_mk]] <- dt_edit_sr

          list_dt_average_sr <- rv_anameth$list_dt_average_sr
          list_dt_average_sr_origin <- rv_anameth$list_dt_average_sr_origin

          # Mark as edited when the SRs are changed.
          bool_edited <- compare_list_dt_average_sr(list_dt_average_sr_origin, list_dt_average_sr)
          if(bool_edited){
            name_stutter <- dt_anameth_edit$Name[dt_anameth_edit$Category == "Stutter data"]

            if(substring(name_stutter, nchar(name_stutter) - 8) != " (edited)"){
              dt_anameth_edit$Name[dt_anameth_edit$Category == "Stutter data"] <- paste0(name_stutter, " (edited)")
              dt_info_anameth[["Stutter data"]][index_select] <- paste0(name_stutter, " (edited)")
            }
          }

          dt_anameth_sr <- create_dt_anameth_sr(list_dt_average_sr, names(list_dt_average_sr))
          dt_anameth_edit <- rbind(dt_anameth_edit[!is.element(dt_anameth_edit$Category, "SR")], dt_anameth_sr)
        }

        rv_anameth$alart_anameth_edit <- alart_anameth_edit

        if(bool_alart_anameth){
          rv_anameth$alart_anameth_edit <- "Error: Some cells are not numeric!"
        }else{
          write.csv(dt_anameth_edit, paste0(path_pack, "/extdata/anameth/", anameth_edit, ".csv"), row.names = FALSE)
          rv_anameth$dt_info_anameth <- dt_info_anameth

          # Reset reactive values for editing the analysis method
          rv_anameth$dt_anameth_edit <- NULL
          rv_anameth$anameth_edit <- NULL
          rv_anameth$list_dt_average_sr_origin <- NULL
          rv_anameth$list_dt_average_sr <- NULL
          rv_anameth$dt_average_sr <- NULL
          rv_anameth$current_mk <- NULL
          rv_anameth$alart_anameth_edit <- character(0)

          removeModal()
        }
      })

      # --- Delete an analysis method --- #

      observeEvent(input$act_anameth_delete, {
        index_select <- input$display_info_amameth_rows_selected
        if(length(index_select) == 1){
          showModal(modalDialog(
            title = "Message",
            "Are you sure you want to delete this method?",
            footer = tagList(
              actionButton(session$ns("act_anameth_delete_done"), "Yes"),
              modalButton("No")
            )
          ))
        }else{
          showModal(modalDialog(title = "Error", "Select an analysis method!", easyClose = TRUE, footer = NULL))
        }
      })

      observeEvent(input$act_anameth_delete_done, {
        index_select <- input$display_info_amameth_rows_selected
        dt_info_anameth <- rv_anameth$dt_info_anameth
        file.remove(paste0(path_pack, "/extdata/anameth/", dt_info_anameth[["Analysis method"]][index_select], ".csv"))
        rv_anameth$dt_info_anameth <- dt_info_anameth[-index_select, ]

        removeModal()
      })

      # --- Display information on the analysis method --- #

      output$display_info_amameth <- renderDataTable({
        dt_info_anameth <- rv_anameth$dt_info_anameth
        datatable(
          dt_info_anameth,
          selection = list(mode = "single", target = "row"),
          options = list(dom = "t", iDisplayLength = nrow(dt_info_anameth), ordering = FALSE, autoWidth = TRUE, scrollY = "650px", fixedHeader = TRUE),
          rownames = FALSE
        )
      })
    }
  )
}
