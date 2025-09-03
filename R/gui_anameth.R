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

#' create_dt_info_anameth
#'
#' @description The function to create the data.table for information on the analysis methods
#' @param path_pack The package path
create_dt_info_anameth <- function(path_pack){
  fn_anameth <- list.files(paste0(path_pack, "/extdata/anameth"))
  n_anameth <- length(fn_anameth)

  if(n_anameth != 0){
    info_kit <- info_stutter <- info_param <- info_dropin <- rep("", n_anameth)
    for(i in 1:n_anameth){
      dt_anameth_1 <- fread(paste0(path_pack, "/extdata/anameth/", fn_anameth[i]))
      info_kit[i] <- dt_anameth_1$Name[dt_anameth_1$Category == "Kit"]
      info_stutter[i] <- dt_anameth_1$Name[dt_anameth_1$Category == "Stutter data"]
      info_param[i] <- dt_anameth_1$Name[dt_anameth_1$Category == "Parameter data"]
      info_dropin[i] <- dt_anameth_1$Name[dt_anameth_1$Category == "Drop-in data"]
    }

    dt_info_anameth <- data.table("Analysis method" = gsub(".{4}$", "", fn_anameth),
                                  "Kit" = info_kit,
                                  "Stutter data" = info_stutter,
                                  "Parameter data" = info_param,
                                  "Drop-in data" = info_dropin)
    return(dt_info_anameth)
  }else{
    return(NULL)
  }
}

#' create_dt_anameth_sr
#'
#' @description The function to create the data.table for stutter ratios as a part of the analysis method
#' @param list_dt_average_sr A list of the data.table for the average stutter ratios
#' @param markers_auto A vector of the autosomal markers
create_dt_anameth_sr <- function(list_dt_average_sr, markers_auto){
  vec_nrow_dt_average_sr <- sapply(list_dt_average_sr, nrow)
  nrow_sr_all <- sum(vec_nrow_dt_average_sr)
  sr_subcategory <- rep("", nrow_sr_all)
  sr_allele <- rep("", nrow_sr_all)
  sr_alcor_bsr <- rep("", nrow_sr_all)
  sr_bsr <- rep("", nrow_sr_all)
  sr_alcor_fsr <- rep("", nrow_sr_all)
  sr_fsr <- rep("", nrow_sr_all)
  sr_alcor_dsr <- rep("", nrow_sr_all)
  sr_dsr <- rep("", nrow_sr_all)
  sr_alcor_m2sr <- rep("", nrow_sr_all)
  sr_m2sr <- rep("", nrow_sr_all)

  index_sr_start <- 1
  for(i in 1:n_mk){
    mk <- markers_auto[i]
    dt_average_sr <- list_dt_average_sr[[mk]]

    index_assign <- index_sr_start:(index_sr_start + vec_nrow_dt_average_sr[i] - 1)
    sr_subcategory[index_assign] <- mk
    sr_allele[index_assign] <- dt_average_sr[["Allele"]]
    sr_alcor_bsr[index_assign] <- dt_average_sr[["Corrected allele (BSR)"]]
    sr_bsr[index_assign] <- dt_average_sr[["Average BSR"]]
    sr_alcor_fsr[index_assign] <- dt_average_sr[["Corrected allele (FSR)"]]
    sr_fsr[index_assign] <- dt_average_sr[["Average FSR"]]
    sr_alcor_dsr[index_assign] <- dt_average_sr[["Corrected allele (DSR)"]]
    sr_dsr[index_assign] <- dt_average_sr[["Average DSR"]]
    sr_alcor_m2sr[index_assign] <- dt_average_sr[["Corrected allele (M2SR)"]]
    sr_m2sr[index_assign] <- dt_average_sr[["Average M2SR"]]

    index_sr_start <- index_sr_start + vec_nrow_dt_average_sr[i]
  }

  dt_anameth_sr <- data.table(Category = rep("SR", nrow_sr_all),
                              Subcategory = sr_subcategory,
                              Name = rep("", nrow_sr_all),
                              Value = rep(NA, nrow_sr_all),
                              Allele = sr_allele,
                              "Corrected allele (BSR)" = sr_alcor_bsr,
                              "Average BSR" = sr_bsr,
                              "Corrected allele (FSR)" = sr_alcor_fsr,
                              "Average FSR" = sr_fsr,
                              "Corrected allele (DSR)" = sr_alcor_dsr,
                              "Average DSR" = sr_dsr,
                              "Corrected allele (M2SR)" = sr_alcor_m2sr,
                              "Average M2SR" = sr_m2sr)

  return(dt_anameth_sr)
}

#' create_default_dt_anameth
#'
#' @description The function to create the default data.table for a new analysis method
#' @param path_pack The package path of KongohPlus
#' @param kit_name A kit name
#' @param stutter_name The name of the stutter data
#' @param param_name The name of the parameter data
#' @param dropin_name The name of the drop-in data
create_default_dt_anameth <- function(path_pack, kit_name, stutter_name, param_name, dropin_name){
  dt_kit <- fread(paste0(path_pack, "/extdata/kit/", kit_name, ".csv"))
  markers_auto <- unique(dt_kit[["Marker"]][dt_kit[["Chromosome"]] == "Autosomal"])
  n_mk <- length(markers_auto)

  load(paste0(path_pack, "/extdata/stutter/", stutter_name, ".RData"))
  load(paste0(path_pack, "/extdata/param/", param_name, ".RData"))
  load(paste0(path_pack, "/extdata/dropin/", dropin_name, ".RData"))

  default_dt_anameth <- data.table(Category = c("Kit", "Stutter data", "Parameter data", "Drop-in data",
                                                rep("AT", n_mk),
                                                rep("Genotype filter", 6),
                                                rep("LSAE", n_mk),
                                                rep("Peak height", 2),
                                                rep("Drop-in", 3)),
                                   Subcategory = c(rep("", 4),
                                                   markers_auto,
                                                   "Maximum BSR", "Maximum FSR", "Maximum DSR", "Maximum M2SR", "Maximum number of drop-in", "Maximum height of drop-in",
                                                   names(expect_ae),
                                                   "CV for allele", "CV for stutter",
                                                   "Frequency", "Shape parameter", "Rate parameter"),
                                   Name = c(kit_name, stutter_name, param_name, dropin_name, rep("", 2 * n_mk + 11)),
                                   Value = c(rep(NA, 4),
                                             rep(100, n_mk),
                                             0.3, 0.2, 0.05, 0.1, 2, 150,
                                             expect_ae,
                                             cv_allele, cv_stutter,
                                             freq_dropin, shape_di, rate_di
                                             ),
                                   Allele = rep(NA, 2 * n_mk + 15),
                                   "Corrected allele (BSR)" = rep(NA, 2 * n_mk + 15),
                                   "Average BSR" = rep(NA, 2 * n_mk + 15),
                                   "Corrected allele (FSR)" = rep(NA, 2 * n_mk + 15),
                                   "Average FSR" = rep(NA, 2 * n_mk + 15),
                                   "Corrected allele (DSR)" = rep(NA, 2 * n_mk + 15),
                                   "Average DSR" = rep(NA, 2 * n_mk + 15),
                                   "Corrected allele (M2SR)" = rep(NA, 2 * n_mk + 15),
                                   "Average M2SR" = rep(NA, 2 * n_mk + 15))

  default_dt_anameth_sr <- create_dt_anameth_sr(list_dt_average_sr, markers_auto)

  default_dt_anameth <- rbind(default_dt_anameth, default_dt_anameth_sr)

  return(default_dt_anameth)
}

#' compare_list_dt_average_sr
#'
#' @description The function to compare 'list_dt_average_sr_origin' (before editing) with 'list_dt_average_sr' (after editing)
#' @param list_dt_average_sr_origin A list of the data.table for average stutter ratios in each marker (before editing)
#' @param list_dt_average_sr A list of the data.table for average stutter ratios in each marker (after editing)
compare_list_dt_average_sr <- function(list_dt_average_sr_origin, list_dt_average_sr){
  markers <- names(list_dt_average_sr_origin)
  bool_edited <- FALSE

  for(i in 1:length(markers)){
    mk <- markers[i]
    dt_average_sr_origin <- list_dt_average_sr_origin[[mk]]
    dt_average_sr <- list_dt_average_sr[[mk]]

    if(!setequal(dt_average_sr_origin$Allele, dt_average_sr$Allele)){
      bool_edited <- TRUE
      break
    }else if(any(dt_average_sr_origin[["Average BSR"]] != dt_average_sr[["Average BSR"]])){
      bool_edited <- TRUE
      break
    }else if(any(dt_average_sr_origin[["Average FSR"]] != dt_average_sr[["Average FSR"]])){
      bool_edited <- TRUE
      break
    }else if(any(dt_average_sr_origin[["Average DSR"]] != dt_average_sr[["Average DSR"]])){
      bool_edited <- TRUE
      break
    }else if(any(dt_average_sr_origin[["Average M2SR"]] != dt_average_sr[["Average M2SR"]])){
      bool_edited <- TRUE
      break
    }
  }

  return(bool_edited)
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

      # --- Define a new analysis method --- #

      observeEvent(input$act_anameth_new, {
        kit_all <- list.files(paste0(path_pack, "/extdata/kit"))
        kit_all <- gsub(".{4}$", "", kit_all)

        stutter_all <- list.files(paste0(path_pack, "/extdata/stutter"))
        stutter_all <- gsub(".{6}$", "", stutter_all)

        param_all <- list.files(paste0(path_pack, "/extdata/param"))
        param_all <- gsub(".{6}$", "", param_all)

        dropin_all <- list.files(paste0(path_pack, "/extdata/dropin"))
        dropin_all <- gsub(".{6}$", "", dropin_all)

        showModal(modalDialog(
          title = "Define a new analysis method",
          textInput(session$ns("anameth_name"), label = "Enter a name of the analysis method"),
          br(),
          selectInput(session$ns("kit_name"), label = "Select kit", choices = kit_all, selected = kit_all[1], selectize = FALSE),
          br(),
          selectInput(session$ns("stutter_name"), label = "Select stutter data", choices = stutter_all, selected = stutter_all[1], selectize = FALSE),
          br(),
          selectInput(session$ns("param_name"), label = "Select parameter data", choices = param_all, selected = param_all[1], selectize = FALSE),
          br(),
          selectInput(session$ns("dropin_name"), label = "Select drop-in data", choices = dropin_all, selected = dropin_all[1], selectize = FALSE),
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
        param_name <- input$param_name
        dropin_name <- input$dropin_name

        if(all(!is.null(kit_name), !is.null(stutter_name), !is.null(param_name), !is.null(dropin_name))){
          load(paste0(path_pack, "/extdata/stutter/", stutter_name, ".RData"))
          load(paste0(path_pack, "/extdata/param/", param_name, ".RData"))
          load(paste0(path_pack, "/extdata/dropin/", dropin_name, ".RData"))

          if(kit_name == name_kit_stutter && kit_name == name_kit_param && kit_name == name_kit_dropin){
            dt_tmp <- data.table("Analysis method" = anameth_name,
                                 "Kit" = kit_name,
                                 "Stutter data" = stutter_name,
                                 "Parameter data" = param_name,
                                 "Drop-in data" = dropin_name)
            rv_anameth$dt_info_anameth <- rbind(rv_anameth$dt_info_anameth, dt_tmp)

            default_dt_anameth <- create_default_dt_anameth(path_pack, kit_name, stutter_name, param_name, dropin_name)
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
        if(isTruthy(index_select)){
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
        if(isTruthy(index_select)){

          # Load the selected analysis method
          dt_info_anameth <- rv_anameth$dt_info_anameth
          rv_anameth$anameth_edit <- anameth_edit <- dt_info_anameth[["Analysis method"]][index_select]
          rv_anameth$dt_anameth_edit <- dt_anameth_edit <- fread(paste0(path_pack, "/extdata/anameth/", anameth_edit, ".csv"))

          markers_auto <- unique(dt_anameth_edit$Subcategory[dt_anameth_edit$Category == "SR"])
          rv_anameth$current_mk <- markers_auto[1]
          n_mk <- length(markers_auto)
          list_dt_average_sr <- vector(mode = "list", length = n_mk)
          names(list_dt_average_sr) <- markers_auto
          for(i in 1:n_mk){
            mk <- markers_auto[i]
            index_target <- sort(intersect(which(dt_anameth_edit$Category == "SR"), which(dt_anameth_edit$Subcategory == mk)))

            dt_average_sr <- data.table("Allele" = dt_anameth_edit[["Allele"]][index_target],
                                        "Corrected allele (BSR)" = dt_anameth_edit[["Corrected allele (BSR)"]][index_target],
                                        "Average BSR" = dt_anameth_edit[["Average BSR"]][index_target],
                                        "Corrected allele (FSR)" = dt_anameth_edit[["Corrected allele (FSR)"]][index_target],
                                        "Average FSR" = dt_anameth_edit[["Average FSR"]][index_target],
                                        "Corrected allele (DSR)" = dt_anameth_edit[["Corrected allele (DSR)"]][index_target],
                                        "Average DSR" = dt_anameth_edit[["Average DSR"]][index_target],
                                        "Corrected allele (M2SR)" = dt_anameth_edit[["Corrected allele (M2SR)"]][index_target],
                                        "Average M2SR" = dt_anameth_edit[["Average M2SR"]][index_target])
            list_dt_average_sr[[i]] <- dt_average_sr
          }
          rv_anameth$list_dt_average_sr_origin <- rv_anameth$list_dt_average_sr <- list_dt_average_sr
          rv_anameth$dt_average_sr <- list_dt_average_sr[[1]]

          # Create UI
          showModal(modalDialog(
            title = "Edit the analysis method",
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

        removeModal()
      })

      # UI for general information
      output$dt_edit_general <- renderDataTable({
        dt_anameth_edit <- rv_anameth$dt_anameth_edit
        anameth_edit <- rv_anameth$anameth_edit

        if(!is.null(dt_anameth_edit)){
          dt_edit_general <- data.table(Category = c("Analysis method", "Kit", "Stutter data", "Parameter data", "Drop-in data"),
                                        Name = c(anameth_edit,
                                                 dt_anameth_edit$Name[dt_anameth_edit$Category == "Kit"],
                                                 dt_anameth_edit$Name[dt_anameth_edit$Category == "Stutter data"],
                                                 dt_anameth_edit$Name[dt_anameth_edit$Category == "Parameter data"],
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

      # UI for editing LSAE
      output$dt_edit_lsae <- renderRHandsontable({
        dt_anameth_edit <- rv_anameth$dt_anameth_edit

        if(!is.null(dt_anameth_edit)){
          dt_edit_lsae <- data.table(Marker = dt_anameth_edit$Subcategory[dt_anameth_edit$Category == "LSAE"],
                                     LSAE = dt_anameth_edit$Value[dt_anameth_edit$Category == "LSAE"])

          rhandsontable(dt_edit_lsae, rowHeaders = FALSE, height = 550) %>%
            hot_table(fillHandle = FALSE, rowHeights = 30) %>%
            hot_col("Marker", type = "text", readOnly = TRUE, width = 150) %>%
            hot_col("LSAE", type = "numeric", format = "0.000000", width = 150)
        }
      })

      # UI for editing gamma parameters
      output$dt_edit_gamma <- renderRHandsontable({
        dt_anameth_edit <- rv_anameth$dt_anameth_edit

        if(!is.null(dt_anameth_edit)){
          dt_edit_gamma <- data.table(Parameter = dt_anameth_edit$Subcategory[dt_anameth_edit$Category == "Peak height"],
                                      Value = dt_anameth_edit$Value[dt_anameth_edit$Category == "Peak height"])

          rhandsontable(dt_edit_gamma, rowHeaders = FALSE) %>%
            hot_table(fillHandle = FALSE, rowHeights = 30) %>%
            hot_col("Parameter", type = "text", readOnly = TRUE, width = 150) %>%
            hot_col("Value", type = "numeric", format = "0.000000", width = 150)
        }
      })

      # UI for editing drop-in parameters
      output$dt_edit_dropin <- renderRHandsontable({
        dt_anameth_edit <- rv_anameth$dt_anameth_edit

        if(!is.null(dt_anameth_edit)){
          dt_edit_dropin <- data.table(Parameter = dt_anameth_edit$Subcategory[dt_anameth_edit$Category == "Drop-in"],
                                       Value = dt_anameth_edit$Value[dt_anameth_edit$Category == "Drop-in"])

          rhandsontable(dt_edit_dropin, rowHeaders = FALSE) %>%
            hot_table(fillHandle = FALSE, rowHeights = 30) %>%
            hot_col("Parameter", type = "text", readOnly = TRUE, width = 150) %>%
            hot_col("Value", type = "numeric", format = "0.000000", width = 150)
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
          rhandsontable(dt_average_sr, rowHeaders = FALSE, height = 500) %>%
            hot_table(fillHandle = FALSE, rowHeights = 30) %>%
            hot_col(c("Allele"), type = "text", readOnly = TRUE, width = 50) %>%
            hot_col(c("Corrected allele (BSR)", "Corrected allele (FSR)", "Corrected allele (DSR)", "Corrected allele (M2SR)"), type = "numeric", format = "0.000", readOnly = TRUE, width = 175) %>%
            hot_col(c("Average BSR", "Average FSR", "Average DSR", "Average M2SR"), type = "numeric", format = "0.0000", width = 120)
        }
      })

      observeEvent(input$select_mk, {
        select_mk <- input$select_mk
        current_mk <- rv_anameth$current_mk

        # Update edited information
        if(!is.null(input$dt_edit_sr) && !is.null(select_mk) && !is.null(current_mk) && select_mk != rv_anameth$current_mk){
          dt_edit_sr <- hot_to_r(input$dt_edit_sr)

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

        if(!is.null(input$dt_edit_at)){
          dt_edit_at <- hot_to_r(input$dt_edit_at)
          dt_anameth_edit$Value[dt_anameth_edit$Category == "AT"] <- dt_edit_at$AT
        }

        if(!is.null(input$dt_edit_gtfltr)){
          dt_edit_gtfltr <- hot_to_r(input$dt_edit_gtfltr)
          dt_anameth_edit$Value[dt_anameth_edit$Category == "Genotype filter"] <- dt_edit_gtfltr$Value
        }

        if(!is.null(input$dt_edit_lsae)){
          dt_edit_lsae <- hot_to_r(input$dt_edit_lsae)

          # Mark as edited when the LSAE parameters or gamma parameters are changed.
          if(any(round(dt_anameth_edit$Value[dt_anameth_edit$Category == "LSAE"], 6) != round(dt_edit_lsae$LSAE, 6))){
            name_param <- dt_anameth_edit$Name[dt_anameth_edit$Category == "Parameter data"]

            if(substring(name_param, nchar(name_param) - 8) != " (edited)"){
              dt_anameth_edit$Name[dt_anameth_edit$Category == "Parameter data"] <- paste0(name_param, " (edited)")
              dt_info_anameth[["Parameter data"]][index_select] <- paste0(name_param, " (edited)")
            }
          }

          dt_anameth_edit$Value[dt_anameth_edit$Category == "LSAE"] <- dt_edit_lsae$LSAE
        }

        if(!is.null(input$dt_edit_gamma)){
          dt_edit_gamma <- hot_to_r(input$dt_edit_gamma)

          # Mark as edited when the LSAE parameters or gamma parameters are changed.
          if(any(round(dt_anameth_edit$Value[dt_anameth_edit$Category == "Peak height"], 6) != round(dt_edit_gamma$Value, 6))){
            name_param <- dt_anameth_edit$Name[dt_anameth_edit$Category == "Parameter data"]

            if(substring(name_param, nchar(name_param) - 8) != " (edited)"){
              dt_anameth_edit$Name[dt_anameth_edit$Category == "Parameter data"] <- paste0(name_param, " (edited)")
              dt_info_anameth[["Parameter data"]][index_select] <- paste0(name_param, " (edited)")
            }
          }

          dt_anameth_edit$Value[dt_anameth_edit$Category == "Peak height"] <- dt_edit_gamma$Value
        }

        if(!is.null(input$dt_edit_dropin)){
          dt_edit_dropin <- hot_to_r(input$dt_edit_dropin)

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

        if(!is.null(input$dt_edit_sr)){
          dt_edit_sr <- hot_to_r(input$dt_edit_sr)

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

        write.csv(dt_anameth_edit, paste0(path_pack, "/extdata/anameth/", anameth_edit, ".csv"), row.names = FALSE)
        rv_anameth$dt_info_anameth <- dt_info_anameth

        # Reset reactive values for editing the analysis method
        rv_anameth$dt_anameth_edit <- NULL
        rv_anameth$anameth_edit <- NULL
        rv_anameth$list_dt_average_sr_origin <- NULL
        rv_anameth$list_dt_average_sr <- NULL
        rv_anameth$dt_average_sr <- NULL
        rv_anameth$current_mk <- NULL

        removeModal()
      })

      # --- Delete an analysis method --- #

      observeEvent(input$act_anameth_delete, {
        index_select <- input$display_info_amameth_rows_selected
        if(isTruthy(index_select)){
          dt_info_anameth <- rv_anameth$dt_info_anameth
          file.remove(paste0(path_pack, "/extdata/anameth/", dt_info_anameth[["Analysis method"]][index_select], ".csv"))
          rv_anameth$dt_info_anameth <- dt_info_anameth[-index_select, ]
        }else{
          showModal(modalDialog(title = "Error", "Select an analysis method!", easyClose = TRUE, footer = NULL))
        }
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
