#' anameth_ui
#'
#' @description The function to create the ui module for the analysis method
anameth_ui <- function(id){
  ns <- NS(id)

  tabPanel("Analysis method",
           titlePanel("Analysis method"),
           br(),
           fluidRow(
             column(7,
                    tags$head(
                      tags$style(HTML("table {table-layout: fixed;")),
                      tags$style(".datatables .display {margin-left: 0;}")
                    ),
                    actionButton(ns("act_anameth_new"), "Define a new analysis method"),
                    actionButton(ns("act_anameth_delete"), "Delete"),
                    actionButton(ns("act_anameth_copy"), "Copy"),
                    br(),
                    br(),
                    dataTableOutput(ns("display_info_amameth"))
             ),
             column(5,
                    fluidRow(
                      column(4,
                             h4(textOutput(ns("title_at"))),
                             rHandsontableOutput(ns("anameth_edit_at"))
                      ),
                      column(8,
                             h4(textOutput(ns("title_al_gt_prob"))),
                             rHandsontableOutput(ns("anameth_edit_algt_prob")),
                             br(),
                             h4(textOutput(ns("title_gt_fltr"))),
                             rHandsontableOutput(ns("anameth_edit_gt_fltr")),
                             br(),
                             h4(textOutput(ns("title_dropin"))),
                             rHandsontableOutput(ns("anameth_edit_dropin")),
                             br(),
                             h4(textOutput(ns("title_mcmc"))),
                             rHandsontableOutput(ns("anameth_edit_mcmc")),
                             br(),
                             uiOutput(ns("act_anameth_save"))
                      )
                    )
             )
           )
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
    info_param <- info_stutter <- info_kit <- rep("", n_anameth)
    for(i in 1:n_anameth){
      dt_anameth_1 <- fread(paste0(path_pack, "/extdata/anameth/", fn_anameth[i]))
      info_param[i] <- dt_anameth_1$Name[dt_anameth_1$Category == "Parameter"]
      info_stutter[i] <- dt_anameth_1$Name[dt_anameth_1$Category == "Stutter file"]
      info_kit[i] <- dt_anameth_1$Name[dt_anameth_1$Category == "Kit"]
    }

    dt_info_anameth <- data.table("Analysis method" = gsub(".{4}$", "", fn_anameth),
                                  "Parameter data" = info_param,
                                  "Stutter file" = info_stutter,
                                  "Kit" = info_kit)
    return(dt_info_anameth)
  }else{
    return(NULL)
  }
}

#' create_default_dt_anameth
#'
#' @description The function to create the default data.table for a new analysis method
#' @param path_pack The package path of KongohPlus
#' @param parameter_name A parameter name
#' @param stutter_file_name A stutter file name
#' @param kit_name A kit name
create_default_dt_anameth <- function(path_pack, parameter_name, stutter_file_name, kit_name){
  dt_kit <- fread(paste0(path_pack, "/extdata/kit/", kit_name, ".csv"))
  markers_auto <- unique(dt_kit[["Marker"]][dt_kit[["Chr"]] == "Auto"])
  n_mk <- length(markers_auto)

  default_dt_anameth <- data.table(Category = c("Parameter", "Stutter file", "Kit",
                                                rep("AT", n_mk),
                                                rep("Genotype probability", 1),
                                                rep("Genotype filter", 6),
                                                rep("Drop-in", 3),
                                                rep("MCMC", 3)),
                                   Subcategory = c("", "", "",
                                                   markers_auto,
                                                   "Theta",
                                                   "Maximum BSR", "Maximum FSR", "Maximum DSR", "Maximum M2SR", "Maximum number of drop-in", "Maximum height of drop-in",
                                                   "Probability of drop-in", "Shape parameter", "Rate parameter",
                                                   "Warmup", "Sampling", "Chains"),
                                   Name = c(parameter_name, stutter_file_name, kit_name, rep("", n_mk + 14)),
                                   Value = c(rep(NA, 3),
                                             rep(100, n_mk),
                                             0,
                                             0.3, 0.15, 0.05, 0.05, 2, 300,
                                             0.0037081, 14.8712, 0.2029315,
                                             1200, 300, 4))
  return(default_dt_anameth)
}

#' anameth_server
#'
#' @description The function to create the server module for the analysis method
anameth_server <- function(id, path_pack){
  moduleServer(
    id,
    function(input, output, session){
      # Define the initial reactive values
      rv_anameth <- reactiveValues()
      rv_anameth$dt_info_anameth <- create_dt_info_anameth(path_pack)

      # Define a new analysis method
      observeEvent(input$act_anameth_new, {
        param_all <- list.files(paste0(path_pack, "/extdata/param"))
        param_all <- gsub(".{6}$", "", param_all)

        showModal(modalDialog(
          title = "Define a new analysis method",
          textInput(session$ns("anameth_name"), label = "Enter a name of the analysis method"),
          br(),
          selectInput(session$ns("parameter_name"), label = h4("Select parameter data"), choices = param_all, selected = param_all[1], selectize = FALSE),
          footer = tagList(
            actionButton(session$ns("act_add_anameth_new"), "Add the analysis method"),
            modalButton("Close")
          )
        ))
      })

      observeEvent(input$act_add_anameth_new, {
        anameth_name <- input$anameth_name
        parameter_name <- input$parameter_name

        if(!is.null(parameter_name)){
          load(paste0(path_pack, "/extdata/param/", parameter_name, ".RData"))
          stutter_file_name <- gsub(".{6}$", "", fn_stutter_file)

          dt_info_anameth <- rv_anameth$dt_info_anameth
          anameth_all <- c(dt_info_anameth[["Analysis method"]], anameth_name)
          parameter_all <- c(dt_info_anameth[["Parameter data"]], parameter_name)
          stutter_file_all <- c(dt_info_anameth[["Stutter file"]], stutter_file_name)
          kit_all <- c(dt_info_anameth[["Kit"]], kit_name)

          dt_info_anameth <- data.table("Analysis method" = anameth_all,
                                        "Parameter data" = parameter_all,
                                        "Stutter file" = stutter_file_all,
                                        "Kit" = kit_all)
          rv_anameth$dt_info_anameth <- dt_info_anameth

          default_dt_anameth <- create_default_dt_anameth(path_pack, parameter_name, stutter_file_name, kit_name)
          write.csv(default_dt_anameth, paste0(path_pack, "/extdata/anameth/", anameth_name, ".csv"), row.names = FALSE)

          output$title_at <- renderText({NULL})
          output$title_al_gt_prob <- renderText({NULL})
          output$title_gt_fltr <- renderText({NULL})
          output$title_dropin <- renderText({NULL})
          output$title_mcmc <- renderText({NULL})
          output$anameth_edit_at <- renderRHandsontable({NULL})
          output$anameth_edit_algt_prob <- renderRHandsontable({NULL})
          output$anameth_edit_gt_fltr <- renderRHandsontable({NULL})
          output$anameth_edit_dropin <- renderRHandsontable({NULL})
          output$anameth_edit_mcmc <- renderRHandsontable({NULL})
          output$act_anameth_save <- renderUI({NULL})
          removeModal()
        }
      })

      # Delete an analysis method
      observeEvent(input$act_anameth_delete, {
        index_select <- input$display_info_amameth_rows_selected
        if(isTruthy(index_select)){
          dt_info_anameth <- rv_anameth$dt_info_anameth
          file.remove(paste0(path_pack, "/extdata/anameth/", dt_info_anameth[["Analysis method"]][index_select], ".csv"))
          rv_anameth$dt_info_anameth <- dt_info_anameth[-index_select, ]

          output$title_at <- renderText({NULL})
          output$title_al_gt_prob <- renderText({NULL})
          output$title_gt_fltr <- renderText({NULL})
          output$title_dropin <- renderText({NULL})
          output$title_mcmc <- renderText({NULL})
          output$anameth_edit_at <- renderRHandsontable({NULL})
          output$anameth_edit_algt_prob <- renderRHandsontable({NULL})
          output$anameth_edit_gt_fltr <- renderRHandsontable({NULL})
          output$anameth_edit_dropin <- renderRHandsontable({NULL})
          output$anameth_edit_mcmc <- renderRHandsontable({NULL})
          output$act_anameth_save <- renderUI({NULL})
        }else{
          showModal(modalDialog(title = "Error", "Select an analysis method!", easyClose = TRUE, footer = NULL))
        }
      })

      # Copy an analysis method
      observeEvent(input$act_anameth_copy, {
        index_select <- input$display_info_amameth_rows_selected
        if(isTruthy(index_select)){
          # Update information on the analysis method
          dt_info_anameth <- rv_anameth$dt_info_anameth
          dt_info_anameth_copy <- copy(dt_info_anameth)
          dt_info_anameth_copy <- dt_info_anameth_copy[index_select, ]
          anameth_select <- dt_info_anameth_copy[["Analysis method"]]
          dt_info_anameth_copy[["Analysis method"]] <- paste0(anameth_select, "_copy")
          dt_info_anameth <- rbind(dt_info_anameth, dt_info_anameth_copy)
          rv_anameth$dt_info_anameth <- dt_info_anameth

          # Make a copy of the selected analysis method
          dt_anameth_copy <- fread(paste0(path_pack, "/extdata/anameth/", anameth_select, ".csv"))
          write.csv(dt_anameth_copy, paste0(path_pack, "/extdata/anameth/", dt_info_anameth_copy[["Analysis method"]], ".csv"), row.names = FALSE)

          output$title_at <- renderText({NULL})
          output$title_al_gt_prob <- renderText({NULL})
          output$title_gt_fltr <- renderText({NULL})
          output$title_dropin <- renderText({NULL})
          output$title_mcmc <- renderText({NULL})
          output$anameth_edit_at <- renderRHandsontable({NULL})
          output$anameth_edit_algt_prob <- renderRHandsontable({NULL})
          output$anameth_edit_gt_fltr <- renderRHandsontable({NULL})
          output$anameth_edit_dropin <- renderRHandsontable({NULL})
          output$anameth_edit_mcmc <- renderRHandsontable({NULL})
          output$act_anameth_save <- renderUI({NULL})
        }else{
          showModal(modalDialog(title = "Error", "Select an analysis method!", easyClose = TRUE, footer = NULL))
        }
      })

      # Display information on the analysis method
      output$display_info_amameth <- renderDataTable({
        dt_info_anameth <- rv_anameth$dt_info_anameth
        datatable(
          dt_info_anameth,
          selection = list(mode = "single", target = "row"),
          #extensions = c("FixedColumns", "FixedHeader"),
          options = list(dom = "t", iDisplayLength = nrow(dt_info_anameth), ordering = FALSE, autoWidth = TRUE, scrollY = "650px", #scrollX = TRUE, scrollY = "650px",
                         fixedHeader = TRUE), #fixedColumns = list(leftColumns = 1, rightColumns = 0),
                         #columnDefs = list(list(targets = c(0, 2), width = "300px"),
                         #                  list(targets = c(1, 3), width = "200px"))),
          rownames = FALSE
        )
      })

      observeEvent(input$display_info_amameth_rows_selected, {
        index_select <- input$display_info_amameth_rows_selected
        if(isTruthy(index_select)){
          dt_info_anameth <- rv_anameth$dt_info_anameth
          anameth_select <- dt_info_anameth[["Analysis method"]][index_select]
          dt_anameth_edit <- fread(paste0(path_pack, "/extdata/anameth/", anameth_select, ".csv"))

          # Set the parameter name
          parameter_name <- dt_anameth_edit$Name[dt_anameth_edit$Category == "Parameter"]

          # Set the stutter file
          fn_stutter_file <- dt_anameth_edit$Name[dt_anameth_edit$Category == "Stutter file"]

          # Set the kit name
          kit_name <- dt_anameth_edit$Name[dt_anameth_edit$Category == "Kit"]

          output$title_at <- renderText({"Analytical threshold"})
          output$title_al_gt_prob <- renderText({"Genotype probability"})
          output$title_gt_fltr <- renderText({"Genotype filter"})
          output$title_dropin <- renderText({"Drop-in"})
          output$title_mcmc <- renderText({"MCMC"})

          output$anameth_edit_at <- renderRHandsontable({
            dt_at <- data.table(Marker = dt_anameth_edit$Subcategory[dt_anameth_edit$Category == "AT"],
                                AT = dt_anameth_edit$Value[dt_anameth_edit$Category == "AT"])
            rhandsontable(dt_at, rowHeaders = FALSE) %>%
              hot_table(fillHandle = FALSE, rowHeights = 30) %>%
              hot_col("Marker", type = "text", readOnly = TRUE, width = 100) %>%
              hot_col("AT", type = "numeric", format = "0", width = 100)
          })

          output$anameth_edit_algt_prob <- renderRHandsontable({
            dt_algt_prob <- data.table(Subcategory = dt_anameth_edit$Subcategory[dt_anameth_edit$Category == "Genotype probability"],
                                       Value = dt_anameth_edit$Value[dt_anameth_edit$Category == "Genotype probability"])
            rhandsontable(dt_algt_prob, rowHeaders = FALSE, colHeaders = FALSE) %>%
              hot_table(fillHandle = FALSE, rowHeights = 30) %>%
              hot_col(col = 1, type = "text", readOnly = TRUE, width = 300) %>%
              hot_col(col = 2, type = "numeric", width = 100)
          })

          output$anameth_edit_gt_fltr <- renderRHandsontable({
            dt_gt_fltr <- data.table(Subcategory = dt_anameth_edit$Subcategory[dt_anameth_edit$Category == "Genotype filter"],
                                     Value = dt_anameth_edit$Value[dt_anameth_edit$Category == "Genotype filter"])
            rhandsontable(dt_gt_fltr , rowHeaders = FALSE, colHeaders = FALSE) %>%
              hot_table(fillHandle = FALSE, rowHeights = 30) %>%
              hot_col(col = 1, type = "text", readOnly = TRUE, width = 300) %>%
              hot_col(col = 2, type = "numeric", width = 100)
          })

          output$anameth_edit_dropin <- renderRHandsontable({
            dt_dropin <- data.table(Subcategory = dt_anameth_edit$Subcategory[dt_anameth_edit$Category == "Drop-in"],
                                    Value = dt_anameth_edit$Value[dt_anameth_edit$Category == "Drop-in"])
            rhandsontable(dt_dropin, rowHeaders = FALSE, colHeaders = FALSE) %>%
              hot_table(fillHandle = FALSE, rowHeights = 30) %>%
              hot_col(col = 1, type = "text", readOnly = TRUE, width = 300) %>%
              hot_col(col = 2, type = "numeric", width = 100)
          })

          output$anameth_edit_mcmc <- renderRHandsontable({
            dt_mcmc <- data.table(Subcategory = dt_anameth_edit$Subcategory[dt_anameth_edit$Category == "MCMC"],
                                  Value = dt_anameth_edit$Value[dt_anameth_edit$Category == "MCMC"])
            rhandsontable(dt_mcmc, rowHeaders = FALSE, colHeaders = FALSE) %>%
              hot_table(fillHandle = FALSE, rowHeights = 30) %>%
              hot_col(col = 1, type = "text", readOnly = TRUE, width = 300) %>%
              hot_col(col = 2, type = "numeric", width = 100)
          })

          output$act_anameth_save <- renderUI({
            actionButton(session$ns("act_anameth_save"), "Save")
          })
        }
      })

      # Save the analysis method
      observeEvent(input$act_anameth_save, {
        # Extract settings
        index_select <- input$display_info_amameth_rows_selected
        dt_info_anameth <- rv_anameth$dt_info_anameth
        anameth_select <- dt_info_anameth[["Analysis method"]][index_select]
        parameter_name <- dt_info_anameth[["Parameter data"]][index_select]
        fn_stutter_file <- dt_info_anameth[["Stutter file"]][index_select]
        kit_name <- dt_info_anameth[["Kit"]][index_select]

        dt_at <- hot_to_r(input$anameth_edit_at)
        dt_algt_prob <- hot_to_r(input$anameth_edit_algt_prob)
        dt_gt_fltr <- hot_to_r(input$anameth_edit_gt_fltr)
        dt_dropin <- hot_to_r(input$anameth_edit_dropin)
        dt_mcmc <- hot_to_r(input$anameth_edit_mcmc)

        # Create a data.table for the new analysis method
        dt_anameth <- data.table(Category = c("Parameter", "Stutter file", "Kit",
                                              rep("AT", nrow(dt_at)),
                                              rep("Genotype probability", nrow(dt_algt_prob)),
                                              rep("Genotype filter", nrow(dt_gt_fltr)),
                                              rep("Drop-in", nrow(dt_dropin)),
                                              rep("MCMC", nrow(dt_mcmc))),
                                 Subcategory = c("", "", "",
                                                 dt_at[["Marker"]],
                                                 dt_algt_prob[["Subcategory"]],
                                                 dt_gt_fltr[["Subcategory"]],
                                                 dt_dropin[["Subcategory"]],
                                                 dt_mcmc[["Subcategory"]]),
                                 Name = c(parameter_name, fn_stutter_file, kit_name, rep("", nrow(dt_at) + nrow(dt_algt_prob) + nrow(dt_gt_fltr) + nrow(dt_dropin) + nrow(dt_mcmc))),
                                 Value = c(rep(NA, 3),
                                           dt_at[["AT"]],
                                           dt_algt_prob[["Value"]],
                                           dt_gt_fltr[["Value"]],
                                           dt_dropin[["Value"]],
                                           dt_mcmc[["Value"]]))
        file.remove(paste0(path_pack, "/extdata/anameth/", anameth_select, ".csv"))
        write.csv(dt_anameth, paste0(path_pack, "/extdata/anameth/", anameth_select, ".csv"), row.names = FALSE)

        # Update information on the analysis methods
        rv_anameth$dt_info_anameth <- create_dt_info_anameth(path_pack)

        showModal(modalDialog(title = "Message", "Save completed successfully.", easyClose = TRUE, footer = NULL))
      }, ignoreInit = TRUE)
    }
  )
}
