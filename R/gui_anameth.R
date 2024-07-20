#' anameth_ui
#'
#' @description The function to create the ui module for the analysis method
anameth_ui <- function(id){
  ns <- NS(id)

  tabPanel("Analysis method",
           titlePanel("Analysis method"),
           br(),
           sidebarLayout(
             sidebarPanel(width = 2,
                          actionButton(ns("act_anameth_new"), "New"),
                          br(),
                          br(),
                          actionButton(ns("act_anameth_edit"), "Edit"),
                          br(),
                          br(),
                          actionButton(ns("act_anameth_delete"), "Delete"),
                          br(),
                          br(),
                          actionButton(ns("act_anameth_copy"), "Copy")
             ),
             mainPanel(width = 10,
                       tags$head(
                         tags$style(HTML("table {table-layout: fixed;")),
                         tags$style(".datatables .display {margin-left: 0;}")
                       ),
                       dataTableOutput(ns("display_info_amameth"))
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
    info_kit <- info_param <- info_at <- info_allele_prob <- rep("", n_anameth)
    for(i in 1:n_anameth){
      dt_anameth_1 <- fread(paste0(path_pack, "/extdata/anameth/", fn_anameth[i]))

      # Set information on the kit
      info_kit[i] <- dt_anameth_1$Value[dt_anameth_1$Category == "Kit"]

      # Set information on the parameters
      info_param[i] <- dt_anameth_1$Value[dt_anameth_1$Category == "Parameter"]

      # Set information on the AT
      ats <- as.numeric(dt_anameth_1$Value[dt_anameth_1$Category == "AT"])
      at_unique <- unique(ats)
      if(length(at_unique) == 1){
        info_at[i] <- as.character(at_unique)
      }else{
        info_at[i] <- paste0(min(at_unique), "-", max(at_unique))
      }

      # Set information on the allele probability
      info_allele_prob[i] <- dt_anameth_1$Value[dt_anameth_1$Subcategory == "Method"]
    }

    dt_info_anameth <- data.table(Name = gsub(".{4}$", "", fn_anameth),
                                  Kit = info_kit,
                                  Parameter = info_param,
                                  AT = info_at,
                                  Allele_Probability = info_allele_prob)
    return(dt_info_anameth)
  }else{
    return(NULL)
  }
}

#' anameth_server
#'
#' @description The function to create the server module for the analysis method
anameth_server <- function(id, path_pack){
  moduleServer(
    id,
    function(input, output, session){

      numericColumn <- function(len, col, val, ...){
        inputs <- character(len)
        for(i in seq_len(len)){
          inputs[i] <- as.character(
            numericInput(paste0("numeric_", col, "_", i), label = NULL, value = val[i], ...)
          )
        }
        return(inputs)
      }

      js <- function(dtid, col_numeric){
        code <- vector("list", 1)
        code[[1]] <- c(
          sprintf("$('body').on('change', '[id^=numeric_%d_]', function() {", col_numeric),
          "  var id = this.getAttribute('id');",
          sprintf("  var i = parseInt(/numeric_%d_(\\d+)/.exec(id)[1]);", col_numeric),
          "  var value = $(this).val();",
          sprintf("  var info = [{row: i, col: %d, value: value}];", col_numeric),
          sprintf("  Shiny.setInputValue('%s', info);", session$ns(sprintf("%s_cell_edit:DT.cellInfo", dtid))),
          "});"
        )
        do.call(c, code)
      }

      ######################################
      # Define the initial reactive values #
      ######################################

      rv_anameth <- reactiveValues()
      rv_anameth$dt_info_anameth <- create_dt_info_anameth(path_pack)
      rv_anameth$dt_anameth <- NULL

      rv_anameth$anameth_new_kit <- NULL
      rv_anameth$dt_at_new <- NULL
      rv_anameth$dt_at_edit <- NULL

      ##############################
      # Make a new analysis method #
      ##############################

      observeEvent(input$act_anameth_new, {
        showModal(modalDialog(
          title = "Create a new analysis method",
          textInput(session$ns("anameth_new_input_name"), label = "Enter a name of the new analysis method", value = NULL),
          tabsetPanel(id = session$ns("tabset_new_anameth"),
                      type = "tabs",
                      tabPanel("Parameter & Kit",
                               br(),
                               uiOutput(session$ns("anameth_new_parameter")),
                               br(),
                               h4("Kit"),
                               textOutput(session$ns("anameth_new_kit"))
                      ),
                      tabPanel("AT",
                               br(),
                               dataTableOutput(session$ns("anameth_new_at"))
                      ),
                      tabPanel("Allele probability",
                               br(),
                               selectInput(session$ns("anameth_new_alprob"), label = h4("Method"), choices = c("Dirichlet", "Observed"), selectize = FALSE),
                               br(),
                               numericInput(session$ns("anameth_new_maf"), label = h4("Minimum allele frequency"), value = 0),
                               br(),
                               numericInput(session$ns("anameth_new_theta"), label = h4("Theta correction"), value = 0)
                      ),
                      tabPanel("Genotype filter",
                               br(),
                               numericInput(session$ns("anameth_new_gt_fltr_bsr"), label = h4("Maximum back stutter ratio"), value = 0.3),
                               br(),
                               numericInput(session$ns("anameth_new_gt_fltr_fsr"), label = h4("Maximum forward stutter ratio"), value = 0.15),
                               br(),
                               numericInput(session$ns("anameth_new_gt_fltr_dsr"), label = h4("Maximum double-back stutter ratio"), value = 0.05),
                               br(),
                               numericInput(session$ns("anameth_new_gt_fltr_m2sr"), label = h4("Maximum minus 2-nt stutter ratio"), value = 0.05),
                               br(),
                               numericInput(session$ns("anameth_new_gt_fltr_n_di"), label = h4("Maximum number of drop-in peaks per locus"), value = 2),
                               br(),
                               numericInput(session$ns("anameth_new_gt_fltr_h_di"), label = h4("Maximum height of drop-in peaks"), value = 300)
                      ),
                      tabPanel("Drop-in parameter",
                               br(),
                               numericInput(session$ns("anameth_new_di_prob"), label = h4("Probability of drop-in"), value = 0.00335),
                               br(),
                               numericInput(session$ns("anameth_new_di_shape"), label = h4("Shape parameter"), value = 1.22),
                               br(),
                               numericInput(session$ns("anameth_new_di_rate"), label = h4("Rate parameter"), value = 0.0503),
                               br(),
                               numericInput(session$ns("anameth_new_di_mt"), label = h4("Minimum threshold"), value = 50)
                      )
          ),
          footer = tagList(
            actionButton(session$ns("act_new_anameth_save"), "Save"),
            modalButton("Cancel")
          ),
          size = "l"
        ))

        output$anameth_new_parameter <- renderUI({
          param_all <- list.files(paste0(path_pack, "/extdata/param"))
          param_all <- gsub(".{6}$", "", param_all)
          selectInput(session$ns("anameth_new_parameter"), label = h4("Parameter data"), choices = param_all, selected = param_all[1], selectize = FALSE)
        })

        observeEvent(input$anameth_new_parameter, {
          anameth_new_parameter <- input$anameth_new_parameter
          if(is.null(anameth_new_parameter)){
            output$anameth_new_kit <- renderText({"NA"})
          }else{
            load(paste0(path_pack, "/extdata/param/", anameth_new_parameter, ".Rdata"))

            # Set the kit name
            rv_anameth$anameth_new_kit <- anameth_new_kit <- dt_par$Kit[1]
            output$anameth_new_kit <- renderText({paste0(anameth_new_kit)})

            # Extract autosomal markers
            dt_kit_select <- fread(paste0(path_pack, "/extdata/kit/", anameth_new_kit, ".csv"))
            dt_kit_select <- change_data_type_kit(dt_kit_select)
            markers <- unique(dt_kit_select$Marker[dt_kit_select$Chr == "Auto"])
            n_mk <- length(markers)

            # Create a data.table for the AT
            dt_at_new <- data.table(Marker = markers) # For saved data
            dt_at_new_display <- copy(dt_at_new) # For displayed data
            dt_at_new$AT <- 100
            rv_anameth$dt_at_new <- dt_at_new
            dt_at_new_display$AT <- numericColumn(n_mk, 2, val = dt_at_new$AT, width = "75%")

            # Display a data.table for the AT
            output$anameth_new_at <- renderDataTable(server = FALSE, {
              datatable(
                dt_at_new_display,
                selection = "none",
                options = list(dom = "t", iDisplayLength = nrow(dt_at_new_display), ordering = FALSE, autoWidth = TRUE),
                rownames = FALSE,
                escape = FALSE,
                callback = JS(js("anameth_new_at", 2))
              )
            })
          }
        })

        # Change the state of the numericInput for the MAF
        observeEvent(input$anameth_new_alprob, {
          if(input$anameth_new_alprob == "Dirichlet distribution"){
            disable("anameth_new_maf")
          }else{
            enable("anameth_new_maf")
          }
        })

        # Update AT
        observeEvent(input[["anameth_new_at_cell_edit"]], {
          info <- input[["anameth_new_at_cell_edit"]]
          dt_at_new <- rv_anameth$dt_at_new
          dt_at_new[[info[1, "col"]]][info[1, "row"]] <- info[1, "value"]
          rv_anameth$dt_at_new <- dt_at_new
        })

        # Save the new analysis method
        observeEvent(input$act_new_anameth_save, {
          # Extract settings
          anameth_new_input_name <- input$anameth_new_input_name
          anameth_new_parameter <- input$anameth_new_parameter
          anameth_new_kit <- rv_anameth$anameth_new_kit
          dt_at_new <- rv_anameth$dt_at_new
          markers <- dt_at_new$Marker
          ats <- dt_at_new$AT
          anameth_new_alprob <- input$anameth_new_alprob
          anameth_new_maf <- input$anameth_new_maf
          anameth_new_theta <- input$anameth_new_theta
          anameth_new_gt_fltr_bsr <- input$anameth_new_gt_fltr_bsr
          anameth_new_gt_fltr_fsr <- input$anameth_new_gt_fltr_fsr
          anameth_new_gt_fltr_dsr <- input$anameth_new_gt_fltr_dsr
          anameth_new_gt_fltr_m2sr <- input$anameth_new_gt_fltr_m2sr
          anameth_new_gt_fltr_n_di <- input$anameth_new_gt_fltr_n_di
          anameth_new_gt_fltr_h_di <- input$anameth_new_gt_fltr_h_di
          anameth_new_di_prob <- input$anameth_new_di_prob
          anameth_new_di_shape <- input$anameth_new_di_shape
          anameth_new_di_rate <- input$anameth_new_di_rate
          anameth_new_di_mt <- input$anameth_new_di_mt

          # Create a data.table for the new analysis method
          dt_anameth_new <- data.table(Category = c("Parameter", "Kit", rep("AT", length(markers)), rep("Allele_probability", 3), rep("Genotype_filter", 6), rep("Drop-in", 4)),
                                       Subcategory = c("", "", markers, "Method", "Maf", "Theta", "BSR", "FSR", "DSR", "M2SR", "N_Dropin", "H_Dropin", "Probability", "Shape", "Rate", "MT"),
                                       Value = c(anameth_new_parameter, anameth_new_kit, ats,
                                                 anameth_new_alprob, anameth_new_maf, anameth_new_theta,
                                                 anameth_new_gt_fltr_bsr, anameth_new_gt_fltr_fsr, anameth_new_gt_fltr_dsr, anameth_new_gt_fltr_m2sr, anameth_new_gt_fltr_n_di, anameth_new_gt_fltr_h_di,
                                                 anameth_new_di_prob, anameth_new_di_shape, anameth_new_di_rate, anameth_new_di_mt))
          write.csv(dt_anameth_new, paste0(path_pack, "/extdata/anameth/", anameth_new_input_name, ".csv"), row.names = FALSE)

          # Update information on the analysis methods
          rv_anameth$dt_info_anameth <- create_dt_info_anameth(path_pack)

          removeModal()
        }, ignoreInit = TRUE)

      })

      ###########################
      # Edit an analysis method #
      ###########################

      observeEvent(input$act_anameth_edit, {
        pos_select <- input$display_info_amameth_rows_selected
        if(isTruthy(pos_select)){
          # Load the selected analysis method
          dt_info_anameth <- rv_anameth$dt_info_anameth
          anameth_select <- dt_info_anameth$Name[pos_select]
          rv_anameth$dt_anameth <- dt_anameth_select <- fread(paste0(path_pack, "/extdata/anameth/", anameth_select, ".csv"))

          # Set the parameter name
          parameter_name <- dt_anameth_select$Value[dt_anameth_select$Category == "Parameter"]

          # Set the kit name
          kit_name <- dt_anameth_select$Value[dt_anameth_select$Category == "Kit"]

          # Extract values for the allele probability
          alprob_method <- dt_anameth_select$Value[dt_anameth_select$Category == "Allele_probability" & dt_anameth_select$Subcategory == "Method"]
          alprob_maf <- as.numeric(dt_anameth_select$Value[dt_anameth_select$Category == "Allele_probability" & dt_anameth_select$Subcategory == "Maf"])
          alprob_theta <- as.numeric(dt_anameth_select$Value[dt_anameth_select$Category == "Allele_probability" & dt_anameth_select$Subcategory == "Theta"])

          # Extract values for the genotype filters
          gtfltr_bsr <- dt_anameth_select$Value[dt_anameth_select$Category == "Genotype_filter" & dt_anameth_select$Subcategory == "BSR"]
          gtfltr_fsr <- dt_anameth_select$Value[dt_anameth_select$Category == "Genotype_filter" & dt_anameth_select$Subcategory == "FSR"]
          gtfltr_dsr <- dt_anameth_select$Value[dt_anameth_select$Category == "Genotype_filter" & dt_anameth_select$Subcategory == "DSR"]
          gtfltr_m2sr <- dt_anameth_select$Value[dt_anameth_select$Category == "Genotype_filter" & dt_anameth_select$Subcategory == "M2SR"]
          gtfltr_n_di <- dt_anameth_select$Value[dt_anameth_select$Category == "Genotype_filter" & dt_anameth_select$Subcategory == "N_Dropin"]
          gtfltr_h_di <- dt_anameth_select$Value[dt_anameth_select$Category == "Genotype_filter" & dt_anameth_select$Subcategory == "H_Dropin"]

          # Extract values for the drop-in parameters
          di_prob <- dt_anameth_select$Value[dt_anameth_select$Category == "Drop-in" & dt_anameth_select$Subcategory == "Probability"]
          di_shape <- dt_anameth_select$Value[dt_anameth_select$Category == "Drop-in" & dt_anameth_select$Subcategory == "Shape"]
          di_rate <- dt_anameth_select$Value[dt_anameth_select$Category == "Drop-in" & dt_anameth_select$Subcategory == "Rate"]
          di_mt <- dt_anameth_select$Value[dt_anameth_select$Category == "Drop-in" & dt_anameth_select$Subcategory == "MT"]

          showModal(modalDialog(
            title = "Edit an analysis method",
            textInput(session$ns("anameth_edit_input_name"), label = "The name of the analysis method", value = anameth_select),
            tabsetPanel(id = session$ns("tabset_edit_anameth"),
                        type = "tabs",
                        tabPanel("Parameter & Kit",
                                 br(),
                                 h4("Parameter data"),
                                 p(paste0(parameter_name, ".Rdata")),
                                 br(),
                                 h4("Kit"),
                                 p(kit_name)
                        ),
                        tabPanel("AT",
                                 br(),
                                 dataTableOutput(session$ns("anameth_edit_at"))
                        ),
                        tabPanel("Allele probability",
                                 br(),
                                 selectInput(session$ns("anameth_edit_alprob"), label = h4("Method"), choices = c("Dirichlet", "Observed"), selected = alprob_method, selectize = FALSE),
                                 br(),
                                 numericInput(session$ns("anameth_edit_maf"), label = h4("Minimum allele frequency"), value = alprob_maf),
                                 br(),
                                 numericInput(session$ns("anameth_edit_theta"), label = h4("Theta correction"), value = alprob_theta)
                        ),
                        tabPanel("Genotype filter",
                                 br(),
                                 numericInput(session$ns("anameth_edit_gt_fltr_bsr"), label = h4("Maximum back stutter ratio"), value = gtfltr_bsr),
                                 br(),
                                 numericInput(session$ns("anameth_edit_gt_fltr_fsr"), label = h4("Maximum forward stutter ratio"), value = gtfltr_fsr),
                                 br(),
                                 numericInput(session$ns("anameth_edit_gt_fltr_dsr"), label = h4("Maximum double-back stutter ratio"), value = gtfltr_dsr),
                                 br(),
                                 numericInput(session$ns("anameth_edit_gt_fltr_m2sr"), label = h4("Maximum minus 2-nt stutter ratio"), value = gtfltr_m2sr),
                                 br(),
                                 numericInput(session$ns("anameth_edit_gt_fltr_n_di"), label = h4("Maximum number of drop-in peaks per locus"), value = gtfltr_n_di),
                                 br(),
                                 numericInput(session$ns("anameth_edit_gt_fltr_h_di"), label = h4("Maximum height of drop-in peaks"), value = gtfltr_h_di)
                        ),
                        tabPanel("Drop-in parameter",
                                 br(),
                                 numericInput(session$ns("anameth_edit_di_prob"), label = h4("Probability of drop-in"), value = di_prob),
                                 br(),
                                 numericInput(session$ns("anameth_edit_di_shape"), label = h4("Shape parameter"), value = di_shape),
                                 br(),
                                 numericInput(session$ns("anameth_edit_di_rate"), label = h4("Rate parameter"), value = di_rate),
                                 br(),
                                 numericInput(session$ns("anameth_edit_di_mt"), label = h4("Minimum threshold"), value = di_mt)
                        )
            ),
            footer = tagList(
              actionButton(session$ns("act_edit_anameth_save"), "Save"),
              modalButton("Cancel")
            ),
            size = "l"
          ))
        }else{
          showModal(modalDialog(title = "Error", "Select an analysis method!", easyClose = TRUE, footer = NULL))
        }
      })

      output$anameth_edit_at <- renderDataTable(server = FALSE, {
        # Extract autosomal markers
        dt_anameth_select <- rv_anameth$dt_anameth
        markers <- dt_anameth_select$Subcategory[dt_anameth_select$Category == "AT"]
        n_mk <- length(markers)

        # Create a data.table for the AT
        dt_at_edit <- data.table(Marker = markers) # For saved data
        dt_at_edit_display <- copy(dt_at_edit) # For displayed data
        dt_at_edit$AT <- as.numeric(dt_anameth_select$Value[dt_anameth_select$Category == "AT"])
        rv_anameth$dt_at_edit <- dt_at_edit
        dt_at_edit_display$AT <- numericColumn(n_mk, 2, val = dt_at_edit$AT, width = "75%")

        datatable(
          dt_at_edit_display,
          selection = "none",
          options = list(dom = "t", iDisplayLength = nrow(dt_at_edit_display), ordering = FALSE, autoWidth = TRUE),
          rownames = FALSE,
          escape = FALSE,
          callback = JS(js("anameth_edit_at", 2))
        )
      })

      # Update AT
      observeEvent(input[["anameth_edit_at_cell_edit"]], {
        info <- input[["anameth_edit_at_cell_edit"]]
        dt_at_edit <- rv_anameth$dt_at_edit
        dt_at_edit[[info[1, "col"]]][info[1, "row"]] <- info[1, "value"]
        rv_anameth$dt_at_edit <- dt_at_edit
      })

      # Save the analysis method
      observeEvent(input$act_edit_anameth_save, {
        # Extract settings
        pos_select <- input$display_info_amameth_rows_selected
        dt_info_anameth <- rv_anameth$dt_info_anameth
        anameth_select <- dt_info_anameth$Name[pos_select]

        anameth_edit_input_name <- input$anameth_edit_input_name
        dt_anameth_select <- rv_anameth$dt_anameth
        parameter_name <- dt_anameth_select$Value[dt_anameth_select$Category == "Parameter"]
        kit_name <- dt_anameth_select$Value[dt_anameth_select$Category == "Kit"]

        dt_at_edit <- rv_anameth$dt_at_edit
        markers <- dt_at_edit$Marker
        ats <- dt_at_edit$AT

        # Create a data.table for the new analysis method
        dt_anameth_new <- data.table(Category = c("Parameter", "Kit", rep("AT", length(markers)), rep("Allele_probability", 3), rep("Genotype_filter", 6), rep("Drop-in", 4)),
                                     Subcategory = c("", "", markers, "Method", "Maf", "Theta", "BSR", "FSR", "DSR", "M2SR", "N_Dropin", "H_Dropin", "Probability", "Shape", "Rate", "MT"),
                                     Value = c(parameter_name, kit_name, ats,
                                               input$anameth_edit_alprob, input$anameth_edit_maf, input$anameth_edit_theta,
                                               input$anameth_edit_gt_fltr_bsr, input$anameth_edit_gt_fltr_fsr, input$anameth_edit_gt_fltr_dsr, input$anameth_edit_gt_fltr_m2sr,
                                               input$anameth_edit_gt_fltr_n_di, input$anameth_edit_gt_fltr_h_di,
                                               input$anameth_edit_di_prob, input$anameth_edit_di_shape, input$anameth_edit_di_rate, input$anameth_edit_di_mt))
        file.remove(paste0(path_pack, "/extdata/anameth/", anameth_select, ".csv"))
        write.csv(dt_anameth_new, paste0(path_pack, "/extdata/anameth/", anameth_edit_input_name, ".csv"), row.names = FALSE)

        # Update information on the analysis methods
        rv_anameth$dt_info_anameth <- create_dt_info_anameth(path_pack)

        removeModal()
      }, ignoreInit = TRUE)

      #############################
      # Delete an analysis method #
      #############################

      observeEvent(input$act_anameth_delete, {
        pos_select <- input$display_info_amameth_rows_selected
        if(isTruthy(pos_select)){
          showModal(modalDialog(
            title = "Delete an analysis method",
            "Click Ok button to delete the selected analysis method.",
            footer = tagList(
              actionButton(session$ns("act_delete_anameth_ok"), "Ok"),
              modalButton("Cancel")
            ),
            size = "l"
          ))
        }else{
          showModal(modalDialog(title = "Error", "Select an analysis method!", easyClose = TRUE, footer = NULL))
        }
      })

      observeEvent(input$act_delete_anameth_ok, {
        pos_select <- input$display_info_amameth_rows_selected
        dt_info_anameth <- rv_anameth$dt_info_anameth
        file.remove(paste0(path_pack, "/extdata/anameth/", dt_info_anameth$Name[pos_select], ".csv"))
        rv_anameth$dt_info_anameth <- dt_info_anameth[-pos_select, ]
        removeModal()
      }, ignoreInit = TRUE)

      ###########################
      # Copy an analysis method #
      ###########################

      observeEvent(input$act_anameth_copy, {
        pos_select <- input$display_info_amameth_rows_selected
        if(isTruthy(pos_select)){
          # Update information on the analysis method
          dt_info_anameth <- rv_anameth$dt_info_anameth
          dt_info_anameth_copy <- copy(dt_info_anameth)
          dt_info_anameth_copy <- dt_info_anameth_copy[pos_select, ]
          anameth_select <- dt_info_anameth_copy$Name
          dt_info_anameth_copy$Name <- paste0(anameth_select, "_copy")
          dt_info_anameth <- rbind(dt_info_anameth, dt_info_anameth_copy)
          rv_anameth$dt_info_anameth <- dt_info_anameth

          # Make a copy of the selected analysis method
          dt_anameth_copy <- fread(paste0(path_pack, "/extdata/anameth/", anameth_select, ".csv"))
          write.csv(dt_anameth_copy, paste0(path_pack, "/extdata/anameth/", dt_info_anameth_copy$Name, ".csv"), row.names = FALSE)
        }else{
          showModal(modalDialog(title = "Error", "Select an analysis method!", easyClose = TRUE, footer = NULL))
        }
      })

      ##############################################
      # Display information on the analysis method #
      ##############################################

      output$display_info_amameth <- renderDataTable({
        dt_info_anameth <- rv_anameth$dt_info_anameth
        datatable(
          dt_info_anameth,
          colnames = c("Analysis method", "Kit", "Parameter file", "AT", "Allele probability"),
          selection = list(mode = "single", target = "row"),
          extensions = c("FixedColumns", "FixedHeader"),
          options = list(dom = "t", iDisplayLength = nrow(dt_info_anameth), ordering = FALSE, autoWidth = TRUE, scrollX = TRUE, scrollY = "650px",
                         fixedHeader = TRUE, fixedColumns = list(leftColumns = 1, rightColumns = 0),
                         columnDefs = list(list(targets = c(0, 2), width = "300px"),
                                           list(targets = c(1, 3, 4), width = "200px"))),
          rownames = FALSE
        )
      })
    }
  )
}
