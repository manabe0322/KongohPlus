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
                       dataTableOutput(ns("dt_info_anameth"))
             )
           )
  )
}

#' anameth_server
#'
#' @description The function to create the server module for the analysis method
anameth_server <- function(id, path_pack){
  moduleServer(
    id,
    function(input, output, session){

      numericColumn <- function(len, col, ...){
        inputs <- character(len)
        for(i in seq_len(len)){
          inputs[i] <- as.character(
            numericInput(paste0("numeric_", col, "_", i), label = NULL, ...)
          )
        }
        return(inputs)
      }

      ######################################
      # Define the initial reactive values #
      ######################################

      rv_anameth <- reactiveValues()
      rv_anameth$dt_info_anameth <- create_dt_info_anameth(path_pack)
      rv_anameth$kit_par <- NULL
      rv_anameth$dt_at_new <- NULL

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
                               uiOutput(session$ns("anameth_new_param")),
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
                               numericInput(session$ns("anameth_new_gt_fltr_m2sr"), label = h4("Maximum minus 2-nt stutter ratio"), value = 0.05)
                      ),
                      tabPanel("Drop-in",
                               br(),
                               numericInput(session$ns("anameth_new_di_prob"), label = h4("Probability of drop-in"), value = 0.00335),
                               br(),
                               numericInput(session$ns("anameth_new_di_alpha"), label = h4("Shape parameter"), value = 1.22),
                               br(),
                               numericInput(session$ns("anameth_new_di_beta"), label = h4("Rate parameter"), value = 0.0503),
                               br(),
                               numericInput(session$ns("anameth_new_di_min_h"), label = h4("Minimum drop-in peak height"), value = 50),
                               br(),
                               numericInput(session$ns("anameth_new_di_max_h"), label = h4("Maximum drop-in peak height"), value = 300),
                               br(),
                               numericInput(session$ns("anameth_new_di_max_n"), label = h4("Maximum number of drop-in peaks"), value = 2)
                      )
          ),
          footer = tagList(
            actionButton(session$ns("act_new_anameth_save"), "Save"),
            modalButton("Cancel")
          ),
          size = "l"
        ))
      }, ignoreInit = TRUE)

      output$anameth_new_param <- renderUI({
        param_all <- list.files(paste0(path_pack, "/extdata/param"))
        param_all <- gsub(".{6}$", "", param_all)
        selectInput(session$ns("anameth_new_param"), label = h4("Parameter data"), choices = param_all, selected = param_all[1], selectize = FALSE)
      })

      observeEvent(input$anameth_new_param, {
        anameth_new_param <- input$anameth_new_param
        if(is.null(anameth_new_param)){
          output$anameth_new_kit <- renderText({"NA"})
        }else{
          load(paste0(path_pack, "/extdata/param/", anameth_new_param, ".Rdata"))

          # Set the kit name
          rv_anameth$kit_par <- kit_par <- dt_par$Kit[1]
          output$anameth_new_kit <- renderText({paste0(kit_par)})

          # Extract autosomal markers
          dt_kit_select <- fread(paste0(path_pack, "/extdata/kit/", kit_par, ".csv"))
          markers <- unique(dt_kit_select$Marker[dt_kit_select$Chr == "Auto"])
          n_mk <- length(markers)

          # Create a data.table for the AT
          dt_at_new <- data.table(Marker = markers) # For saved data
          dt_at_new_display <- copy(dt_at_new) # For displayed data
          dt_at_new$AT <- 100
          rv_anameth$dt_at_new <- dt_at_new
          dt_at_new_display$AT <- numericColumn(n_mk, 2, value = 100, width = "75%")

          # Display a data.table for the AT
          output$anameth_new_at <- renderDataTable(server = FALSE, {
            datatable(
              dt_at_new_display,
              selection = "none",
              options = list(dom = "t", iDisplayLength = nrow(dt_at_new_display), ordering = FALSE, autoWidth = TRUE),
              rownames = FALSE,
              escape = FALSE
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

      observeEvent(input$act_new_anameth_save, {
        # Extract settings
        anameth_new_input_name <- input$anameth_new_input_name
        anameth_new_param <- input$anameth_new_param
        kit_par <- rv_anameth$kit_par
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
        anameth_new_di_prob <- input$anameth_new_di_prob
        anameth_new_di_alpha <- input$anameth_new_di_alpha
        anameth_new_di_beta <- input$anameth_new_di_beta
        anameth_new_di_min_h <- input$anameth_new_di_min_h
        anameth_new_di_max_h <- input$anameth_new_di_max_h
        anameth_new_di_max_n <- input$anameth_new_di_max_n

        # Create a data.table for the new analysis method
        dt_anameth_new <- data.table(Category = c("Parameter", "Kit", rep("AT", length(markers)), rep("Allele_probability", 3), rep("Genotype_filter", 4), rep("Drop-in", 6)),
                                     Subcategory = c("", "", markers, "Method", "Maf", "Theta", "BSR", "FSR", "DSR", "M2SR", "Probability", "Shape", "Rate", "Min_height", "Max_height", "Max_number"),
                                     Value = c(anameth_new_param, kit_par, ats,
                                               anameth_new_alprob, anameth_new_maf, anameth_new_theta,
                                               anameth_new_gt_fltr_bsr, anameth_new_gt_fltr_fsr, anameth_new_gt_fltr_dsr, anameth_new_gt_fltr_m2sr,
                                               anameth_new_di_prob, anameth_new_di_alpha, anameth_new_di_beta, anameth_new_di_min_h, anameth_new_di_max_h, anameth_new_di_max_n))
        write.csv(dt_anameth_new, paste0(path_pack, "/extdata/anameth/", anameth_new_input_name, ".csv"), row.names = FALSE)

        # Update information on the analysis methods
        rv_anameth$dt_info_anameth <- create_dt_info_anameth(path_pack)

        removeModal()
      })





      output$dt_info_anameth <- renderDataTable({
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
