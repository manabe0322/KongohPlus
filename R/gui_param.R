#' est_par_setting_ui
#'
#' @description The function to create the ui module for setting to estimate parameters
est_par_setting_ui <- function(id){
  ns <- NS(id)

  tabPanel("Setting",
           br(),

           tags$style(
             HTML(
               "td > div.form-group {
                  margin: 0px;
                  height: 30px;
                  display: flex;
                  align-items: center;
                }
                td:nth-child(1) {
                  width: 100px;
                }
                td:nth-child(2) {
                  width: 50px;
                }
                td:nth-child(3) {
                  width: 100px;
                }
                td:nth-child(4) {
                  width: 150px;
                }
                td:nth-child(5) {
                  width: 150px;
                }
                td:nth-child(6) {
                  width: 150px;
                }
                td:nth-child(7) {
                  width: 150px;
                }"
             )
           ),

           sidebarLayout(
             sidebarPanel(width = 3,
                          fileInput(ns("exp_data"), label = h4("Experimental data"), accept = ".csv"),
                          fileInput(ns("gt_data"), label = h4("Genotype data"), accept = ".csv"),
                          uiOutput(ns("seq_data")),
                          uiOutput(ns("kit_par")),
                          br(),
                          actionButton(ns("act_advanced"), "Advanced setting"),
                          br(),
                          br(),
                          br(),
                          actionButton(ns("act_analyze_est_par"), label = "Analyze", class = "btn btn-primary btn-lg")
             ),
             mainPanel(width = 9,
                       dataTableOutput(ns("dt_kit_set"))
             )
           )
  )
}

est_par_result_ui <- function(id){
  ns <- NS(id)

  tabPanel("Result",
           br(),
           tabsetPanel(id = ns("tabset_est_par_result"),
                       type = "tabs",
                       tabPanel("AE",
                                sidebarLayout(
                                  sidebarPanel(width = 2,
                                               actionButton(ns("act_trace_ae"), label = "Trace plot"),
                                               br(),
                                               br(),
                                               br(),
                                               actionButton(ns("act_plot_ae"), label = "Data plot")
                                  ),
                                  mainPanel(width = 10,
                                            dataTableOutput(ns("dt_result_ae"))
                                  )
                                )
                       ),
                       tabPanel("BSR"),
                       tabPanel("FSR"),
                       tabPanel("DSR"),
                       tabPanel("M2SR")
           )
  )
}

#' est_par_ui
#'
#' @description The function to create the ui module for estimating parameters
est_par_ui <- function(id){
  ns <- NS(id)

  tabPanel("Estimate parameters",
           titlePanel("Estimate parameters"),
           br(),
           tabsetPanel(id = ns("tabset_est_par"),
                       type = "tabs",
                       est_par_setting_ui(ns("est_par_setting")),
                       est_par_result_ui(ns("est_par_result"))
           )
  )
}

#' est_par_setting_server
#'
#' @description The function to create the server module for setting to estimate parameters
est_par_setting_server <- function(id, path_pack, parent_session){
  moduleServer(
    id,
    function(input, output, session){
      selectColumn <- function(len, col, choices, select_default, ...){
        inputs <- character(len)
        for(i in seq_len(len)){
          inputs[i] <- as.character(
            selectizeInput(paste0("select_", col, "_", i), label = NULL, choices = choices, selected = select_default[i], ...)
          )
        }
        return(inputs)
      }

      numericColumn <- function(len, col, ...){
        inputs <- character(len)
        for(i in seq_len(len)){
          inputs[i] <- as.character(
            numericInput(paste0("numeric_", col, "_", i), label = NULL, ...)
          )
        }
        return(inputs)
      }

      js <- function(dtid, col_numeric, cols_check){
        code <- vector("list", 1 + length(cols_check))
        code[[1]] <- c(
          sprintf("$('body').on('change', '[id^=numeric_%d_]', function() {", col_numeric),
          "  var id = this.getAttribute('id');",
          sprintf("  var i = parseInt(/numeric_%d_(\\d+)/.exec(id)[1]);", col_numeric),
          "  var value = $(this).val();",
          sprintf("  var info = [{row: i, col: %d, value: value}];", col_numeric),
          sprintf("  Shiny.setInputValue('%s', info);", session$ns(sprintf("%s_cell_edit:DT.cellInfo", dtid))),
          "});"
        )

        for(i in seq_along(cols_check)){
          col <- cols_check[i]
          code[[i + 1]] <- c(
            sprintf("$('body').on('change', '[id^=select_%d_]', function() {", col),
            "  var id = this.getAttribute('id');",
            sprintf("  var i = parseInt(/select_%d_(\\d+)/.exec(id)[1]);", col),
            "  var value = $(this).val();",
            sprintf("  var info = [{row: i, col: %d, value: value}];", col),
            sprintf("  Shiny.setInputValue('%s', info);", session$ns(sprintf("%s_cell_edit:DT.cellInfo", dtid))),
            "});"
          )
        }
        do.call(c, code)
      }

      # Load experimental data
      react_dt_exp <- reactive({
        file_input <- input$exp_data
        if(is.null(file_input)){
          return(NULL)
        }else{
          dt_exp <- fread(file_input$datapath)
          dt_exp <- change_data_type_exp(dt_exp)
          return(dt_exp)
        }
      })

      # Load genotype data
      react_dt_gt <- reactive({
        file_input <- input$gt_data
        if(is.null(file_input)){
          return(NULL)
        }else{
          dt_gt <- fread(file_input$datapath)
          dt_gt <- change_data_type_gt(dt_gt)
          return(dt_gt)
        }
      })

      # Load sequence data
      react_dt_seq <- reactive({
        file_input <- input$seq_data
        if(is.null(file_input)){
          return(NULL)
        }else{
          dt_seq <- fread(paste0(path_pack, "/extdata/sequence_data/", file_input))
          dt_seq <- change_data_type_seq(dt_seq)
          return(dt_seq)
        }
      })

      # Kit setting
      rv_dt_kit_set <- reactiveVal()

      react_dt_kit_set0 <- reactive({
        # Extract autosomal markers
        kit <- input$kit_par
        dt_kit_info <- fread(paste0(path_pack, "/extdata/kit/", kit, ".csv"))
        setkey(dt_kit_info, Sex_Chr)
        dt_kit_auto <- dt_kit_info[""]
        markers_par <- dt_kit_auto[, Marker]
        dyes_par <- dt_kit_auto[, Dye]

        # Make a data.table for kit-setting
        dt_kit_set0 <- data.table(Marker = markers_par, Dye = dyes_par)
        return(dt_kit_set0)
      })

      react_dt_kit_set1 <- reactive({
        dt_kit_set1 <- react_dt_kit_set0()
        dt_kit_set1$AT <- 50
        dt_kit_set1$Model_BSR <- "Locus-specific"
        dt_kit_set1$Model_FSR <- "Multiple loci together"
        dt_kit_set1$Model_DSR <- "Multiple loci together"
        dt_kit_set1$Model_M2SR <- "Not consider"
        return(dt_kit_set1)
      })

      react_dt_kit_set2 <- reactive({
        dt_kit_set2 <- react_dt_kit_set0()
        n_mk <- nrow(dt_kit_set2)
        dt_kit_set2$numeric_at <- numericColumn(n_mk, 3, value = 50, width = "75%")
        dt_kit_set2$select_bsr <- selectColumn(n_mk, 4, choices = c("Not consider", "Locus-specific", "Multiple loci together"), select_default = rep("Locus-specific", n_mk), width = "90%")
        dt_kit_set2$select_fsr <- selectColumn(n_mk, 5, choices = c("Not consider", "Locus-specific", "Multiple loci together"), select_default = rep("Multiple loci together", n_mk), width = "90%")
        dt_kit_set2$select_dsr <- selectColumn(n_mk, 6, choices = c("Not consider", "Locus-specific", "Multiple loci together"), select_default = rep("Multiple loci together", n_mk), width = "90%")
        dt_kit_set2$select_m2sr <- selectColumn(n_mk, 7, choices = c("Not consider", "Locus-specific", "Multiple loci together"), select_default = rep("Not consider", n_mk), width = "90%")
        return(dt_kit_set2)
      })

      output$seq_data <- renderUI({
        ns <- session$ns
        seq_files <- list.files(paste0(path_pack, "/extdata/sequence_data"))
        selectInput(ns("seq_data"), label = h4("Sequence data"), choices = seq_files, selected = seq_files[1], multiple = FALSE)
      })

      output$kit_par <- renderUI({
        ns <- session$ns
        kit_names <- list.files(paste0(path_pack, "/extdata/kit"))
        kit_names <- gsub(".csv", "", kit_names)
        selectInput(ns("kit_par"), label = h4("Kit"), choices = kit_names, selected = kit_names[1], multiple = FALSE)
      })

      observeEvent(input$act_advanced, {
        showModal(modalDialog(
          title = "Advanced setting",
          fluidRow(
            column(6,
                   h4("MCMC conditions"),
                   numericInput(session$ns("warmup"), label = "Warmup iterations", value = 2000),
                   numericInput(session$ns("sampling"), label = "Sampling iterations", value = 2000),
                   numericInput(session$ns("chains"), label = "Number of chains", value = 4)
            ),
            column(6,
                   h4("Pull-up"),
                   checkboxInput(session$ns("pu_auto"), label = "Automatically removed", value = TRUE),
                   numericInput(session$ns("pu_dist"), label = "Maximum pull-up distance (bp)", value = 1),
                   h4("Stutter model"),
                   numericInput(session$ns("motif_lus"), label = "Minimum number of motifs for considering the LUS model", value = 2),
                   numericInput(session$ns("motif_lus"), label = "Minimum number of motifs for considering the multi-sequence model", value = 3)
            )
          ),
          footer = tagList(
            disabled(actionButton(session$ns("act_advance_save"), "Save")),
            modalButton("Cancel")
          ),
          size = "l"
        ))
      })

      observeEvent(input$kit_par, ignoreInit = TRUE, {
        rv_dt_kit_set(react_dt_kit_set1())

        output$dt_kit_set <- renderDataTable(server = FALSE, {
          datatable(
            isolate(react_dt_kit_set2()),
            colnames = c("Marker", "Dye", "AT", "BSR model", "FSR model", "DSR model", "M2SR model"),
            selection = "none",
            options = list(iDisplayLength = 50, ordering = FALSE, autowidth = TRUE),
            rownames = FALSE,
            escape = FALSE,
            callback = JS(js("dt_kit_set", 3, c(4, 5, 6, 7)))
          )
        })
      })

      observeEvent(input[["dt_kit_set_cell_edit"]], {
        info <- input[["dt_kit_set_cell_edit"]]
        dt_kit_set <- rv_dt_kit_set()
        dt_kit_set[[info[1, "col"]]][info[1, "row"]] <- info[1, "value"]
        rv_dt_kit_set(dt_kit_set)

        print(rv_dt_kit_set())
      })

      observeEvent(input$act_analyze_est_par, {
        # Fix data.table
        dt_exp <- react_dt_exp()
        dt_gt <- react_dt_gt()
        dt_seq <- react_dt_seq()
        dt_kit_set <- rv_dt_kit_set()

        # Check input files and setting
        error_message <- check_est_par(dt_exp, dt_gt, dt_seq, dt_kit_set)
        if(error_message != ""){
          showModal(modalDialog(title = "Error", HTML(error_message), easyClose = TRUE, footer = NULL))
        }else{
          results_est_par <- estimate_parameters(dt_exp, dt_gt, dt_seq, dt_kit_set, n_warmup, n_sampling, n_chains, n_parallel_chains)

          updateTabsetPanel(session = parent_session, inputId = "tabset_est_par", selected = "Result")
        }
      })

      react_list <- list(react_dt_exp = react_dt_exp, react_dt_gt = react_dt_gt)
      return(react_list)
    }
  )
}

est_par_result_server <- function(id){
  moduleServer(
    id,
    function(input, output, session){
    }
  )
}

#' est_par_server
#'
#' @description The function to create the server module for estimating parameters
est_par_server <- function(id, path_pack){
  moduleServer(
    id,
    function(input, output, session){
      react_list <- est_par_setting_server("est_par_setting", path_pack, session)
      est_par_result_server("est_par_result")
    }
  )
}
