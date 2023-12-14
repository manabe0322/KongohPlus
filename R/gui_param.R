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

  tabPanel("Result")
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
      checkboxColumn <- function(len, col, ...){
        inputs <- character(len)
        for(i in seq_len(len)){
          inputs[i] <- as.character(
            checkboxInput(paste0("checkb_", col, "_", i), label = NULL, ...)
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

      js <- function(dtid, col_numeric, cols_check, ns = identity){
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
            sprintf("$('body').on('click', '[id^=checkb_%d_]', function() {", col),
            "  var id = this.getAttribute('id');",
            sprintf("  var i = parseInt(/checkb_%d_(\\d+)/.exec(id)[1]);", col),
            "  var value = $(this).prop('checked');",
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

      rv_dt_kit_set <- reactiveVal()

      react_dt_kit_set0 <- reactive({
        # Extract autosomal markers
        kit <- input$kit_par
        dt_kit_info <- fread(paste0(path_pack, "/extdata/kit/", kit, ".csv"))
        setkey(dt_kit_info, Sex_Chr)
        dt_kit_auto <- dt_kit_info[""]
        markers_par <- dt_kit_auto[, Marker]
        n_mk <- length(markers_par)

        # Make a data.table for kit-setting
        dt_kit_set0 <- data.table(Marker = markers_par)
        return(dt_kit_set0)
      })

      react_dt_kit_set1 <- reactive({
        dt_kit_set1 <- react_dt_kit_set0()
        dt_kit_set1$AT <- 50
        dt_kit_set1$Bool_BSR <- FALSE
        dt_kit_set1$Bool_FSR <- FALSE
        dt_kit_set1$Bool_DSR <- FALSE
        return(dt_kit_set1)
      })

      react_dt_kit_set2 <- reactive({
        dt_kit_set2 <- react_dt_kit_set0()
        n_mk <- nrow(dt_kit_set2)
        dt_kit_set2$numeric_at <- numericColumn(n_mk, 2, value = 50, width = "50%")
        dt_kit_set2$check_bsr <- checkboxColumn(n_mk, 3)
        dt_kit_set2$check_fsr <- checkboxColumn(n_mk, 4)
        dt_kit_set2$check_dsr <- checkboxColumn(n_mk, 5)
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

      observeEvent(input$kit_par, ignoreInit = TRUE, {
        rv_dt_kit_set(react_dt_kit_set1())

        output$dt_kit_set <- renderDataTable(server = FALSE, {
          datatable(
            isolate(react_dt_kit_set2()),
            colnames = c("Marker", "AT", "Marker-specific BSR", "Marker-specific FSR", "Marker-specific DSR"),
            selection = "none",
            options = list(iDisplayLength = 50, ordering = FALSE, autowidth = FALSE),
            rownames = FALSE,
            escape = FALSE,
            callback = JS(js("dt_kit_set", 2, c(3, 4, 5)))
          )
        })
      })

      observeEvent(input[["dt_kit_set_cell_edit"]], {
        info <- input[["dt_kit_set_cell_edit"]]
        options(warn = -1) # To suppress the message: Warning in coerceValue(v, data[i, j, drop = TRUE]) : The data type is not supported: logical
        rv_dt_kit_set(editData(rv_dt_kit_set(), info))
        options(warn = 0)
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
      # Check input files
      #error_message <- check_file_est_par(dt_exp, dt_gt, dt_seq, path_pack)
      #if(error_message != ""){
      #  showModal(modalDialog(title = "Error", error_message, easyClose = TRUE, footer = NULL))
      #}else{
      #updateTabsetPanel(session = parent_session, inputId = "tabset_est_par", selected = "est_par_setting")
      #}
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
