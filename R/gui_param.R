#' est_par_setting_ui
#'
#' @description The function to create the ui module for setting to estimate parameters
est_par_setting_ui <- function(id){
  ns <- NS(id)

  tabPanel("Setting",
           br(),
           sidebarLayout(
             sidebarPanel(width = 3,
                          fileInput(ns("exp_data"), label = h4("Experimental data"), accept = ".csv"),
                          fileInput(ns("gt_data"), label = h4("Genotype data"), accept = ".csv"),
                          uiOutput(ns("seq_data")),
                          uiOutput(ns("kit_par")),
                          actionButton(ns("act_advanced"), "Advanced setting"),
                          br(),
                          br(),
                          actionButton(ns("act_analyze_est_par"), label = "Analyze", class = "btn btn-primary btn-lg")
             ),
             mainPanel(width = 9,
                       verbatimTextOutput("value1"),
                       dataTableOutput(ns("dt_par_setting")), tags$script(HTML('$(document).on("click", "input", function () {
                                                             var checkboxes = document.getElementsByName("BSR");
                                                             var checkboxesChecked = [];
                                                             for (var i=0; i<checkboxes.length; i++) {
                                                             if (checkboxes[i].checked) {
                                                             checkboxesChecked.push(checkboxes[i].value);
                                                             }
                                                             }
                                                             Shiny.onInputChange("checked_rows",checkboxesChecked);  })'))

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

      output$value1 <- renderPrint({ input$checked_rows })

      observeEvent(input$kit_par, ignoreInit = TRUE, {
        output$dt_par_setting <- renderDataTable(server = FALSE, {
          # Extract autosomal markers
          kit <- input$kit_par
          dt_kit_info <- fread(paste0(path_pack, "/extdata/kit/", kit, ".csv"))
          setkey(dt_kit_info, Sex_Chr)
          dt_kit_auto <- dt_kit_info[""]
          markers_par <- dt_kit_auto[, Marker]
          n_mk <- length(markers_par)

          # Define the default checked information
          bsr_checked <- rep('checked=\"checked\"', n_mk)
          fsr_checked <- rep('', n_mk)
          fsr_checked[markers_par == "D22S1045"] <- 'checked=\"checked\"'
          dsr_checked <- rep('', n_mk)

          # Make a data.table for kit-setting
          dt_kit_setting <- data.table(Marker = markers_par)
          dt_kit_setting$AT = glue('<input type="number" name="AT" value="{rep(50, n_mk)}"><br>')
          dt_kit_setting$BSR = glue('<input type="checkbox" name="BSR" {bsr_checked} value="{paste0("bsr_mk", 1:n_mk)}"><br>')
          dt_kit_setting$FSR = glue('<input type="checkbox" name="FSR" {fsr_checked} value="{paste0("fsr_mk", 1:n_mk)}"><br>')
          dt_kit_setting$DSR = glue('<input type="checkbox" name="DSR" {dsr_checked} value="{paste0("dsr_mk", 1:n_mk)}"><br>')

          datatable(
            dt_kit_setting,
            colnames = c("Marker", "AT (RFU)", "Marker-specific BSR", "Marker-specific FSR", "Marker-specific DSR"),
            selection = "none",
            options = list(iDisplayLength = 50, ordering = FALSE, autowidth = FALSE),
            rownames = FALSE,
            escape = FALSE
          )
        })
      })

      observeEvent(input$act_analyze_est_par, {
        print(input$checked_rows)
        updateTabsetPanel(session = parent_session, inputId = "tabset_est_par", selected = "Result")
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
