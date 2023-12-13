#' est_par_setting_ui
#'
#' @description The function to create the ui module for setting to estimate parameters
est_par_setting_ui <- function(id){
  ns <- NS(id)

  tabPanel("Setting",
           br(),

           #           tags$style(
           #             HTML(
           #               "td > div.form-group > div.checkbox {
           #                  margin:0px;
           #                  bottom:10px;
           #                }
           #                td > div.form-group {
           #                  position:absolute;
           #                }"
           #             )
           #           ),

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
        n_mk <- nrow(dt_kit_set1)
        at_html <- bool_bsr_html <- bool_fsr_html <- bool_dsr_html <- rep("", n_mk)
        for(i in 1:n_mk){
          at_html[i] <- paste0("<input id='", session$ns(paste0("at_", i)), "' class='shiny-bound-input' value = '", 50, "' type='number' step = '1'>")
          bool_bsr_html[i] <- paste0("<input id='", session$ns(paste0("check_bsr_", i)), "' class='shiny-bound-input' value = '", i, "' type='checkbox'>")
          bool_fsr_html[i] <- paste0("<input id='", session$ns(paste0("check_bsr_", i)), "' class='shiny-bound-input' value = '", i, "' type='checkbox'>")
          bool_dsr_html[i] <- paste0("<input id='", session$ns(paste0("check_bsr_", i)), "' class='shiny-bound-input' value = '", i, "' type='checkbox'>")
        }
        dt_kit_set1$AT <- at_html
        dt_kit_set1$Bool_BSR <- bool_bsr_html
        dt_kit_set1$Bool_FSR <- bool_fsr_html
        dt_kit_set1$Bool_DSR <- bool_dsr_html
        return(dt_kit_set1)
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
            isolate(rv_dt_kit_set()),
            colnames = c("Marker", "AT", "Marker-specific BSR", "Marker-specific FSR", "Marker-specific DSR"),
            selection = "none",
            options = list(iDisplayLength = 50, ordering = FALSE, autowidth = FALSE,
                           drawCallback = JS('function(settings) {
                                             Shiny.bindAll(this.api().table().node());}')),
            rownames = FALSE,
            escape = FALSE
          )
        })
      })

      observeEvent(input$act_analyze_est_par, {
        isolate(dt_kit_set <- rv_dt_kit_set())
        markers <- dt_kit_set[, Marker]
        n_mk <- length(markers)
        at_all_mk <- rep(0, n_mk)
        bool_bsr_all_mk <- bool_fsr_all_mk <- bool_dsr_all_mk <- rep(FALSE, n_mk)

        for(i in 1:n_mk){
          at_all_mk[i] <- input[[paste0("at_", i)]]
          bool_bsr_all_mk[i] <- input[[paste0("check_bsr_", i)]]
          bool_fsr_all_mk[i] <- input[[paste0("check_fsr_", i)]]
          bool_dsr_all_mk[i] <- input[[paste0("check_dsr_", i)]]
        }

        dt_test <- data.table(Marker = markers,
                              AT = at_all_mk,
                              Bool_BSR = bool_bsr_all_mk,
                              Bool_FSR = bool_fsr_all_mk,
                              Bool_DSR = bool_dsr_all_mk)

        print(dt_test)
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
