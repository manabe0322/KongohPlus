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
                          selectInput(ns("seq_data"), label = h4("Sequence data"), choices = c("sequence_data_default.csv"), selected = c("sequence_data_default.csv"), multiple = FALSE),
                          uiOutput(ns("kit_par")),
                          actionButton(ns("act_advanced"), "Advanced setting")
             ),
             mainPanel(width = 9,
                       dataTableOutput(ns("dt_par_setting"))
             )
           ),
           actionButton(ns("act_analyze_est_par"), label = "Analyze", class = "btn btn-primary btn-lg")
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

      output$kit_par <- renderUI({
        ns <- session$ns
        selectInput(ns("kit_par"), label = h4("Kit"), choices = c("GlobalFiler", "Identifiler"), selected = "GlobalFiler", multiple = FALSE)
      })

      output$dt_par_setting <- renderDataTable({
        datatable(
          data.table(Marker = c("A", "B", "C"),
                     BSR = c("Allele", "LUS", "Multi-seq")),
          rownames = FALSE
        )
      })

      observeEvent(input$act_analyze_est_par, {
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
