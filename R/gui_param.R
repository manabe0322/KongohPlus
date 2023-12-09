#' est_par_load_ui
#'
#' @description The function to create the ui module for loading data to estimate parameters
est_par_load_ui <- function(id){
  ns <- NS(id)

  tabPanel("Load",
           br(),
           fluidRow(
             column(4,
                    fileInput(ns("exp_data"), label = h4("Experimental data"), accept = ".csv")
             ),
             column(4,
                    fileInput(ns("gt_data"), label = h4("Genotype data"), accept = ".csv"),
             ),
             column(4,
                    fileInput(ns("seq_data"), label = h4("Sequence data"), accept = ".csv")
             )
           ),
           br(),
           actionButton(ns("act_est_par_next"), label = "Next", class = "btn btn-primary btn-lg")
  )
}

est_par_setting_ui <- function(id){
  ns <- NS(id)

  tabPanel("Setting")
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
           tabsetPanel(
             est_par_load_ui(ns("est_par_load")),
             est_par_setting_ui(ns("est_par_setting")),
             est_par_result_ui(ns("est_par_result"))
           )
  )
}

#' est_par_load_server
#'
#' @description The function to create the server module for loading data to estimate parameters
est_par_load_server <- function(id){
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
          dt_gt <- change_data_type_exp(dt_gt)
          return(dt_gt)
        }
      })

      # Load sequence data
      react_dt_seq <- reactive({
        file_input <- input$seq_data
        if(is.null(file_input)){
          return(NULL)
        }else{
          dt_seq <- fread(file_input$datapath)
          dt_seq <- change_data_type_exp(dt_seq)
          return(dt_seq)
        }
      })

      observeEvent(input$act_est_par_next){
        # Check input files
      }
    }
  )
}

est_par_setting_server <- function(id){
  moduleServer(
    id,
    function(input, output, session){}
  )
}

est_par_result_server <- function(id){
  moduleServer(
    id,
    function(input, output, session){}
  )
}

#' est_par_server
#'
#' @description The function to create the server module for estimating parameters
est_par_server <- function(id){
  moduleServer(
    id,
    function(input, output, session){
      est_par_load_server("est_par_load")
      est_par_setting_server("est_par_setting")
      est_par_result_server("est_par_result")
    }
  )
}
