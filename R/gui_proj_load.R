#' load_proj_ui
#'
#' @description The function to create the ui module for loading the project
load_proj_ui <- function(id){
  ns <- NS(id)

  tabPanel("Load project",
           useWaiter(),
           h2("Load project"),
           br(),
           fileInput(ns("file_proj"), label = "Select a project file", accept = ".RData"),
           br(),
           actionButton(ns("act_load_proj"), label = "Load project", class = "btn btn-primary btn-lg")
  )
}

#' load_proj_server
#'
#' @description The function to create the server module for loading the project
load_proj_server <- function(id, session_top){
  moduleServer(
    id,
    function(input, output, session){
      rv_load_proj <- reactiveValues()
      rv_load_proj$data_deconvo <- NULL
      rv_load_proj$data_lr <- NULL
      rv_load_proj$data_deconvo_hp <- NULL
      rv_load_proj$data_deconvo_hd <- NULL
      rv_load_proj$lr_rand <- NULL

      load_proj <- reactive({
        file_input <- input$file_proj
        if(is.null(file_input)){
          return(NULL)
        }else{
          return("OK")
        }
      })

      observeEvent(input$act_load_proj, {
        proj <- load_proj()
        if(is.null(proj)){
          showModal(modalDialog(title = "Error", "Select a project file!", easyClose = TRUE, footer = NULL))
        }else{
          waiter_show(html = spin_3k(), color = "white")
          load(input$file_proj$datapath)
          rv_load_proj$data_deconvo <- data_deconvo
          rv_load_proj$data_lr <- data_lr
          rv_load_proj$data_deconvo_hp <- data_deconvo_hp
          rv_load_proj$data_deconvo_hd <- data_deconvo_hd
          rv_load_proj$lr_rand <- lr_rand

          if(is.null(data_deconvo)){
            if(is.null(data_lr)){
              updateNavbarPage(session = session_top, "navbar", selected = "File")
            }else{
              updateNavbarPage(session = session_top, "navbar", selected = "Likelihood ratio")
            }
          }else{
            if(is.null(data_lr)){
              updateNavbarPage(session = session_top, "navbar", selected = "Deconvolution")
            }else{
              updateNavbarPage(session = session_top, "navbar", selected = "Likelihood ratio")
            }
          }
          waiter_hide()
        }
      }, ignoreInit = TRUE)

      return(rv_load_proj)
    }
  )
}
