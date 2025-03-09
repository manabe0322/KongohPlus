#' save_proj_ui
#'
#' @description The function to create the ui module for saving the project
save_proj_ui <- function(id){
  ns <- NS(id)

  tabPanel("Save project",
           useShinyjs(),
           useWaiter(),
           h2("Save project"),
           br(),
           textInput(ns("name_proj"), label = "Enter the project name", value = NULL),
           br(),
           downloadButton(ns("download_proj"), label = "Save as", class = "btn btn-primary btn-lg")
  )
}

#' save_proj_server
#'
#' @description The function to create the server module for saving the project
save_proj_server <- function(id, rv_deconvo, rv_lr, rv_hd_true){
  moduleServer(
    id,
    function(input, output, session){
      data_deconvo <- rv_deconvo$data_deconvo
      dt_ref <- rv_lr$dt_ref
      fn_ref <- rv_lr$fn_ref
      data_lr_list <- rv_lr$data_lr_list
      dt_lr_calc_record <- rv_lr$dt_lr_calc_record
      lr_rand <- rv_hd_true$lr_rand

      output$download_proj <- downloadHandler(
        filename = function(){
          if(isTruthy(input$name_proj)){
            paste0(input$name_proj, ".RData")
          }else{
            paste0(gsub(" ", "_", format(as.POSIXct(Sys.time()), "%Y-%m-%d %H%M%S")), "_KongohPlus_project.RData")
          }
        },

        content = function(file){
          waiter_show(html = spin_3k(), color = "white")
          save(list = c("data_deconvo",
                        "dt_ref",
                        "fn_ref",
                        "data_lr_list",
                        "dt_lr_calc_record",
                        "lr_rand"),
               file = file)
          waiter_hide()
        }
      )
    }
  )
}
