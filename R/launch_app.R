#' KongohPlus
#'
#' @description Main window
#' @usage KongohPlus()
#' @export
KongohPlus <- function(){
  ver_soft <- packageVersion("KongohPlus")
  path_pack <- path.package("KongohPlus", quiet = FALSE)
  options(shiny.maxRequestSize = 500 * 1024^2)

  ui <- fluidPage(theme = shinytheme("cerulean"),
                  useShinyjs(),
                  useShinyFeedback(),
                  tags$style(type = "text/css", "body{padding-top: 65px;}"),

                  navbarPage(title = paste0("KongohPlus ver. ", ver_soft),
                             id = "navbar",
                             position = c("fixed-top"),

                             file_ui("file"),

                             deconvo_ui("deconvo"),

                             lr_ui("lr"),

                             report_ui("report"),

                             navbarMenu("Project",
                                        tabPanel("New project",
                                                 useShinyjs(),
                                                 h2("New project"),
                                                 br(),
                                                 actionButton("act_new_proj", label = "New project", class = "btn btn-primary btn-lg")
                                        ),
                                        load_proj_ui("load_proj"),
                                        save_proj_ui("save_proj")
                             ),

                             navbarMenu("Setting",
                                        anameth_ui("anameth"),
                                        kit_ui("kit"),
                                        size_ui("size")
                             ),

                             navbarMenu("Experimental parameters",
                                        seq_data_ui("seq_data"),
                                        est_sr_ui("est_sr"),
                                        edit_sr_ui("edit_sr"),
                                        est_par_ui("est_par")
                             ),

                             navbarMenu("Tools",
                                        tabPanel("Example files"),
                                        hd_true_ui("hd_true"),
                                        tabPanel("Validation mode")
                             ),

                             tabPanel("Manual")
                  )
  )

  server <- function(input, output, session){
    rv_file <- file_server("file", session)
    rv_deconvo <- deconvo_server("deconvo", path_pack, session, rv_file)
    rv_lr <- lr_server("lr", rv_file, rv_deconvo)
    rv_hd_true <- hd_true_server("hd_true", path_pack, rv_deconvo)

    anameth_server("anameth", path_pack)
    kit_server("kit", path_pack)
    size_server("size", path_pack)

    seq_data_server("seq_data", path_pack)
    rv_est_sr <- est_sr_server("est_sr", path_pack)
    edit_sr_server("edit_sr", path_pack)
    est_par_server("est_par", path_pack)

    observeEvent(input$act_new_proj, {
      refresh()
    })

    rv_load_proj <- load_proj_server("load_proj", session)
    observe({
      req(rv_load_proj)
      req(rv_file)
      req(rv_deconvo)
      req(rv_lr)
      req(rv_hd_true)

      data_deconvo <- rv_load_proj$data_deconvo
      if(!is.null(data_deconvo)){
        rv_deconvo$data_deconvo <- data_deconvo
        rv_file$dt_csp <- data_deconvo$dt_csp
        rv_file$dt_af <- data_deconvo$dt_af
      }

      dt_ref <- rv_load_proj$dt_ref
      if(!is.null(dt_ref)){
        rv_file$dt_ref <- dt_ref
        rv_lr$dt_ref <- dt_ref
        rv_lr$fn_ref <- rv_load_proj$fn_ref
        rv_lr$refs <- rv_load_proj$refs
      }

      data_lr_list <- rv_load_proj$data_lr_list
      if(length(data_lr_list) != 0){
        rv_lr$data_lr_list <- data_lr_list
        rv_lr$dt_lr_calc_record <- rv_load_proj$dt_lr_calc_record
      }

      lr_rand <- rv_load_proj$lr_rand
      if(!is.null(lr_rand)){
        rv_hd_true$lr_rand <- lr_rand
      }
    })

    observe({
      req(rv_file)
      req(rv_deconvo)
      req(rv_lr)
      req(rv_hd_true)
      save_proj_server("save_proj", rv_file, rv_deconvo, rv_lr, rv_hd_true)
    })
  }

  shinyApp(ui, server)
}
