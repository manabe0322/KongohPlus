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
                                        info_al_ui("info_al")
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
    bool_load_proj <- reactiveVal(FALSE)

    rv_input <- file_server("file", session, bool_load_proj)
    rv_deconvo <- deconvo_server("deconvo", path_pack, session, rv_input)
    rv_lr <- lr_server("lr", rv_input, rv_deconvo)
    rv_hd_true <- hd_true_server("hd_true", path_pack, rv_deconvo)

    anameth_server("anameth", path_pack)
    kit_server("kit", path_pack)
    info_al_server("info_al", path_pack)

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
      req(rv_input)
      req(rv_deconvo)
      req(rv_lr)
      req(rv_hd_true)

      data_deconvo <- rv_load_proj$data_deconvo
      if(!is.null(data_deconvo)){
        bool_load_proj(TRUE)
        rv_deconvo$data_deconvo <- data_deconvo
        rv_input$dt_csp <- data_deconvo$dt_csp
        rv_input$fn_csp <- data_deconvo$fn_csp
        rv_input$dt_af <- data_deconvo$dt_af
        rv_input$fn_af <- data_deconvo$fn_af
        rv_input$hnc_from <- min(data_deconvo$hnc_all)
        rv_input$hnc_to <- max(data_deconvo$hnc_all)
        rv_input$fn_anameth <- data_deconvo$fn_anameth
        rv_input$dt_ignore_mk <- data_deconvo$dt_ignore_mk
      }

      dt_ref <- rv_load_proj$dt_ref
      if(!is.null(dt_ref)){
        rv_input$dt_ref <- dt_ref
        rv_lr$dt_ref <- dt_ref
        rv_lr$fn_ref <- rv_load_proj$fn_ref
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
      req(rv_deconvo)
      req(rv_lr)
      req(rv_hd_true)
      save_proj_server("save_proj", rv_deconvo, rv_lr, rv_hd_true)
      report_server("report", path_pack, rv_deconvo)
    })
  }

  shinyApp(ui, server)
}
