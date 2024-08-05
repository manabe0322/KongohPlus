#' KongohPlus
#'
#' @description Main window
#' @usage KongohPlus()
#' @export
KongohPlus <- function(){
  ver_soft <- packageVersion("KongohPlus")
  path_pack <- path.package("KongohPlus", quiet = FALSE)
  options(shiny.maxRequestSize = 500 * 1024^2)

  # Install and activate pandoc
  pandoc_install()
  pandoc_activate()

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
                                        sr_ui("sr"),
                                        kit_ui("kit"),
                                        size_ui("size")
                             ),

                             navbarMenu("Tools",
                                        tabPanel("Example files"),
                                        est_par_ui("est_par"),
                                        tabPanel("Validation mode")
                             ),

                             tabPanel("Manual")
                  )
  )

  server <- function(input, output, session){
    rv_file <- file_server("file", session)
    rv_deconvo <- deconvo_server("deconvo", path_pack, rv_file)
    rv_lr <- lr_server("lr", path_pack, rv_file, rv_deconvo)

    est_par_server("est_par", path_pack)
    anameth_server("anameth", path_pack)
    sr_server("sr", path_pack)
    kit_server("kit", path_pack)
    size_server("size", path_pack)

    observeEvent(input$act_new_proj, {
      refresh()
    })

    rv_load_proj <- load_proj_server("load_proj", session)
    observe({
      req(rv_load_proj)

      rv_file$dt_csp <- rv_load_proj$dt_csp
      rv_file$dt_ref <- rv_load_proj$dt_ref
      rv_file$dt_af <- rv_load_proj$dt_af

      data_deconvo <- rv_load_proj$data_deconvo
      if(!is.null(data_deconvo)){
        rv_deconvo$data_deconvo <- data_deconvo
      }

      data_lr <- rv_load_proj$data_lr
      if(!is.null(data_lr)){
        rv_lr$data_lr <- data_lr
        rv_lr$data_deconvo_hp <- rv_load_proj$data_deconvo_hp
        rv_lr$data_deconvo_hd <- rv_load_proj$data_deconvo_hd
      }
    })

    observe({
      req(rv_file)
      req(rv_deconvo)
      req(rv_lr)
      save_proj_server("save_proj", rv_file, rv_deconvo, rv_lr)
    })
  }

  shinyApp(ui, server)
}
