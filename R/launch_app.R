#' KongohPlus
#'
#' @description Main window
#' @param num_threads Number of threads for parallelization
#' @usage KongohPlus()
#' @export
KongohPlus <- function(num_threads = NULL){

  # --- Information on the package 'KongohPlus' --- #

  ver_soft <- packageVersion("KongohPlus")
  path_pack <- path.package("KongohPlus", quiet = FALSE)

  # --- Set the number of threads for parallelization --- #

  if(is.null(num_threads)){
    num_threads <- max(detectCores() - 2, 1)
  }else{
    num_threads <- max(round(num_threads, 0), 1)
  }
  Sys.setenv(OMP_NUM_THREADS = num_threads)

  # --- Maximum file size for shiny --- #

  options(shiny.maxRequestSize = 500 * 1024^2)

  ui <- fluidPage(theme = shinytheme("cerulean"),
                  useShinyjs(),
                  useShinyFeedback(),
                  tags$style(type = "text/css", "body{padding-top: 65px;}"),
                  tags$style(HTML(".gt_table { font-family: inherit !important; }")),

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
                                        kit_ui("kit"),
                                        info_al_ui("info_al"),
                                        est_sr_ui("est_sr"),
                                        est_par_ui("est_par"),
                                        est_dropin_ui("est_dropin"),
                                        anameth_ui("anameth")
                             ),

                             navbarMenu("Tools",
                                        example_file_ui("example_file"),
                                        view_sr_ui("view_sr"),
                                        view_lsae_ui("view_lsae"),
                                        view_gamma_ui("view_gamma"),
                                        view_dropin_ui("view_dropin"),
                                        hd_true_random_ui("hd_true_random"),
                                        hd_true_importance_ui("hd_true_importance")
                             ),

                             manual_ui("manual")
                  )
  )

  server <- function(input, output, session){

    # --- Analysis --- #

    rv_input <- file_server("file", session)
    rv_deconvo <- deconvo_server("deconvo", path_pack, session, rv_input)
    rv_lr <- lr_server("lr", rv_input, rv_deconvo)
    report_server("report", path_pack, rv_deconvo, rv_lr)

    # --- Tools --- #

    example_file_server("example_file", path_pack)
    view_sr_server("view_sr", path_pack)
    view_lsae_server("view_lsae", path_pack)
    view_gamma_server("view_gamma", path_pack)
    view_dropin_server("view_dropin", path_pack)
    rv_hd_true_random <- hd_true_random_server("hd_true_random", path_pack, rv_deconvo)
    rv_hd_true_importance <- hd_true_importance_server("hd_true_importance", path_pack, rv_deconvo)

    # --- Manual --- #

    manual_server("manual", ver_soft, path_pack)

    # --- Setting --- #

    kit_server("kit", path_pack)
    info_al_server("info_al", path_pack)
    rv_est_sr <- est_sr_server("est_sr", path_pack)
    est_par_server("est_par", path_pack)
    est_dropin_server("est_dropin", path_pack)
    anameth_server("anameth", path_pack)

    # --- New project --- #

    observeEvent(input$act_new_proj, {
      refresh()
    })

    # --- Load project --- #

    rv_load_proj <- load_proj_server("load_proj", session)
    observe({
      req(rv_load_proj)
      req(rv_input)
      req(rv_deconvo)
      req(rv_lr)
      req(rv_hd_true_random)
      req(rv_hd_true_importance)

      rv_deconvo$data_deconvo <- data_deconvo <- rv_load_proj$data_deconvo
      rv_deconvo$bool_perform_deconvo <- FALSE

      if(!is.null(data_deconvo)){
        rv_input$dt_csp <- data_deconvo$dt_csp
        rv_input$fn_csp <- data_deconvo$fn_csp
        rv_input$dt_af <- data_deconvo$dt_af
        rv_input$fn_af <- data_deconvo$fn_af
      }

      data_lr <- rv_load_proj$data_lr
      if(length(data_lr) != 0){
        rv_input$dt_ref <- rv_lr$dt_ref <- data_lr$dt_ref
        rv_input$fn_ref <- rv_lr$fn_ref <- data_lr$fn_ref
        rv_lr$refs <- data_lr$refs
        rv_lr$lr_list <- data_lr$lr_list
        rv_lr$dt_lr_calc_record <- data_lr$dt_lr_calc_record
      }

      hd_true_random <- rv_load_proj$hd_true_random
      if(!is.null(hd_true_random)){
        rv_hd_true_random$lr_rand <- hd_true_random$lr_rand
      }

      hd_true_importance <- rv_load_proj$hd_true_importance
      if(!is.null(hd_true_importance)){
        rv_hd_true_importance$lr_importance <- hd_true_importance$lr_importance
        rv_hd_true_importance$bias <- hd_true_importance$bias
        rv_hd_true_importance$expect_lr_hd <- hd_true_importance$expect_lr_hd
      }
    })

    # --- Save project --- #

    observe({
      req(rv_deconvo)
      req(rv_lr)
      req(rv_hd_true_random)
      req(rv_hd_true_importance)
      save_proj_server("save_proj", rv_deconvo, rv_lr, rv_hd_true_random, rv_hd_true_importance)
    })
  }

  shinyApp(ui, server)
}
