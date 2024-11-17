#' est_sr_setting_ui
#'
#' @description The function to create the ui module to create the setting tab for estimating average stutter ratios
est_sr_setting_ui <- function(id){
  ns <- NS(id)

  tabPanel("Setting",
           br(),
           sidebarLayout(
             sidebarPanel(width = 3,
                          fileInput(ns("exp_data"), label = h4("Experimental data"), accept = ".csv", buttonLabel = "Load"),
                          fileInput(ns("gt_data"), label = h4("Genotype data"), accept = ".csv", buttonLabel = "Load"),
                          uiOutput(ns("seq_data")),
                          uiOutput(ns("kit_par")),
                          br(),
                          actionButton(ns("act_advance_par"), "Advanced setting"),
                          br(),
                          br(),
                          br(),
                          actionButton(ns("act_analyze_est_par"), label = "Analyze", class = "btn btn-primary btn-lg")
             ),
             mainPanel(width = 9,
                       tags$style(HTML("
                         .handsontable .htCore td {
                           vertical-align: middle !important;
                           text-align: center !important;
                         }"
                       )),
                       rHandsontableOutput(ns("dt_setting_est_sr")),
                       br()
             )
           )
  )
}

#' create_init_dt_setting_est_sr
#'
#' @description The function to create the initial data.table in the setting tab for estimating average stutter ratios
#' @param kit_name A kit name
#' @param path_pack The package path of KongohPlus
create_init_dt_setting_est_sr <- function(kit_name, path_pack){
  dt_kit <- fread(paste0(path_pack, "/extdata/kit/", kit_name, ".csv"))
  dt_kit <- change_data_type_kit(dt_kit)
  dt_kit_auto <- extract_kit_auto(dt_kit)
  markers <- unique(dt_kit_auto$Marker)
  dyes <- extract_dye(dt_kit_auto, markers)

  dt_setting_est_sr <- data.table(Marker = markers, Dye = dyes)
  dt_setting_est_sr$AT <- 50
  dt_setting_est_sr$Model_BSR <- "Allele"
  dt_setting_est_sr$Model_FSR <- "Uniform (Multiple loci)"
  dt_setting_est_sr$Model_DSR <- "Uniform (Multiple loci)"
  dt_setting_est_sr$Model_M2SR <- "Not consider"
  dt_setting_est_sr$X_Multi <- NA

  return(dt_setting_est_sr)
}

#' est_sr_setting_server
#'
#' @description The function to create the server module to create the setting tab for estimating average stutter ratios
est_sr_setting_server <- function(id, path_pack){
  moduleServer(
    id,
    function(input, output, session){
      rv_est_sr <- reactiveValues()
      rv_est_sr$dt_exp <- NULL
      rv_est_sr$dt_gt <- NULL
      rv_est_sr$dt_seq <- NULL
      rv_est_sr$dt_kit <- NULL
      rv_est_sr$dt_setting_est_sr <- NULL
      rv_est_sr$results_est_sr <- NULL

      # Load experimental data
      observeEvent(input$exp_data, {
        file_input <- input$exp_data
        if(!is.null(file_input)){
          dt_exp <- fread(file_input$datapath)
          dt_exp <- change_data_type_exp(dt_exp)
          rv_est_sr$dt_exp <- dt_exp
        }
      })

      # Load genotype data
      observeEvent(input$gt_data, {
        file_input <- input$gt_data
        if(!is.null(file_input)){
          dt_gt <- fread(file_input$datapath)
          dt_gt <- change_data_type_gt(dt_gt)
          rv_est_sr$dt_gt <- dt_gt
        }
      })

      # Select sequence data
      output$seq_data <- renderUI({
        ns <- session$ns
        seq_files <- list.files(paste0(path_pack, "/extdata/seq"))
        selectInput(ns("seq_data"), label = h4("Sequence data"), choices = seq_files, selected = seq_files[1], multiple = FALSE, selectize = FALSE)
      })

      observeEvent(input$seq_data, {
        file_input <- input$seq_data
        if(!is.null(file_input)){
          dt_seq <- fread(paste0(path_pack, "/extdata/seq/", file_input))
          dt_seq <- change_data_type_seq(dt_seq)
          rv_est_sr$dt_seq <- dt_seq
        }
      })

      # Select kit
      output$kit_par <- renderUI({
        kit_names <- list.files(paste0(path_pack, "/extdata/kit"))
        kit_names <- gsub(".csv", "", kit_names)
        selectInput(session$ns("kit_name"), label = h4("Kit"), choices = kit_names, selected = kit_names[1], multiple = FALSE, selectize = FALSE)
      })

      observeEvent(input$kit_name, ignoreInit = TRUE, {
        rv_est_sr$dt_setting_est_sr <- dt_setting_est_sr <- create_init_dt_setting_est_sr(input$kit_name, path_pack)

        output$dt_setting_est_sr <- renderRHandsontable({
          rhandsontable(dt_setting_est_sr, rowHeaders = FALSE) %>%
            hot_table(fillHandle = FALSE, rowHeights = 30) %>%
            hot_col(c("Marker"), type = "text", readOnly = TRUE, width = 100) %>%
            hot_col(c("Dye"), type = "text", readOnly = TRUE, width = 50) %>%
            hot_col(c("AT"), type = "numeric", format = "0", width = 50) %>%
            hot_col(c("Model_BSR", "Model_FSR", "Model_DSR", "Model_M2SR"), type = "dropdown", source = c("Allele", "LUS", "Multi-seq", "Uniform", "Uniform (Multiple loci)", "Not consider"), width = 200) %>%
            hot_col(c("X_Multi"), type = "numeric", format = "0", width = 50)
        })
      })

      # Advanced setting


      return(rv_est_sr)
    }
  )
}
