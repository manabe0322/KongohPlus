#' est_par_result_lsae_ui
#'
#' @description The function to create the ui module to create the result tab for estimating LSAE parameters
est_par_result_lsae_ui <- function(id){
  ns <- NS(id)

  tabPanel(
    "Result (LSAE)",
    br(),
    fluidRow(
      column(
        3,
        actionButton(ns("act_saveas_est_par_1"), label = "Save as", class = "btn btn-primary btn-lg"),
        br(),
        wellPanel(
          h4("Expected LSAE"),
          dataTableOutput(ns("dt_expect_ae"))
        )
      ),
      column(
        9,
        wellPanel(
          h4("Experimental LSAE data"),
          uiOutput(ns("ui_plot_lsae")),
          br(),
          dataTableOutput(ns("dt_lsae_data"))
        )
      )
    )
  )
}

#' est_par_result_lsae_server
#'
#' @description The function to create the server module to create the result tab for estimating LSAE parameters
est_par_result_lsae_server <- function(id, path_pack, rv_est_par){
  moduleServer(
    id,
    function(input, output, session){
      rv_est_par_result_lsae <- reactiveValues()
      rv_est_par_result_lsae$dt_prop_ae_data_use <- NULL
      rv_est_par_result_lsae$dt_expect_ae <- NULL
      rv_est_par_result_lsae$df_select_lsae <- NULL

      observe({
        req(rv_est_par)
        data_est_par <- rv_est_par$data_est_par

        if(!is.null(data_est_par)){
          rv_est_par_result_lsae$dt_prop_ae_data_use <- data_est_par$dt_prop_ae_data_use
          rv_est_par_result_lsae$dt_expect_ae <- data_est_par$dt_expect_ae
          rv_est_par_result_lsae$df_select_lsae <- NULL
          enable("act_saveas_est_par_1")
        }else{
          rv_est_par_result_lsae$dt_prop_ae_data_use <- NULL
          rv_est_par_result_lsae$dt_expect_ae <- NULL
          rv_est_par_result_lsae$df_select_lsae <- NULL
          disable("act_saveas_est_par_1")
        }
      })

      observeEvent(input$act_saveas_est_par_1, {
        req(rv_est_sr)
        data_est_par <- rv_est_par$data_est_par

        if(!is.null(data_est_par)){
          showModal(modalDialog(
            title = "Save parameter data",
            textInput(session$ns("fn_est_par"), label = "Enter the file name"),
            footer = tagList(
              disabled(actionButton(session$ns("act_saveas_est_par_2"), "Save")),
              modalButton("Close")
            )
          ))
        }

        else{
          showModal(modalDialog(title = "Error", "No parameter data.", easyClose = TRUE, footer = NULL))
        }
      })

      observeEvent(input$fn_est_par, {
        fn_est_par <- input$fn_est_par
        par_file_names <- list.files(paste0(path_pack, "/extdata/param"))
        if(is.element(paste0(fn_est_par, ".RData"), par_file_names)){
          showFeedbackDanger(inputId = "fn_est_par", text = "A file with the same name already exists.")
          disable("act_saveas_est_par_2")
        }else if(nchar(fn_est_par) == 0){
          hideFeedback("fn_est_par")
          disable("act_saveas_est_par_2")
        }else{
          hideFeedback("fn_est_par")
          enable("act_saveas_est_par_2")
        }
      })

      observeEvent(input$act_saveas_est_par_2, {
        req(rv_est_par)
        data_est_par <- rv_est_par$data_est_par
        fn_est_par <- input$fn_est_par

        if(!is.null(data_est_par)){
          save(list = c("data_est_par"),
               file = paste0(path_pack, "/extdata/param/", fn_est_par, ".RData"))

          removeModal()
        }

        else{
          removeModal()
          showModal(modalDialog(title = "Error", "No parameter data.", easyClose = TRUE, footer = NULL))
        }
      })

      output$dt_expect_ae <- renderDataTable({
        req(rv_est_par_result_lsae$dt_expect_ae)
        dt_expect_ae <- rv_est_par_result_lsae$dt_expect_ae

        datatable(
          dt_expect_ae,
          selection = "none",
          options = list(dom = "t", iDisplayLength = nrow(dt_expect_ae), ordering = FALSE, autoWidth = FALSE),
          rownames = FALSE
        ) %>%
          formatSignif(columns = c("Alpha", "LSAE"), digits = 4) %>%
          formatStyle(names(dt_expect_ae), lineHeight = "80%")
      })

      output$ui_plot_lsae <- renderUI({
        req(rv_est_par_result_lsae$dt_prop_ae_data_use)
        plotOutput(session$ns("plot_lsae"), width = "100%", height = "300px")
      })

      output$plot_lsae <- renderPlot({
        req(rv_est_par_result_lsae$dt_prop_ae_data_use)
        dt_prop_ae_data_use <- rv_est_par_result_lsae$dt_prop_ae_data_use
        df_select_lsae <- rv_est_par_result_lsae$df_select_lsae  # It can be either null or not.

        n_sample <- nrow(dt_prop_ae_data_use)
        n_mk <- ncol(dt_prop_ae_data_use) - 1
        markers <- colnames(dt_prop_ae_data_use)[-1]

        df <- data.frame(
          Marker = as.vector(sapply(markers, rep, n_sample)),
          LSAE = unlist(dt_prop_ae_data_use[, -c("SampleName"), with = FALSE])
        )
        df$Marker <- factor(df$Marker, levels = markers)

        base_plot <- ggplot(df, aes(x = Marker, y = LSAE)) +
          geom_boxplot(outlier.colour = "#c0c0c0")

        if(!is.null(df_select_lsae)){
          base_plot <- base_plot +
            geom_point(
              data = df_select_lsae,
              aes(x = Marker, y = LSAE),
              color = "#e9bc00"
            )
        }

        base_plot +
          theme(
            plot.title = element_text(size = 14),
            legend.text = element_text(size = 12),
            axis.title.x = element_text(size = 12),
            axis.title.y = element_text(size = 12),
            axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(size = 12)
          )
      })

      output$dt_lsae_data <- renderDataTable({
        req(rv_est_par_result_lsae$dt_prop_ae_data_use)
        dt_prop_ae_data_use <- rv_est_par_result_lsae$dt_prop_ae_data_use

        markers <- colnames(dt_prop_ae_data_use)[-1]

        datatable(
          dt_prop_ae_data_use,
          selection = list(mode = "single", target = "row"),
          options = list(
            iDisplayLength = 10,
            autoWidth = FALSE,
            scrollX = TRUE,
            columnDefs = list(
              list(
                targets = 0,
                width = '300px' # Do not work because of too many columns
              )
            )
          ),
          rownames = FALSE
        ) %>%
          formatSignif(columns = markers, digits = 4) %>%
          formatStyle(names(dt_prop_ae_data_use), lineHeight = "80%")
      })

      observeEvent(input$dt_lsae_data_rows_selected, {
        req(rv_est_par_result_lsae$dt_prop_ae_data_use)
        dt_prop_ae_data_use <- rv_est_par_result_lsae$dt_prop_ae_data_use
        index_select <- input$dt_lsae_data_rows_selected

        if(length(index_select) == 1){
          markers <- colnames(dt_prop_ae_data_use)[-1]

          n_mk <- length(markers)
          selected_lsae <- rep(NA, n_mk)
          for(i in 1:n_mk){
            selected_lsae[i] <- dt_prop_ae_data_use[[markers[i]]][index_select]
          }

          df <- data.frame(
            Marker = markers,
            LSAE = selected_lsae
          )
          df$Marker <- factor(df$Marker, levels = markers)

          rv_est_par_result_lsae$df_select_lsae <- df
        }else{
          rv_est_par_result_lsae$df_select_lsae <- NULL
        }
      })
    }
  )
}
