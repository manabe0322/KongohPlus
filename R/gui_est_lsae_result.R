#' est_lsae_result_ui
#'
#' @description The function to create the ui module to create the result tab for estimating the locus-specific amplification efficiencies
est_lsae_result_ui <- function(id){
  ns <- NS(id)

  tabPanel(
    "Result",
    tags$style("
             .edit-modal .modal-dialog {
               width: 90% !important;   /* Expand up to 90% of the screen width */
               max-width: 2000px;       /* Set maximum width */
             }
             .edit-modal .modal-body {
               max-height: 80vh;        /* Limit height to 80% of the viewport */
               overflow-y: auto;        /* Enable vertical scrolling if needed */
             }
           "),
    br(),
    fluidRow(
      column(
        3,
        actionButton(ns("act_saveas_est_lsae_1"), label = "Save as", class = "btn btn-primary btn-lg"),
        br(),
        br(),
        wellPanel(
          h4("Expected LSAE"),
          dataTableOutput(ns("dt_expect_ae")),
          br(),
          uiOutput(ns("button_all_mle"))
        ),
        wellPanel(
          h4("Warnings"),
          htmlOutput(ns("warning_message")),
          br()
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

#' est_lsae_result_server
#'
#' @description The function to create the server module to create the result tab for estimating the locus-specific amplification efficiencies
est_lsae_result_server <- function(id, path_pack, rv_est_lsae){
  moduleServer(
    id,
    function(input, output, session){
      rv_est_lsae_result <- reactiveValues()
      rv_est_lsae_result$warning_message <- NULL
      rv_est_lsae_result$dt_prop_ae_data_use <- NULL
      rv_est_lsae_result$dt_expect_ae <- NULL
      rv_est_lsae_result$opt_list <- NULL
      rv_est_lsae_result$index_converged <- NULL
      rv_est_lsae_result$min_alpha <- NULL
      rv_est_lsae_result$max_alpha <- NULL
      rv_est_lsae_result$df_select_lsae <- NULL

      observe({
        req(rv_est_lsae)
        data_est_lsae <- rv_est_lsae$data_est_lsae

        rv_est_lsae_result$warning_message <- NULL
        rv_est_lsae_result$dt_prop_ae_data_use <- NULL
        rv_est_lsae_result$dt_expect_ae <- NULL
        rv_est_lsae_result$opt_list <- NULL
        rv_est_lsae_result$index_converged <- NULL
        rv_est_lsae_result$min_alpha <- NULL
        rv_est_lsae_result$max_alpha <- NULL
        rv_est_lsae_result$df_select_lsae <- NULL

        if(!is.null(data_est_lsae)){
          rv_est_lsae_result$warning_message <- data_est_lsae$warning_message
          rv_est_lsae_result$dt_prop_ae_data_use <- data_est_lsae$dt_prop_ae_data_use
          rv_est_lsae_result$dt_expect_ae <- data_est_lsae$dt_expect_ae
          rv_est_lsae_result$opt_list <- data_est_lsae$opt_list
          rv_est_lsae_result$index_converged <- data_est_lsae$index_converged
          rv_est_lsae_result$min_alpha <- data_est_lsae$min_alpha
          rv_est_lsae_result$max_alpha <- data_est_lsae$max_alpha
          enable("act_saveas_est_lsae_1")
        }else{
          disable("act_saveas_est_lsae_1")
        }
      })

      observeEvent(input$act_saveas_est_lsae_1, {
        req(rv_est_lsae)
        data_est_lsae <- rv_est_lsae$data_est_lsae

        if(!is.null(data_est_lsae)){
          showModal(modalDialog(
            title = "Save data for locus-specific amplification efficiencies",
            textInput(session$ns("fn_est_lsae"), label = "Enter the file name"),
            footer = tagList(
              disabled(actionButton(session$ns("act_saveas_est_lsae_2"), "Save")),
              modalButton("Close")
            )
          ))
        }

        else{
          showModal(modalDialog(title = "Error", "No data.", easyClose = TRUE, footer = NULL))
        }
      })

      observeEvent(input$fn_est_lsae, {
        fn_est_lsae <- input$fn_est_lsae
        lsae_file_names <- list.files(paste0(path_pack, "/extdata/lsae"))
        if(is.element(paste0(fn_est_lsae, ".RData"), lsae_file_names)){
          showFeedbackDanger(inputId = "fn_est_lsae", text = "A file with the same name already exists.")
          disable("act_saveas_est_lsae_2")
        }else if(nchar(fn_est_lsae) == 0){
          hideFeedback("fn_est_lsae")
          disable("act_saveas_est_lsae_2")
        }else{
          hideFeedback("fn_est_lsae")
          enable("act_saveas_est_lsae_2")
        }
      })

      observeEvent(input$act_saveas_est_lsae_2, {
        req(rv_est_lsae)
        data_est_lsae <- rv_est_lsae$data_est_lsae
        fn_est_lsae <- input$fn_est_lsae

        if(!is.null(data_est_lsae)){
          save(list = c("data_est_lsae"),
               file = paste0(path_pack, "/extdata/lsae/", fn_est_lsae, ".RData"))

          removeModal()
        }

        else{
          removeModal()
          showModal(modalDialog(title = "Error", "No data.", easyClose = TRUE, footer = NULL))
        }
      })

      output$dt_expect_ae <- renderDataTable({
        req(rv_est_lsae_result$dt_expect_ae)
        dt_expect_ae <- rv_est_lsae_result$dt_expect_ae

        datatable(
          dt_expect_ae,
          selection = "none",
          options = list(dom = "t", iDisplayLength = nrow(dt_expect_ae), ordering = FALSE, autoWidth = FALSE),
          rownames = FALSE
        ) %>%
          formatSignif(columns = c("Alpha", "LSAE"), digits = 4) %>%
          formatStyle(names(dt_expect_ae), lineHeight = "80%")
      })

      output$button_all_mle <- renderUI({
        req(rv_est_lsae_result$dt_expect_ae)
        actionButton(session$ns("act_all_mle"), label = "View all MLE attempts")
      })

      observeEvent(input$act_all_mle, {
        showModal(modalDialog(
          title = "All MLE attempts",
          dataTableOutput(session$ns("dt_all_mle")),
          easyClose = TRUE,
          footer = NULL
        ) %>% tagAppendAttributes(class = "edit-modal")
        )
      })

      output$dt_all_mle <- renderDataTable({
        req(rv_est_lsae_result$dt_expect_ae)
        req(rv_est_lsae_result$opt_list)
        req(rv_est_lsae_result$index_converged)
        req(rv_est_lsae_result$min_alpha)
        req(rv_est_lsae_result$max_alpha)

        dt_expect_ae <- rv_est_lsae_result$dt_expect_ae
        opt_list <- rv_est_lsae_result$opt_list
        index_converged <- rv_est_lsae_result$index_converged
        min_alpha <- rv_est_lsae_result$min_alpha
        max_alpha <- rv_est_lsae_result$max_alpha

        n_mk <- length(opt_list[[1]][["par"]])
        count_mle <- length(opt_list)

        mat_all_mle <- matrix(NA, n_mk + 1, count_mle)
        all_scores <- rep(NA, count_mle)

        for(i in 1:count_mle){
          opt_one <- opt_list[[i]]
          all_scores[i] <- score <- opt_one$value
          est_params <- opt_one$par

          alpha <- est_params * (max_alpha - min_alpha) + min_alpha
          lsae <- n_mk * (alpha + 1) / (sum(alpha) + n_mk)

          mat_all_mle[1, i] <- - score
          mat_all_mle[2:(n_mk + 1), i] <- lsae
        }

        dt_all_mle <- as.data.table(mat_all_mle)
        names(dt_all_mle) <- paste0("MLE ", 1:count_mle)

        name_param <- paste0("LSAE ", dt_expect_ae$Marker)
        dt_all_mle[, Parameter := c("Log10 (likelihood)", name_param)]
        setcolorder(dt_all_mle, c("Parameter", setdiff(names(dt_all_mle), "Parameter")))

        name_converged <- paste0("MLE ", index_converged)

        datatable(
          dt_all_mle,
          selection = "none",
          options = list(dom = "t", iDisplayLength = nrow(dt_all_mle), ordering = FALSE, autowidth = TRUE),
          rownames = FALSE
        ) %>%
          formatSignif(columns = paste0("MLE ", 1:count_mle), digits = 6) %>%
          formatStyle(names(dt_all_mle), lineHeight = "80%") %>%
          formatStyle(name_converged, backgroundColor = "lightblue")
      })

      output$warning_message <- renderUI({
        req(rv_est_lsae_result$warning_message)
        HTML(rv_est_lsae_result$warning_message)
      })

      output$ui_plot_lsae <- renderUI({
        req(rv_est_lsae_result$dt_prop_ae_data_use)
        plotOutput(session$ns("plot_lsae"), width = "100%", height = "300px")
      })

      output$plot_lsae <- renderPlot({
        req(rv_est_lsae_result$dt_prop_ae_data_use)
        dt_prop_ae_data_use <- rv_est_lsae_result$dt_prop_ae_data_use
        df_select_lsae <- rv_est_lsae_result$df_select_lsae  # It can be either null or not.

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
        req(rv_est_lsae_result$dt_prop_ae_data_use)
        dt_prop_ae_data_use <- rv_est_lsae_result$dt_prop_ae_data_use

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
        req(rv_est_lsae_result$dt_prop_ae_data_use)
        dt_prop_ae_data_use <- rv_est_lsae_result$dt_prop_ae_data_use
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

          rv_est_lsae_result$df_select_lsae <- df
        }else{
          rv_est_lsae_result$df_select_lsae <- NULL
        }
      })
    }
  )
}
