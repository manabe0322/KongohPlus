#' view_lsae_ui
#'
#' @description The function to create the ui module for viewing the locus-specific amplification efficiencies
#' @param id The ui ID
view_lsae_ui <- function(id){
  ns <- NS(id)

  tabPanel(
    "View the locus-specific amplification efficiencies",
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
    titlePanel("View the locus-specific amplification efficiencies"),
    br(),
    fluidRow(
      column(
        3,
        wellPanel(
          fluidRow(
            column(3, br(), shinyFilesButton(ns("load_lsae_data"), label = "Load", title = "Select a file", multiple = FALSE)),
            column(9, uiOutput(ns("fn_lsae_data")))
          )
        ),
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

#' view_lsae_server
#'
#' @description The function to create the server module for viewing the locus-specific amplification efficiencies
#' @param id The server ID
#' @param path_pack The path of the package KongohPlus
view_lsae_server <- function(id, path_pack){
  moduleServer(
    id,
    function(input, output, session){
      rv_view_lsae <- reactiveValues()
      rv_view_lsae$warning_message <- NULL
      rv_view_lsae$fn_lsae_data <- NULL
      rv_view_lsae$dt_prop_ae_data_use <- NULL
      rv_view_lsae$dt_expect_ae <- NULL
      rv_view_lsae$dt_all_mle <- NULL
      rv_view_lsae$index_mle_best <- NULL
      rv_view_lsae$index_converged <- NULL
      rv_view_lsae$df_select_lsae <- NULL

      shinyFileChoose(input, "load_lsae_data", roots = c(documents = paste0(path_pack, "/extdata/lsae")), session = session, filetypes = c("RData"))

      observeEvent(input$load_lsae_data, {
        req(input$load_lsae_data)
        path_lsae_data <- parseFilePaths(roots = c(documents = paste0(path_pack, "/extdata/lsae")), input$load_lsae_data)$datapath

        if(length(path_lsae_data) > 0){
          rv_view_lsae$fn_lsae_data <- basename(path_lsae_data)

          load(path_lsae_data)

          rv_view_lsae$warning_message <- data_est_lsae$warning_message
          rv_view_lsae$dt_prop_ae_data_use <- data_est_lsae$dt_prop_ae_data_use
          rv_view_lsae$dt_expect_ae <- data_est_lsae$dt_expect_ae
          rv_view_lsae$dt_all_mle <- data_est_lsae$dt_all_mle
          rv_view_lsae$index_mle_best <- data_est_lsae$index_mle_best
          rv_view_lsae$index_converged <- data_est_lsae$index_converged
        }
      })

      output$fn_lsae_data <- renderUI({
        disabled(textInput(session$ns("fn_lsae_data"), label = "", value = rv_view_lsae$fn_lsae_data))
      })

      output$dt_expect_ae <- renderDataTable({
        req(rv_view_lsae$dt_expect_ae)

        dt_expect_ae <- rv_view_lsae$dt_expect_ae

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
        req(rv_view_lsae$dt_all_mle)
        actionButton(session$ns("act_all_mle"), label = "View all MLE attempts")
      })

      observeEvent(input$act_all_mle, {
        showModal(modalDialog(
          title = "All MLE attempts",

          div(
            span(style = "color:#BFEFFF;", "■"), " : Selected maximum likelihood",
            span(style = "color:#ADD8E6;", "■"), " : Close to the maximum likelihood"
          ),

          dataTableOutput(session$ns("dt_all_mle")),
          easyClose = TRUE,
          footer = NULL
        ) %>% tagAppendAttributes(class = "edit-modal")
        )
      })

      output$dt_all_mle <- renderDataTable({
        req(rv_view_lsae$index_mle_best)
        req(rv_view_lsae$index_converged)
        req(rv_view_lsae$dt_all_mle)

        index_mle_best <- rv_view_lsae$index_mle_best
        index_converged <- rv_view_lsae$index_converged
        dt_all_mle <- rv_view_lsae$dt_all_mle

        name_best <- paste0("MLE ", index_mle_best)
        name_converged <- setdiff(paste0("MLE ", index_converged), name_best)

        datatable(
          dt_all_mle,
          selection = "none",
          options = list(dom = "t", iDisplayLength = nrow(dt_all_mle), ordering = FALSE, autowidth = TRUE),
          rownames = FALSE
        ) %>%
          formatSignif(columns = paste0("MLE ", 1:(ncol(dt_all_mle) - 1)), digits = 6) %>%
          formatStyle(names(dt_all_mle), lineHeight = "80%") %>%
          formatStyle(name_best, backgroundColor = "#BFEFFF") %>%
          formatStyle(name_converged, backgroundColor = "#ADD8E6")
      })

      output$warning_message <- renderUI({
        req(rv_view_lsae$warning_message)
        HTML(rv_view_lsae$warning_message)
      })

      output$ui_plot_lsae <- renderUI({
        req(rv_view_lsae$dt_prop_ae_data_use)
        plotOutput(session$ns("plot_lsae"), width = "100%", height = "300px")
      })

      output$plot_lsae <- renderPlot({
        req(rv_view_lsae$dt_prop_ae_data_use)

        dt_prop_ae_data_use <- rv_view_lsae$dt_prop_ae_data_use
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

        df_select_lsae <- rv_view_lsae$df_select_lsae
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
        req(rv_view_lsae$dt_prop_ae_data_use)

        dt_prop_ae_data_use <- rv_view_lsae$dt_prop_ae_data_use
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
        index_select <- input$dt_lsae_data_rows_selected
        if(length(index_select) == 1){
          dt_prop_ae_data_use <- rv_view_lsae$dt_prop_ae_data_use
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

          rv_view_lsae$df_select_lsae <- df
        }else{
          rv_view_lsae$df_select_lsae <- NULL
        }
      })
    }
  )
}
