#' view_gamma_ui
#'
#' @description The function to create the ui module for viewing the gamma scale parameter
#' @param id The ui ID
view_gamma_ui <- function(id){
  ns <- NS(id)

  tabPanel(
    "View the gamma scale parameter",
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
    titlePanel("View the gamma scale parameter"),
    br(),
    fluidRow(
      column(
        3,
        wellPanel(
          fluidRow(
            column(3, br(), shinyFilesButton(ns("load_gamma_data"), label = "Load", title = "Select a file", multiple = FALSE)),
            column(9, uiOutput(ns("fn_gamma_data")))
          )
        ),
        wellPanel(
          h4("Estimated parameters"),
          dataTableOutput(ns("dt_est_gamma")),
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
          h4("Experimental peak height data"),
          fluidRow(
            column(
              4,
              uiOutput(ns("ui_plot_gamma_allele"))
            ),
            column(
              4,
              uiOutput(ns("ui_plot_gamma_stutter"))
            ),
            column(
              4,
              uiOutput(ns("ui_plot_gamma_mix"))
            )
          ),
          br(),
          dataTableOutput(ns("dt_gamma_data"))
        )
      )
    )
  )
}

#' view_gamma_server
#'
#' @description The function to create the server module for viewing the gamma scale parameter
#' @param id The server ID
#' @param path_pack The path of the package KongohPlus
view_gamma_server <- function(id, path_pack){
  moduleServer(
    id,
    function(input, output, session){
      rv_view_gamma <- reactiveValues()
      rv_view_gamma$fn_gamma_data <- NULL
      rv_view_gamma$warning_message <- NULL
      rv_view_gamma$dt_gamma_data <- NULL
      rv_view_gamma$dt_est_gamma <- NULL
      rv_view_gamma$dt_all_mle <- NULL
      rv_view_gamma$index_mle_best <- NULL
      rv_view_gamma$index_converged <- NULL

      rv_view_gamma$df_select_gamma_allele <- NULL
      rv_view_gamma$df_select_gamma_stutter <- NULL
      rv_view_gamma$df_select_gamma_mix <- NULL

      shinyFileChoose(input, "load_gamma_data", roots = c(documents = paste0(path_pack, "/extdata/gamma")), session = session, filetypes = c("RData"))

      observeEvent(input$load_gamma_data, {
        req(input$load_gamma_data)
        path_gamma_data <- parseFilePaths(roots = c(documents = paste0(path_pack, "/extdata/gamma")), input$load_gamma_data)$datapath

        if(length(path_gamma_data) > 0){
          rv_view_gamma$fn_gamma_data <- basename(path_gamma_data)

          load(path_gamma_data)

          rv_view_gamma$warning_message <- data_est_gamma$warning_message
          rv_view_gamma$dt_gamma_data <- data_est_gamma$dt_gamma_data
          rv_view_gamma$dt_est_gamma <- data_est_gamma$dt_est_gamma
          rv_view_gamma$dt_all_mle <- data_est_gamma$dt_all_mle
          rv_view_gamma$index_mle_best <- data_est_gamma$index_mle_best
          rv_view_gamma$index_converged <- data_est_gamma$index_converged
        }
      })

      output$fn_gamma_data <- renderUI({
        disabled(textInput(session$ns("fn_gamma_data"), label = "", value = rv_view_gamma$fn_gamma_data))
      })

      output$dt_est_gamma <- renderDataTable({
        req(rv_view_gamma$dt_est_gamma)

        dt_est_gamma <- rv_view_gamma$dt_est_gamma

        datatable(
          dt_est_gamma,
          selection = "none",
          options = list(dom = "t", iDisplayLength = nrow(dt_est_gamma), ordering = FALSE, autoWidth = FALSE, fixedHeader = TRUE),
          rownames = FALSE
        ) %>%
          formatSignif(columns = c("Value"), digits = 4) %>%
          formatStyle(names(dt_est_gamma), lineHeight = "80%")
      })

      output$button_all_mle <- renderUI({
        req(rv_view_gamma$dt_all_mle)
        actionButton(session$ns("act_all_mle"), label = "View all MLE attempts")
      })

      observeEvent(input$act_all_mle, {
        showModal(modalDialog(
          title = "All MLE attempts",

          div(
            span(style = "color:#BFEFFF;", "\u25A0"), " : Selected maximum likelihood",
            span(style = "color:#ADD8E6;", "\u25A0"), " : Close to the maximum likelihood"
          ),

          dataTableOutput(session$ns("dt_all_mle")),
          easyClose = TRUE,
          footer = NULL
        ) %>% tagAppendAttributes(class = "edit-modal")
        )
      })

      output$dt_all_mle <- renderDataTable({
        req(rv_view_gamma$index_mle_best)
        req(rv_view_gamma$index_converged)
        req(rv_view_gamma$dt_all_mle)

        index_mle_best <- rv_view_gamma$index_mle_best
        index_converged <- rv_view_gamma$index_converged
        dt_all_mle <- rv_view_gamma$dt_all_mle

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
        req(rv_view_gamma$warning_message)
        HTML(rv_view_gamma$warning_message)
      })

      output$ui_plot_gamma_allele <- renderUI({
        dt_gamma_data <- rv_view_gamma$dt_gamma_data

        if(!is.null(dt_gamma_data)){
          plotOutput(session$ns("plot_gamma_allele"), width = "100%", height = "300px")
        }
      })

      output$plot_gamma_allele <- renderPlot({
        dt_gamma_data <- rv_view_gamma$dt_gamma_data

        if(!is.null(dt_gamma_data)){
          index_target <- which(dt_gamma_data[["Peak type"]] == "Allele")
          n_data <- length(index_target)
          ephs <- dt_gamma_data[["Expected peak height"]][index_target]
          ophs <- dt_gamma_data[["Observed peak height"]][index_target]
          ats <- dt_gamma_data[["AT"]][index_target]
          quantile95 <- dt_gamma_data[["95% quantile"]][index_target]
          quantile5 <- dt_gamma_data[["5% quantile"]][index_target]

          index_sort <- order(ephs)

          df_data <- data.frame(Expected = ephs[index_sort], Observed = ophs[index_sort])
          df_95 <- data.frame(Expected = ephs[index_sort], Observed = quantile95[index_sort])
          df_5 <- data.frame(Expected = ephs[index_sort], Observed = quantile5[index_sort])
          df_select_gamma_allele <- rv_view_gamma$df_select_gamma_allele

          color_values <- c("Data" = "#c0c0c0", "5% and 95% quantile" = "#00BFC4")

          base_plot <- ggplot() +
            geom_point(data = df_data, aes(Expected, Observed, color = "Data"), na.rm = TRUE) +
            geom_line(data = df_95, aes(Expected, Observed, color = "5% and 95% quantile"), na.rm = TRUE) +
            geom_line(data = df_5, aes(Expected, Observed, color = "5% and 95% quantile"), na.rm = TRUE)

          if(!is.null(df_select_gamma_allele)){
            base_plot <- base_plot +
              geom_point(data = df_select_gamma_allele, aes(Expected, Observed, color = "Selected data"), na.rm = TRUE)

            color_values <- c(color_values, "Selected data" = "#e9bc00")
          }

          base_plot +
            scale_x_log10(limits = c(min(df_data$Expected[df_data$Expected != 0]), max(df_data$Expected[df_data$Expected != 0]))) +
            scale_y_log10(limits = c(min(df_data$Observed[df_data$Observed != 0]), max(df_data$Observed[df_data$Observed != 0]))) +
            labs(color = NULL, x = "Expected peak height", y = "Observed peak height", title = "Allele") +
            scale_color_manual(values = color_values) +
            theme(
              plot.title = element_text(size = 14),
              legend.text = element_text(size = 12),
              legend.position = "top",
              axis.title.x = element_text(size = 12),
              axis.title.y = element_text(size = 12),
              axis.text.x = element_text(size = 12),
              axis.text.y = element_text(size = 12)
            )
        }
      })

      output$ui_plot_gamma_stutter <- renderUI({
        dt_gamma_data <- rv_view_gamma$dt_gamma_data

        if(!is.null(dt_gamma_data)){
          plotOutput(session$ns("plot_gamma_stutter"), width = "100%", height = "300px")
        }
      })

      output$plot_gamma_stutter <- renderPlot({
        dt_gamma_data <- rv_view_gamma$dt_gamma_data

        if(!is.null(dt_gamma_data)){
          index_target <- which(dt_gamma_data[["Peak type"]] == "Stutter")
          n_data <- length(index_target)
          ephs <- dt_gamma_data[["Expected peak height"]][index_target]
          ophs <- dt_gamma_data[["Observed peak height"]][index_target]
          ats <- dt_gamma_data[["AT"]][index_target]
          quantile95 <- dt_gamma_data[["95% quantile"]][index_target]
          quantile5 <- dt_gamma_data[["5% quantile"]][index_target]

          index_sort <- order(ephs)

          df_data <- data.frame(Expected = ephs[index_sort], Observed = ophs[index_sort])
          df_95 <- data.frame(Expected = ephs[index_sort], Observed = quantile95[index_sort])
          df_5 <- data.frame(Expected = ephs[index_sort], Observed = quantile5[index_sort])
          df_select_gamma_stutter <- rv_view_gamma$df_select_gamma_stutter

          color_values <- c("Data" = "#c0c0c0", "5% and 95% quantile" = "#00BFC4")

          base_plot <- ggplot() +
            geom_point(data = df_data, aes(Expected, Observed, color = "Data"), na.rm = TRUE) +
            geom_line(data = df_95, aes(Expected, Observed, color = "5% and 95% quantile"), na.rm = TRUE) +
            geom_line(data = df_5, aes(Expected, Observed, color = "5% and 95% quantile"), na.rm = TRUE)

          if(!is.null(df_select_gamma_stutter)){
            base_plot <- base_plot +
              geom_point(data = df_select_gamma_stutter, aes(Expected, Observed, color = "Selected data"), na.rm = TRUE)

            color_values <- c(color_values, "Selected data" = "#e9bc00")
          }

          base_plot +
            scale_x_log10(limits = c(min(df_data$Expected[df_data$Expected != 0]), max(df_data$Expected[df_data$Expected != 0]))) +
            scale_y_log10(limits = c(min(df_data$Observed[df_data$Observed != 0]), max(df_data$Observed[df_data$Observed != 0]))) +
            labs(color = NULL, x = "Expected peak height", y = "Observed peak height", title = "Stutter") +
            scale_color_manual(values = color_values) +
            theme(
              plot.title = element_text(size = 14),
              legend.text = element_text(size = 12),
              legend.position = "top",
              axis.title.x = element_text(size = 12),
              axis.title.y = element_text(size = 12),
              axis.text.x = element_text(size = 12),
              axis.text.y = element_text(size = 12)
            )
        }
      })

      output$ui_plot_gamma_mix <- renderUI({
        dt_gamma_data <- rv_view_gamma$dt_gamma_data

        if(!is.null(dt_gamma_data)){
          plotOutput(session$ns("plot_gamma_mix"), width = "100%", height = "300px")
        }
      })

      output$plot_gamma_mix <- renderPlot({
        dt_gamma_data <- rv_view_gamma$dt_gamma_data

        if(!is.null(dt_gamma_data)){
          index_target <- which(dt_gamma_data[["Peak type"]] == "Allele + Stutter")
          n_data <- length(index_target)
          ephs <- dt_gamma_data[["Expected peak height"]][index_target]
          ophs <- dt_gamma_data[["Observed peak height"]][index_target]
          ats <- dt_gamma_data[["AT"]][index_target]
          quantile95 <- dt_gamma_data[["95% quantile"]][index_target]
          quantile5 <- dt_gamma_data[["5% quantile"]][index_target]

          index_sort <- order(ephs)

          df_data <- data.frame(Expected = ephs[index_sort], Observed = ophs[index_sort])
          df_95 <- data.frame(Expected = ephs[index_sort], Observed = quantile95[index_sort])
          df_5 <- data.frame(Expected = ephs[index_sort], Observed = quantile5[index_sort])
          df_select_gamma_mix <- rv_view_gamma$df_select_gamma_mix

          color_values <- c("Data" = "#c0c0c0", "5% and 95% quantile" = "#00BFC4")

          base_plot <- ggplot() +
            geom_point(data = df_data, aes(Expected, Observed, color = "Data"), na.rm = TRUE) +
            geom_line(data = df_95, aes(Expected, Observed, color = "5% and 95% quantile"), na.rm = TRUE) +
            geom_line(data = df_5, aes(Expected, Observed, color = "5% and 95% quantile"), na.rm = TRUE)

          if(!is.null(df_select_gamma_mix)){
            base_plot <- base_plot +
              geom_point(data = df_select_gamma_mix, aes(Expected, Observed, color = "Selected data"), na.rm = TRUE)

            color_values <- c(color_values, "Selected data" = "#e9bc00")
          }

          base_plot +
            scale_x_log10(limits = c(min(df_data$Expected[df_data$Expected != 0]), max(df_data$Expected[df_data$Expected != 0]))) +
            scale_y_log10(limits = c(min(df_data$Observed[df_data$Observed != 0]), max(df_data$Observed[df_data$Observed != 0]))) +
            labs(color = NULL, x = "Expected peak height", y = "Observed peak height", title = "Allele + Stutter") +
            scale_color_manual(values = color_values) +
            theme(
              plot.title = element_text(size = 14),
              legend.text = element_text(size = 12),
              legend.position = "top",
              axis.title.x = element_text(size = 12),
              axis.title.y = element_text(size = 12),
              axis.text.x = element_text(size = 12),
              axis.text.y = element_text(size = 12)
            )
        }
      })

      output$dt_gamma_data <- renderDataTable({
        dt_gamma_data <- rv_view_gamma$dt_gamma_data

        if(!is.null(dt_gamma_data)){
          dt_gamma_data_display <- data.table("SampleName" = factor(dt_gamma_data[["SampleName"]], levels = unique(dt_gamma_data[["SampleName"]]), labels = unique(dt_gamma_data[["SampleName"]])),
                                              "Marker" = factor(dt_gamma_data[["Marker"]], levels = unique(dt_gamma_data[["Marker"]]), labels = unique(dt_gamma_data[["Marker"]])),
                                              "Peak" = dt_gamma_data[["Peak"]],
                                              "Observed peak height" = dt_gamma_data[["Observed peak height"]],
                                              "Expected peak height" = dt_gamma_data[["Expected peak height"]],
                                              "Peak type" = factor(dt_gamma_data[["Peak type"]], levels = unique(dt_gamma_data[["Peak type"]]), labels = unique(dt_gamma_data[["Peak type"]]))
                                              )

          datatable(
            dt_gamma_data_display,
            filter = "top",
            selection = list(mode = "single", target = "row"),
            options = list(
              iDisplayLength = 10,
              autoWidth = FALSE,
              columnDefs = list(
                list(
                  targets = 0,
                  width = '250px'
                ),
                list(
                  targets = c(1, 5),
                  width = '125px'
                ),
                list(
                  targets = c(2, 3, 4),
                  searchable = FALSE
                )
              )
            ),
            rownames = FALSE
          ) %>%
            formatRound(columns = c("Expected peak height"), digits = 2) %>%
            formatStyle(names(dt_gamma_data_display), lineHeight = "80%")
        }
      })

      observeEvent(input$dt_gamma_data_rows_selected, {
        index_select <- input$dt_gamma_data_rows_selected
        if(length(index_select) == 1){
          dt_gamma_data <- rv_view_gamma$dt_gamma_data

          if(dt_gamma_data[["Peak type"]][index_select] == "Allele"){
            rv_view_gamma$df_select_gamma_allele <- data.frame(Expected = dt_gamma_data[["Expected peak height"]][index_select], Observed = dt_gamma_data[["Observed peak height"]][index_select])
            rv_view_gamma$df_select_gamma_stutter <- NULL
            rv_view_gamma$df_select_gamma_mix <- NULL
          }else if(dt_gamma_data[["Peak type"]][index_select] == "Stutter"){
            rv_view_gamma$df_select_gamma_allele <- NULL
            rv_view_gamma$df_select_gamma_stutter <- data.frame(Expected = dt_gamma_data[["Expected peak height"]][index_select], Observed = dt_gamma_data[["Observed peak height"]][index_select])
            rv_view_gamma$df_select_gamma_mix <- NULL
          }else if(dt_gamma_data[["Peak type"]][index_select] == "Allele + Stutter"){
            rv_view_gamma$df_select_gamma_allele <- NULL
            rv_view_gamma$df_select_gamma_stutter <- NULL
            rv_view_gamma$df_select_gamma_mix <- data.frame(Expected = dt_gamma_data[["Expected peak height"]][index_select], Observed = dt_gamma_data[["Observed peak height"]][index_select])
          }
        }
      })
    }
  )
}
