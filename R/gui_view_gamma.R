#' view_gamma_ui
#'
#' @description The function to create the ui module for viewing peak height parameters
view_gamma_ui <- function(id){
  ns <- NS(id)

  tabPanel(
    "View peak height parameters",
    titlePanel("View peak height parameters"),
    br(),
    fluidRow(
      column(
        3,
        wellPanel(
          fluidRow(
            column(3, br(), shinyFilesButton(ns("load_param_data"), label = "Load", title = "Select a file", multiple = FALSE)),
            column(9, uiOutput(ns("fn_param_data")))
          )
        ),
        wellPanel(
          h4("Estimated parameters"),
          dataTableOutput(ns("dt_est_gamma"))
        )
      ),
      column(
        9,
        wellPanel(
          h4("Experimental peak height data (Allele)"),
          fluidRow(
            column(
              6,
              uiOutput(ns("ui_plot_gamma_allele"))
            ),
            column(
              6
            )
          )

        ),
        wellPanel(
          h4("Experimental peak height data (Stutter)"),
          fluidRow(
            column(
              6,
              uiOutput(ns("ui_plot_gamma_stutter"))
            ),
            column(
              6
            )
          )
        ),
        wellPanel(
          h4("Experimental peak height data (Allele + Stutter)"),
          fluidRow(
            column(
              6,
              uiOutput(ns("ui_plot_gamma_mix"))
            ),
            column(
              6
            )
          )
        )
      )
    )
  )
}

#' view_gamma_server
#'
#' @description The function to create the server module for viewing peak height parameters
view_gamma_server <- function(id, path_pack){
  moduleServer(
    id,
    function(input, output, session){
      rv_view_gamma <- reactiveValues()
      rv_view_gamma$fn_param_data <- NULL
      rv_view_gamma$dt_est_gamma <- NULL
      rv_view_gamma$dt_gamma_data <- NULL

      shinyFileChoose(input, "load_param_data", roots = c(documents = paste0(path_pack, "/extdata/param")), session = session, filetypes = c("RData"))

      observeEvent(input$load_param_data, {
        req(input$load_param_data)
        path_param_data <- parseFilePaths(roots = c(documents = paste0(path_pack, "/extdata/param")), input$load_param_data)$datapath

        if(length(path_param_data) > 0){
          rv_view_gamma$fn_param_data <- basename(path_param_data)

          load(path_param_data)

          rv_view_gamma$dt_gamma_data <- fread("C:/Users/manab/Desktop/dt_gamma_data.csv")

          rv_view_gamma$dt_est_gamma <- data.table(Parameter = c("CV for allele", "CV for stutter"),
                                                   Value = c(cv_allele, cv_stutter))
        }
      })

      output$fn_param_data <- renderUI({
        disabled(textInput(session$ns("fn_param_data"), label = "", value = rv_view_gamma$fn_param_data))
      })

      output$dt_est_gamma <- renderDataTable({
        dt_est_gamma <- rv_view_gamma$dt_est_gamma

        if(!is.null(dt_est_gamma)){
          datatable(
            dt_est_gamma,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_est_gamma), ordering = FALSE, autoWidth = FALSE, fixedHeader = TRUE),
            rownames = FALSE
          ) %>%
            formatSignif(columns = c("Value"), digits = 4) %>%
            formatStyle(names(dt_est_gamma), lineHeight = "80%")
        }
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
          index_target <- which(dt_gamma_data[["Peak type"]] == "allele")
          n_data <- length(index_target)
          ephs <- dt_gamma_data[["Expected peak height"]][index_target]
          ophs <- dt_gamma_data[["Observed peak height"]][index_target]
          ats <- dt_gamma_data[["AT"]][index_target]
          quantile95 <- dt_gamma_data[["95% quantile"]][index_target]
          quantile5 <- dt_gamma_data[["5% quantile"]][index_target]

          index_sort <- order(ephs)

          df_data <- data.frame(Expected = log10(ephs[index_sort]), Observed = log10(ophs[index_sort]))
          df_data <- df_data[is.finite(df_data$Observed), ]
          df_95 <- data.frame(Expected = log10(ephs[index_sort]), Observed = log10(quantile95[index_sort]))
          df_5 <- data.frame(Expected = log10(ephs[index_sort]), Observed = log10(quantile5[index_sort]))

          ggplot() +
            geom_point(data = df_data, aes(Expected, Observed, color = "Data")) +
            geom_line(data = df_95, aes(Expected, Observed, color = "5% and 95% quantile")) +
            geom_line(data = df_5, aes(Expected, Observed, color = "5% and 95% quantile")) +
            xlim(min(df_data$Expected), max(df_data$Expected)) +
            ylim(min(df_data$Observed), max(df_data$Observed)) +
            labs(color = NULL) +
            scale_color_manual(values = c("Data" = "#c0c0c0", "5% and 95% quantile" = "#00BFC4"))
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
          index_target <- which(dt_gamma_data[["Peak type"]] == "stutter")
          n_data <- length(index_target)
          ephs <- dt_gamma_data[["Expected peak height"]][index_target]
          ophs <- dt_gamma_data[["Observed peak height"]][index_target]
          ats <- dt_gamma_data[["AT"]][index_target]
          quantile95 <- dt_gamma_data[["95% quantile"]][index_target]
          quantile5 <- dt_gamma_data[["5% quantile"]][index_target]

          index_sort <- order(ephs)

          df_data <- data.frame(Expected = log10(ephs[index_sort]), Observed = log10(ophs[index_sort]))
          df_data <- df_data[is.finite(df_data$Observed), ]
          df_95 <- data.frame(Expected = log10(ephs[index_sort]), Observed = log10(quantile95[index_sort]))
          df_5 <- data.frame(Expected = log10(ephs[index_sort]), Observed = log10(quantile5[index_sort]))

          ggplot() +
            geom_point(data = df_data, aes(Expected, Observed, color = "Data")) +
            geom_line(data = df_95, aes(Expected, Observed, color = "5% and 95% quantile")) +
            geom_line(data = df_5, aes(Expected, Observed, color = "5% and 95% quantile")) +
            xlim(min(df_data$Expected), max(df_data$Expected)) +
            ylim(min(df_data$Observed), max(df_data$Observed)) +
            labs(color = NULL) +
            scale_color_manual(values = c("Data" = "#c0c0c0", "5% and 95% quantile" = "#00BFC4"))
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
          index_target <- which(dt_gamma_data[["Peak type"]] == "mix")
          n_data <- length(index_target)
          ephs <- dt_gamma_data[["Expected peak height"]][index_target]
          ophs <- dt_gamma_data[["Observed peak height"]][index_target]
          ats <- dt_gamma_data[["AT"]][index_target]
          quantile95 <- dt_gamma_data[["95% quantile"]][index_target]
          quantile5 <- dt_gamma_data[["5% quantile"]][index_target]

          index_sort <- order(ephs)

          df_data <- data.frame(Expected = log10(ephs[index_sort]), Observed = log10(ophs[index_sort]))
          df_data <- df_data[is.finite(df_data$Observed), ]
          df_95 <- data.frame(Expected = log10(ephs[index_sort]), Observed = log10(quantile95[index_sort]))
          df_5 <- data.frame(Expected = log10(ephs[index_sort]), Observed = log10(quantile5[index_sort]))

          ggplot() +
            geom_point(data = df_data, aes(Expected, Observed, color = "Data")) +
            geom_line(data = df_95, aes(Expected, Observed, color = "5% and 95% quantile")) +
            geom_line(data = df_5, aes(Expected, Observed, color = "5% and 95% quantile")) +
            xlim(min(df_data$Expected), max(df_data$Expected)) +
            ylim(min(df_data$Observed), max(df_data$Observed)) +
            labs(color = NULL) +
            scale_color_manual(values = c("Data" = "#c0c0c0", "5% and 95% quantile" = "#00BFC4"))
        }
      })
    }
  )
}
