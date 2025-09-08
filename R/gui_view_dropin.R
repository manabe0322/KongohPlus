#' view_dropin_ui
#'
#' @description The function to create the ui module for viewing drop-in parameters
view_dropin_ui <- function(id){
  ns <- NS(id)

  tabPanel(
    "View drop-in parameters",
    titlePanel("View drop-in parameters"),
    br(),
    fluidRow(
      column(
        3,
        wellPanel(
          fluidRow(
            column(3, br(), shinyFilesButton(ns("load_dropin_data"), label = "Load", title = "Select a file", multiple = FALSE)),
            column(9, uiOutput(ns("fn_dropin_data")))
          )
        )
      ),
      column(
        4,
        wellPanel(
          h4("Estimated parameters"),
          dataTableOutput(ns("dt_est_dropin"))
        )
      ),
      column(
        5,
        wellPanel(
          h4("Experimental drop-in data"),
          uiOutput(ns("ui_plot_dropin")),
          br(),
          dataTableOutput(ns("dt_dropin_data"))
        )
      )
    )
  )
}

#' view_dropin_server
#'
#' @description The function to create the server module for viewing drop-in parameters
view_dropin_server <- function(id, path_pack){
  moduleServer(
    id,
    function(input, output, session){
      rv_view_dropin <- reactiveValues()
      rv_view_dropin$fn_dropin_data <- NULL
      rv_view_dropin$dt_est_dropin <- NULL
      rv_view_dropin$dt_dropin_data <- NULL

      shinyFileChoose(input, "load_dropin_data", roots = c(documents = paste0(path_pack, "/extdata/dropin")), session = session, filetypes = c("RData"))

      observeEvent(input$load_dropin_data, {
        req(input$load_dropin_data)
        path_dropin_data <- parseFilePaths(roots = c(documents = paste0(path_pack, "/extdata/dropin")), input$load_dropin_data)$datapath

        if(length(path_dropin_data) > 0){
          rv_view_dropin$fn_dropin_data <- basename(path_dropin_data)

          load(path_dropin_data)

          rv_view_dropin$dt_est_dropin <- data_est_dropin$dt_est_dropin
          rv_view_dropin$dt_dropin_data <- data_est_dropin$dt_dropin_data
        }
      })

      output$fn_dropin_data <- renderUI({
        disabled(textInput(session$ns("fn_dropin_data"), label = "", value = rv_view_dropin$fn_dropin_data))
      })

      output$dt_est_dropin <- renderDataTable({
        dt_est_dropin <- rv_view_dropin$dt_est_dropin

        if(!is.null(dt_est_dropin)){
          datatable(
            dt_est_dropin,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_est_dropin), ordering = FALSE, autowidth = FALSE),
            rownames = FALSE
          )
        }
      })

      output$ui_plot_dropin <- renderUI({
        dt_est_dropin <- rv_view_dropin$dt_est_dropin
        dt_dropin_data <- rv_view_dropin$dt_dropin_data

        if(!is.null(dt_est_dropin) && !is.null(dt_dropin_data) && nrow(dt_dropin_data) != 0){
          plotOutput(session$ns("plot_dropin"), width = "100%", height = "300px")
        }
      })

      output$plot_dropin <- renderPlot({
        dt_est_dropin <- rv_view_dropin$dt_est_dropin
        dt_dropin_data <- rv_view_dropin$dt_dropin_data

        if(!is.null(dt_est_dropin) && !is.null(dt_dropin_data) && nrow(dt_dropin_data) != 0){
          heights <- dt_dropin_data$Height
          shape_di <- dt_est_dropin$Value[dt_est_dropin$Parameter == "Shape parameter"]
          rate_di <- dt_est_dropin$Value[dt_est_dropin$Parameter == "Rate parameter"]

          tb_ph <- table(heights)
          x_val <- 0:round(max(heights) * 1.5, 0)
          y_val <- rep(0, length(x_val))
          y_val[match(as.numeric(names(tb_ph)), x_val)] <- tb_ph

          ggplot() +
            geom_line(data = data.frame(X = x_val, Y = y_val, P = y_val * 0.01), aes(X, P, color = "Observed"), linewidth = 1) +
            stat_function(fun = dgamma, args = list(shape = shape_di, rate = rate_di), aes(color = "Estimated"), linewidth = 1) +
            xlim(0, max(x_val)) +
            labs(color = NULL, x = "Peak height (RFU)") +
            scale_y_continuous(
              name = "Probability density",
              sec.axis = sec_axis(~ . * 100, name = "Count")
            ) +
            scale_color_manual(values = c("Observed" = "#4169e1", "Estimated" = "#c71585")) +
            theme(
              plot.title = element_text(size = 14),
              legend.text = element_text(size = 12),
              axis.title.x = element_text(size = 12),
              axis.title.y = element_text(size = 12),
              axis.text.x = element_text(size = 12),
              axis.text.y = element_text(size = 12)
            )
        }
      })

      output$dt_dropin_data <- renderDataTable({
        dt_dropin_data <- rv_view_dropin$dt_dropin_data

        if(!is.null(dt_dropin_data)){
          datatable(
            dt_dropin_data,
            selection = "none",
            options = list(
              iDisplayLength = 10,
              autoWidth = FALSE
            ),
            rownames = FALSE
          ) %>%
            formatStyle(names(dt_dropin_data), lineHeight = "80%")
        }
      })
    }
  )
}
