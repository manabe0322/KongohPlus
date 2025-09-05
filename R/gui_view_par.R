#' view_par_ui
#'
#' @description The function to create the ui module for viewing LSAE and peak height parameters
view_par_ui <- function(id){
  ns <- NS(id)

  tabPanel(
    "View LSAE & peak height parameters",
    titlePanel("View LSAE & peak height parameters"),
    br(),
    splitLayout(
      cellWidths = c("500px"),
      wellPanel(
        fluidRow(
          column(3, br(), shinyFilesButton(ns("load_param_data"), label = "Load", title = "Select a file", multiple = FALSE)),
          column(9, uiOutput(ns("fn_param_data")))
        )
      )
    ),
    br(),
    fluidRow(
      column(
        2,
        wellPanel(
          h4("Expected LSAE"),
          dataTableOutput(ns("dt_expect_ae"))
        )
      ),
      column(
        4,
        wellPanel(
          h4("Experimental LSAE data"),
          uiOutput(ns("ui_plot_lsae")),
          br(),
          dataTableOutput(ns("dt_lsae_data"))
        )
      ),
      column(
        2,
        wellPanel(
          h4("Peak height parameters"),
          dataTableOutput(ns("dt_est_gamma"))
        )
      ),
      column(
        4,
        wellPanel(
          h4("P-P plot for peak height parameters")
        )
      )
    )
  )
}

#' view_par_server
#'
#' @description The function to create the server module for viewing LSAE and peak height parameters
view_par_server <- function(id, path_pack){
  moduleServer(
    id,
    function(input, output, session){
      rv_view_par <- reactiveValues()
      rv_view_par$fn_param_data <- NULL
      rv_view_par$dt_prop_ae_data_use <- NULL
      rv_view_par$dt_expect_ae <- NULL
      rv_view_par$dt_est_gamma <- NULL
      rv_view_par$df_select_lsae <- NULL

      shinyFileChoose(input, "load_param_data", roots = c(documents = paste0(path_pack, "/extdata/param")), session = session, filetypes = c("RData"))

      observeEvent(input$load_param_data, {
        req(input$load_param_data)
        path_param_data <- parseFilePaths(roots = c(documents = paste0(path_pack, "/extdata/param")), input$load_param_data)$datapath

        if(length(path_param_data) > 0){
          rv_view_par$fn_param_data <- basename(path_param_data)

          load(path_param_data)

          rv_view_par$dt_prop_ae_data_use <- fread("C:/Users/manab/Desktop/dt_prop_ae_data_use.csv")

          rv_view_par$dt_expect_ae <- data.table(Marker = names(expect_ae),
                                                 Alpha = alpha,
                                                 LSAE = expect_ae)

          rv_view_par$dt_est_gamma <- data.table(Parameter = c("CV for allele", "CV for stutter"),
                                                 Value = c(cv_allele, cv_stutter))
        }
      })

      output$fn_param_data <- renderUI({
        disabled(textInput(session$ns("fn_param_data"), label = "", value = rv_view_par$fn_param_data))
      })

      output$dt_expect_ae <- renderDataTable({
        dt_expect_ae <- rv_view_par$dt_expect_ae

        if(!is.null(dt_expect_ae)){
          datatable(
            dt_expect_ae,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_expect_ae), ordering = FALSE, autoWidth = TRUE, scrollY = "600px", fixedHeader = TRUE),
            rownames = FALSE
          ) %>%
            formatSignif(columns = c("Alpha", "LSAE"), digits = 4) %>%
            formatStyle(names(dt_expect_ae), lineHeight = "80%")
        }
      })

      output$dt_est_gamma <- renderDataTable({
        dt_est_gamma <- rv_view_par$dt_est_gamma

        if(!is.null(dt_est_gamma)){
          datatable(
            dt_est_gamma,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_est_gamma), ordering = FALSE, autoWidth = TRUE, fixedHeader = TRUE),
            rownames = FALSE
          ) %>%
            formatSignif(columns = c("Value"), digits = 4) %>%
            formatStyle(names(dt_est_gamma), lineHeight = "80%")
        }
      })

      output$ui_plot_lsae <- renderUI({
        dt_prop_ae_data_use <- rv_view_par$dt_prop_ae_data_use

        if(!is.null(dt_prop_ae_data_use)){
          plotOutput(session$ns("plot_lsae"), width = "100%", height = "300px")
        }
      })

      output$plot_lsae <- renderPlot({
        dt_prop_ae_data_use <- rv_view_par$dt_prop_ae_data_use

        if(!is.null(dt_prop_ae_data_use)){
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

          df_select_lsae <- rv_view_par$df_select_lsae
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
              plot.title = element_text(size = 16),
              legend.text = element_text(size = 14),
              axis.title.x = element_text(size = 14, margin = margin(t = 20)),
              axis.title.y = element_text(size = 14, margin = margin(r = 20)),
              axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5),
              axis.text.y = element_text(size = 14)
            )
        }
      })

      output$dt_lsae_data <- renderDataTable({
        dt_prop_ae_data_use <- rv_view_par$dt_prop_ae_data_use

        if(!is.null(dt_prop_ae_data_use)){
          markers <- colnames(dt_prop_ae_data_use)[-1]

          datatable(
            dt_prop_ae_data_use,
            selection = list(mode = "single", target = "row"),
            extensions = "FixedColumns",
            options = list(dom = "t", iDisplayLength = nrow(dt_prop_ae_data_use), autoWidth = TRUE, scrollX = TRUE, scrollY = "300px",
                           fixedHeader = TRUE, fixedColumns = list(leftColumns = 1)
            ),
            rownames = FALSE
          ) %>%
            formatSignif(columns = markers, digits = 4) %>%
            formatStyle(names(dt_prop_ae_data_use), lineHeight = "80%")
        }
      })

      observeEvent(input$dt_lsae_data_rows_selected, {
        index_select <- input$dt_lsae_data_rows_selected
        if(length(index_select) == 1){
          dt_prop_ae_data_use <- rv_view_par$dt_prop_ae_data_use
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

          rv_view_par$df_select_lsae <- df
        }else{
          rv_view_par$df_select_lsae <- NULL
        }
      })
    }
  )
}
