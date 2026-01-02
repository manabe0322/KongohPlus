#' view_lsae_ui
#'
#' @description The function to create the ui module for viewing the locus-specific amplification efficiencies
view_lsae_ui <- function(id){
  ns <- NS(id)

  tabPanel(
    "View locus-specific amplification efficiencies",
    titlePanel("View locus-specific amplification efficiencies"),
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

#' view_lsae_server
#'
#' @description The function to create the server module for viewing the locus-specific amplification efficiencies
view_lsae_server <- function(id, path_pack){
  moduleServer(
    id,
    function(input, output, session){
      rv_view_lsae <- reactiveValues()
      rv_view_lsae$fn_lsae_data <- NULL
      rv_view_lsae$dt_prop_ae_data_use <- NULL
      rv_view_lsae$dt_expect_ae <- NULL
      rv_view_lsae$df_select_lsae <- NULL

      shinyFileChoose(input, "load_lsae_data", roots = c(documents = paste0(path_pack, "/extdata/lsae")), session = session, filetypes = c("RData"))

      observeEvent(input$load_lsae_data, {
        req(input$load_lsae_data)
        path_lsae_data <- parseFilePaths(roots = c(documents = paste0(path_pack, "/extdata/lsae")), input$load_lsae_data)$datapath

        if(length(path_lsae_data) > 0){
          rv_view_lsae$fn_lsae_data <- basename(path_lsae_data)

          load(path_lsae_data)

          rv_view_lsae$dt_prop_ae_data_use <- data_est_lsae$dt_prop_ae_data_use
          rv_view_lsae$dt_expect_ae <- data_est_lsae$dt_expect_ae
        }
      })

      output$fn_lsae_data <- renderUI({
        disabled(textInput(session$ns("fn_lsae_data"), label = "", value = rv_view_lsae$fn_lsae_data))
      })

      output$dt_expect_ae <- renderDataTable({
        dt_expect_ae <- rv_view_lsae$dt_expect_ae

        if(!is.null(dt_expect_ae)){
          datatable(
            dt_expect_ae,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_expect_ae), ordering = FALSE, autoWidth = FALSE),
            rownames = FALSE
          ) %>%
            formatSignif(columns = c("Alpha", "LSAE"), digits = 4) %>%
            formatStyle(names(dt_expect_ae), lineHeight = "80%")
        }
      })

      output$ui_plot_lsae <- renderUI({
        dt_prop_ae_data_use <- rv_view_lsae$dt_prop_ae_data_use

        if(!is.null(dt_prop_ae_data_use)){
          plotOutput(session$ns("plot_lsae"), width = "100%", height = "300px")
        }
      })

      output$plot_lsae <- renderPlot({
        dt_prop_ae_data_use <- rv_view_lsae$dt_prop_ae_data_use

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
        }
      })

      output$dt_lsae_data <- renderDataTable({
        dt_prop_ae_data_use <- rv_view_lsae$dt_prop_ae_data_use

        if(!is.null(dt_prop_ae_data_use)){
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
        }
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
