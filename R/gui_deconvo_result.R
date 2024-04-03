#' deconvo_result_ui
#'
#' @description The function to create the ui module for results of deconvolution
deconvo_result_ui <- function(id){
  ns <- NS(id)

  tabPanel("Result",
           br(),
           fluidRow(
             column(3,
                    wellPanel(
                      h4("Estimated parameter"),
                      br(),
                      dataTableOutput(ns("dt_mpdeg")),
                      br(),
                      textOutput(ns("tph")),
                      br(),
                      textOutput(ns("lambda")),
                      br(),
                      actionButton(ns("act_ae"), label = "AE")
                    ),
                    wellPanel(
                      h4("Diagnostics"),
                      br(),
                      actionButton(ns("act_rhat"), label = "Rhat"),
                      br(),
                      actionButton(ns("act_traceplot_lp"), label = "Trace plot")
                    ),
                    br(),
                    actionButton(ns("act_proceed_lr"), label = "Likelihood ratio", class = "btn btn-primary btn-lg"),
                    br(),
                    downloadButton(ns("download_report"), "Report")
             ),
             column(6,
                    wellPanel(
                      h4("Probabilistic genotyping"),
                      br(),
                      uiOutput(ns("mk_pg")),
                      br(),
                      dataTableOutput(ns("dt_pg")),
                      br()
                    )
             ),
             column(3,
                    wellPanel(
                      h4("Likelihood"),
                      br(),
                      dataTableOutput(ns("dt_log10like")),
                      br()
                    )
             )
           )
  )
}

#' deconvo_result_server
#'
#' @description The function to create the server module for results of deconvolution
deconvo_result_server <- function(id, path_pack, rv_deconvo){
  moduleServer(
    id,
    function(input, output, session){
      data_deconvo <- rv_deconvo$data_deconvo
      if(!is.null(data_deconvo)){
        dt_csp <- rv_deconvo$dt_csp
        dt_ref <- rv_deconvo$dt_ref
        dt_af <- rv_deconvo$dt_af

        hnc <- data_deconvo$hnc
        markers_auto <- data_deconvo$markers_auto
        dt_pg_list <- data_deconvo$dt_pg_list
        dt_log10like <- data_deconvo$dt_log10like
        dt_mpdeg <- data_deconvo$dt_mpdeg
        tph <- data_deconvo$tph
        lambda <- data_deconvo$lambda
        dt_ae <- data_deconvo$dt_ae
        dt_rhat <- data_deconvo$dt_rhat
        traceplot <- data_deconvo$traceplot

        output$dt_mpdeg <- renderDataTable(server = FALSE, {
          datatable(
            dt_mpdeg,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_mpdeg), ordering = FALSE, autowidth = TRUE),
            rownames = FALSE
          ) %>%
            formatSignif(columns = c("MP", "Deg"), digits = 3)
        })

        output$tph <- renderText({paste0("tph: ", round(tph, 0))})

        output$lambda <- renderText({paste0("lambda: ", paste(round(lambda, 3), collapse = ", "))})

        observeEvent(input$act_ae, {
          showModal(modalDialog(
            title = "AE",
            dataTableOutput(session$ns("dt_ae")),
            easyClose = TRUE,
            footer = NULL
          ))
        })

        output$dt_ae <- renderDataTable(server = FALSE, {
          datatable(
            dt_ae,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_ae), ordering = FALSE, autoWidth = TRUE),
            rownames = FALSE
          ) %>%
            formatRound(columns = c("AE"), digits = 3)
        })

        observeEvent(input$act_rhat, {
          showModal(modalDialog(
            title = "Rhat",
            dataTableOutput(session$ns("dt_rhat")),
            easyClose = TRUE,
            footer = NULL
          ))
        })

        output$dt_rhat <- renderDataTable(server = FALSE, {
          datatable(
            dt_rhat,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_rhat), ordering = FALSE, autoWidth = TRUE),
            rownames = FALSE
          ) %>%
            formatRound(columns = c("Rhat"), digits = 3)
        })

        observeEvent(input$act_traceplot_lp, {
          showModal(modalDialog(
            title = "Trace plot",
            plotOutput(session$ns("traceplot_lp")),
            footer = tagList(
              modalButton("Close")
            ),
            size = "l"
          ))
        })

        output$download_report <- downloadHandler(
          filename = "report.html",
          content = function(file){
            report_tmp <- file.path(tempdir(), "report.Rmd")
            file.copy(paste0(path_pack, "/rmd/report.Rmd"), report_tmp, overwrite = TRUE)
            params <- list(dt_csp = dt_csp, dt_ref = dt_ref, hnc = hnc)
            render(report_tmp, output_file = file, params = params, envir = new.env(parent = globalenv()))
          }
        )

        output$traceplot_lp <- renderPlot({
          traceplot
        })

        output$mk_pg <- renderUI({
          selectInput(session$ns("mk_pg"), label = "Marker", choices = markers_auto, selected = markers_auto[1])
        })

        observeEvent(input$mk_pg, {
          mk_select <- input$mk_pg

          output$dt_pg <- renderDataTable(server = FALSE, {
            datatable(
              dt_pg_list[[mk_select]],
              selection = list(mode = "single", target = "row"),
              options = list(iDisplayLength = 25, ordering = FALSE, autowidth = TRUE,
                             columnDefs = list(list(targets = c(hnc + 1, hnc + 2), visible = FALSE))),
              rownames = FALSE
            ) %>%
              formatSignif(columns = c("Weight"), digits = 3)
          })
        })

        output$dt_log10like <- renderDataTable(server = FALSE, {
          datatable(
            dt_log10like,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_log10like), ordering = FALSE, autowidth = TRUE),
            rownames = FALSE
          ) %>%
            formatSignif(columns = c("Log10_Likelihood"), digits = 3)
        })
      }
    }
  )
}
