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
                      h4("Parameters for contributors"),
                      br(),
                      dataTableOutput(ns("dt_mpdeg")),
                    ),
                    br(),
                    actionButton(ns("act_pardiag"), label = "Parameters and diagnostics"),
                    br(),
                    br(),
                    actionButton(ns("act_traceplot_lp"), label = "Trace plot"),
                    br(),
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
        markers_use <- data_deconvo$markers_use
        dt_pg_list <- data_deconvo$dt_pg_list
        dt_log10like <- data_deconvo$dt_log10like
        dt_mpdeg <- data_deconvo$dt_mpdeg
        dt_mcmc <- data_deconvo$dt_mcmc
        traceplot <- data_deconvo$traceplot

        output$dt_mpdeg <- renderDataTable(server = FALSE, {
          datatable(
            dt_mpdeg,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_mpdeg), ordering = FALSE, autowidth = TRUE),
            rownames = FALSE,
            colnames = c("Contributor", "Mixture proportion", "Degradation")
          ) %>%
            formatSignif(columns = c("MP", "Deg"), digits = 2)
        })

        observeEvent(input$act_pardiag, {
          showModal(modalDialog(
            title = "Parameters and diagnostics",
            dataTableOutput(session$ns("dt_mcmc")),
            easyClose = TRUE,
            footer = NULL,
            size = "l"
          ))
        })

        output$dt_mcmc <- renderDataTable(server = FALSE, {
          datatable(
            dt_mcmc,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_mcmc), ordering = FALSE, autowidth = TRUE,
                           columnDefs = list(list(width = "200px", targets = c(0)),
                                             list(width = "150px", targets = c(1, 2, 3, 4)))),
            rownames = FALSE,
            colnames = c("Parameter", "Contributor", "Marker", "Estimated value", "Rhat")
          ) %>%
            formatSignif(columns = c("Estimated_value", "Rhat"), digits = 3)
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
            params <- list(dt_csp = dt_csp, dt_ref = dt_ref, hnc = hnc, markers_use = markers_use, dt_mpdeg = dt_mpdeg, dt_mcmc = dt_mcmc, dt_pg_list = dt_pg_list, dt_log10like = dt_log10like)
            render(report_tmp, output_file = file, params = params, envir = new.env(parent = globalenv()))
          }
        )

        output$traceplot_lp <- renderPlot({
          traceplot
        })

        output$mk_pg <- renderUI({
          selectInput(session$ns("mk_pg"), label = "Marker", choices = markers_use, selected = markers_use[1], selectize = FALSE)
        })

        observeEvent(input$mk_pg, {
          mk_select <- input$mk_pg

          output$dt_pg <- renderDataTable(server = FALSE, {
            datatable(
              dt_pg_list[[mk_select]],
              selection = list(mode = "single", target = "row"),
              options = list(iDisplayLength = 10, ordering = FALSE, autowidth = TRUE,
                             columnDefs = list(list(targets = c(hnc, hnc + 1, hnc + 2), visible = FALSE))),
              rownames = FALSE
            ) %>%
              formatSignif(columns = c("post_prob"), digits = 3)
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
