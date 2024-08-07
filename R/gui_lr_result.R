#' lr_result_ui
#'
#' @description The function to create the ui module for results of the likelihood ratio
lr_result_ui <- function(id){
  ns <- NS(id)

  tabPanel("Result",
           br(),
           fluidRow(
             column(3,
                    wellPanel(
                      h4("Crime stain profile"),
                      textOutput(ns("fn_csp")),
                      br(),
                      h4("Reference profile"),
                      textOutput(ns("fn_ref")),
                      br(),
                      h4("Allele frequency"),
                      textOutput(ns("fn_af"))
                    ),
                    wellPanel(
                      h4("Parameters for contributors"),
                      br(),
                      p("Hp"),
                      dataTableOutput(ns("dt_mpdeg_hp")),
                      br(),
                      p("Hd"),
                      dataTableOutput(ns("dt_mpdeg_hd"))
                    ),
                    br(),
                    actionButton(ns("act_pardiag_hp"), label = "Parameters and diagnostics (Hp)"),
                    br(),
                    br(),
                    actionButton(ns("act_pardiag_hd"), label = "Parameters and diagnostics (Hd)"),
                    br(),
                    br(),
                    actionButton(ns("act_traceplot_lp_hp"), label = "Trace plot (Hp)"),
                    br(),
                    br(),
                    actionButton(ns("act_traceplot_lp_hd"), label = "Trace plot (Hd)")
             ),
             column(6,
                    wellPanel(
                      h4("Probabilistic genotyping"),
                      br(),
                      uiOutput(ns("hyp_pg")),
                      br(),
                      uiOutput(ns("mk_pg")),
                      br(),
                      dataTableOutput(ns("dt_pg")),
                      br()
                    )
             ),
             column(3,
                    wellPanel(
                      h4("Likelihood ratio"),
                      br(),
                      dataTableOutput(ns("dt_lr")),
                      br()
                    )
             )
           )
  )
}

#' lr_result_server
#'
#' @description The function to create the server module for results of the likelihood ratio
lr_result_server <- function(id, rv_lr){
  moduleServer(
    id,
    function(input, output, session){
      data_lr <- rv_lr$data_lr
      if(!is.null(data_lr)){
        fn_ref <- data_lr$fn_ref
        dt_mpdeg_hp <- data_lr$dt_mpdeg_hp
        dt_mpdeg_hd <- data_lr$dt_mpdeg_hd
        dt_pg_list_hp <- data_lr$dt_pg_list_hp
        dt_pg_list_hd <- data_lr$dt_pg_list_hd
        dt_lr <- data_lr$dt_lr

        data_deconvo_hp <- rv_lr$data_deconvo_hp
        hnc_hp <- data_deconvo_hp$hnc
        dt_mcmc_hp <- data_deconvo_hp$dt_mcmc
        traceplot_hp <- data_deconvo_hp$traceplot
        markers_use <- data_deconvo_hp$markers_use
        fn_csp <- data_deconvo_hp$fn_csp
        fn_af <- data_deconvo_hp$fn_af

        data_deconvo_hd <- rv_lr$data_deconvo_hd
        hnc_hd <- data_deconvo_hd$hnc
        dt_mcmc_hd <- data_deconvo_hd$dt_mcmc
        traceplot_hd <- data_deconvo_hd$traceplot

        rv_lr_result <- reactiveValues()
        rv_lr_result$dt_pg_display <- dt_pg_list_hp

        output$fn_csp <- renderText({paste0(fn_csp)})
        output$fn_ref <- renderText({paste0(fn_ref)})
        output$fn_af <- renderText({paste0(fn_af)})

        output$dt_mpdeg_hp <- renderDataTable(server = FALSE, {
          datatable(
            dt_mpdeg_hp,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_mpdeg_hp), ordering = FALSE, autowidth = TRUE),
            rownames = FALSE,
            colnames = c("Contributor", "Mixture proportion", "Degradation")
          ) %>%
            formatSignif(columns = c("MP", "Deg"), digits = 3)
        })

        output$dt_mpdeg_hd <- renderDataTable(server = FALSE, {
          datatable(
            dt_mpdeg_hd,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_mpdeg_hd), ordering = FALSE, autowidth = TRUE),
            rownames = FALSE,
            colnames = c("Contributor", "Mixture proportion", "Degradation")
          ) %>%
            formatSignif(columns = c("MP", "Deg"), digits = 3)
        })

        observeEvent(input$act_pardiag_hp, {
          showModal(modalDialog(
            title = "Parameters and diagnostics in Hp",
            dataTableOutput(session$ns("dt_mcmc_hp")),
            easyClose = TRUE,
            footer = NULL,
            size = "l"
          ))
        })

        output$dt_mcmc_hp <- renderDataTable(server = FALSE, {
          datatable(
            dt_mcmc_hp,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_mcmc_hp), ordering = FALSE, autowidth = TRUE,
                           columnDefs = list(list(width = "200px", targets = c(0)),
                                             list(width = "100px", targets = c(1, 2, 3, 4)))),
            rownames = FALSE,
            colnames = c("Parameter", "Contributor", "Marker", "Estimated value", "Rhat")
          ) %>%
            formatSignif(columns = c("Estimated_value", "Rhat"), digits = 3)
        })

        observeEvent(input$act_pardiag_hd, {
          showModal(modalDialog(
            title = "Parameters and diagnostics in Hd",
            dataTableOutput(session$ns("dt_mcmc_hd")),
            easyClose = TRUE,
            footer = NULL,
            size = "l"
          ))
        })

        output$dt_mcmc_hd <- renderDataTable(server = FALSE, {
          datatable(
            dt_mcmc_hd,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_mcmc_hd), ordering = FALSE, autowidth = TRUE,
                           columnDefs = list(list(width = "200px", targets = c(0)),
                                             list(width = "100px", targets = c(1, 2, 3, 4)))),
            rownames = FALSE,
            colnames = c("Parameter", "Contributor", "Marker", "Estimated value", "Rhat")
          ) %>%
            formatSignif(columns = c("Estimated_value", "Rhat"), digits = 3)
        })

        observeEvent(input$act_traceplot_lp_hp, {
          showModal(modalDialog(
            title = "Trace plot (Hp)",
            plotOutput(session$ns("traceplot_lp_hp")),
            footer = tagList(
              modalButton("Close")
            ),
            size = "l"
          ))
        })

        output$traceplot_lp_hp <- renderPlot({
          traceplot_hp
        })

        observeEvent(input$act_traceplot_lp_hd, {
          showModal(modalDialog(
            title = "Trace plot (Hd)",
            plotOutput(session$ns("traceplot_lp_hd")),
            footer = tagList(
              modalButton("Close")
            ),
            size = "l"
          ))
        })

        output$traceplot_lp_hd <- renderPlot({
          traceplot_hd
        })

        output$hyp_pg <- renderUI({
          selectInput(session$ns("hyp_pg"), label = "Hypothesis", choices = c("Hp", "Hd"), selected = "Hp", selectize = FALSE)
        })

        output$mk_pg <- renderUI({
          selectInput(session$ns("mk_pg"), label = "Marker", choices = markers_use, selected = markers_use[1], selectize = FALSE)
        })

        observeEvent(input$hyp_pg, {
          if(input$hyp_pg == "Hp"){
            rv_lr_result$dt_pg_display <- dt_pg_list_hp
            hnc <- hnc_hp
          }else{
            rv_lr_result$dt_pg_display <- dt_pg_list_hd
            hnc <- hnc_hd
          }

          mk_select <- input$mk_pg
          dt_pg_display <- rv_lr_result$dt_pg_display

          output$dt_pg <- renderDataTable(server = FALSE, {
            datatable(
              dt_pg_display[[mk_select]],
              selection = list(mode = "single", target = "row"),
              options = list(iDisplayLength = 25, ordering = FALSE, autowidth = TRUE,
                             columnDefs = list(list(targets = c(hnc, hnc + 1, hnc + 2), visible = FALSE))),
              rownames = FALSE
            ) %>%
              formatSignif(columns = c("Probability"), digits = 3)
          })
        })

        observeEvent(input$mk_pg, {
          if(input$hyp_pg == "Hp"){
            hnc <- hnc_hp
          }else{
            hnc <- hnc_hd
          }

          mk_select <- input$mk_pg
          dt_pg_display <- rv_lr_result$dt_pg_display

          output$dt_pg <- renderDataTable(server = FALSE, {
            datatable(
              dt_pg_display[[mk_select]],
              selection = list(mode = "single", target = "row"),
              options = list(iDisplayLength = 25, ordering = FALSE, autowidth = TRUE,
                             columnDefs = list(list(targets = c(hnc, hnc + 1, hnc + 2), visible = FALSE))),
              rownames = FALSE
            ) %>%
              formatSignif(columns = c("Probability"), digits = 3)
          })
        })

        output$dt_lr <- renderDataTable(server = FALSE, {
          datatable(
            dt_lr,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_lr), ordering = FALSE, autowidth = TRUE,
                           columnDefs = list(list(targets = c(1, 2), visible = FALSE))),
            rownames = FALSE
          ) %>%
            formatSignif(columns = c("LR"), digits = 3)
        })
      }
    }
  )
}
