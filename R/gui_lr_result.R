#' lr_result_ui
#'
#' @description The function to create the ui module for results of the likelihood ratio
lr_result_ui <- function(id){
  ns <- NS(id)

  tabPanel("Result",
           br(),
           fluidRow(
             column(2,
                    wellPanel(
                      uiOutput(ns("hp")),
                      uiOutput(ns("noc_lr_hp")),
                      uiOutput(ns("hd")),
                      uiOutput(ns("noc_lr_hd"))
                    ),
                    wellPanel(
                      h4(textOutput(ns("title_warnings"))),
                      htmlOutput(ns("warning_message"))
                    )
             ),
             column(2,
                    wellPanel(
                      h4(textOutput(ns("title_lr"))),
                      dataTableOutput(ns("dt_lr")),
                      br()
                    )
             ),
             column(4,
                    wellPanel(
                      h4(textOutput(ns("title_pg"))),
                      fluidRow(
                        column(6,
                               uiOutput(ns("select_meth_pg")),
                               uiOutput(ns("select_mk_pg"))
                        ),
                        column(6,
                               uiOutput(ns("select_noc_pg")),
                               uiOutput(ns("select_contrib_pg"))
                        )
                      ),
                      dataTableOutput(ns("dt_pg"))
                    )
             ),
             column(4,
                    wellPanel(
                      h4(textOutput(ns("title_par_diag"))),
                      fluidRow(
                        column(6,
                               uiOutput(ns("select_noc_diag"))
                        ),
                        column(6,
                               uiOutput(ns("button_traceplot"))
                        )
                      ),
                      dataTableOutput(ns("dt_diag")),
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
        dt_lr <- data_lr$dt_lr
        refs <- data_lr$refs
        hp <- data_lr$hp
        hd <- data_lr$hd

        data_deconvo_hp <- rv_lr$data_deconvo_hp
        markers_use <- data_deconvo_hp$markers_use
        message_sr <- data_deconvo_hp$message_sr
        hnc_hp <- data_deconvo_hp$hnc
        dt_pg_list_hp <- data_deconvo_hp$dt_pg_list
        traceplot_hp <- data_deconvo_hp$traceplot
        dt_diag_hp <- data_deconvo_hp$dt_diag

        data_deconvo_hd <- rv_lr$data_deconvo_hd
        hnc_hd <- data_deconvo_hd$hnc
        dt_pg_list_hd <- data_deconvo_hd$dt_pg_list
        traceplot_hd <- data_deconvo_hd$traceplot
        dt_diag_hd <- data_deconvo_hd$dt_diag

        rv_tmp <- reactiveValues()
        rv_tmp$update <- FALSE
        rv_tmp$cand_contrib <- NULL

        output$title_warnings <- renderText({"Warnings"})
        output$title_lr <- renderText({"Likelihood ratio"})
        output$title_pg <- renderText({"Probabilistic genotyping"})
        output$title_par_diag <- renderText({"Parameters and diagnostics"})

        output$hp <- renderUI({
          disabled(checkboxGroupInput(session$ns("hp"), label = h4("Prosecutor hypothesis"), choices = refs, selected = hp))
        })

        output$noc_lr_hp <- renderUI({
          disabled(numericInput(session$ns("noc_lr_hp"), label = "Number of contributors", value = hnc_hp))
        })

        output$hd <- renderUI({
          disabled(checkboxGroupInput(session$ns("hd"), label = h4("Defense hypothesis"), choices = refs, selected = hd))
        })

        output$noc_lr_hd <- renderUI({
          disabled(numericInput(session$ns("noc_lr_hd"), label = "Number of contributors", value = hnc_hd))
        })

        output$warning_message <- renderUI({HTML(message_sr)})

        ## Probabilistic genotyping ##

        output$select_meth_pg <- renderUI({
          selectInput(session$ns("select_meth_pg"), label = "Method", choices = c("Per combination", "Per contributor"), selected = "Per combination", selectize = FALSE)
        })

        output$select_noc_pg <- renderUI({
          cand_noc <- sort(unique(c(hnc_hp, hnc_hd)))
          selectInput(session$ns("select_noc_pg"), label = "Number of contributors", choices = cand_noc, selected = cand_noc[1], selectize = FALSE)
        })

        output$select_mk_pg <- renderUI({
          selectInput(session$ns("select_mk_pg"), label = "Marker", choices = markers_use, selected = markers_use[1], selectize = FALSE)
        })

        output$select_contrib_pg <- renderUI({
          selectInput(session$ns("select_contrib_pg"), label = "Contributor", choices = rv_tmp$cand_contrib, selected = rv_tmp$cand_contrib[1], selectize = FALSE)
        })

        observeEvent(input$select_meth_pg, {
          rv_tmp$update <- TRUE
        })

        observeEvent(input$select_noc_pg, {
          rv_tmp$update <- TRUE
        })

        observeEvent(input$select_mk_pg, {
          rv_tmp$update <- TRUE
        })

        observeEvent(input$select_contrib_pg, {
          rv_tmp$update <- TRUE
        })

        observe({
          if(rv_tmp$update){
            hnc <- as.numeric(input$select_noc_pg)
            rv_tmp$cand_contrib <- paste0("U", 1:hnc)
            if(hnc == hnc_hp){
              dt_pg_display <- dt_pg_list_hp
            }else{
              dt_pg_display <- dt_pg_list_hd
            }

            meth_pg <- input$select_meth_pg
            if(meth_pg == "Per combination"){
              disable("select_contrib_pg")

              output$dt_pg <- renderDataTable(server = FALSE, {
                datatable(
                  dt_pg_display[[input$select_mk_pg]],
                  selection = list(mode = "single", target = "row"),
                  options = list(iDisplayLength = 10, ordering = FALSE, autowidth = TRUE, searching = FALSE,
                                 columnDefs = list(list(targets = c(hnc, hnc + 1, hnc + 2), visible = FALSE))),
                  rownames = FALSE
                ) %>%
                  formatSignif(columns = c("Probability"), digits = 3)
              })
            }else{
              enable("select_contrib_pg")

              dt_pg_select <- dt_pg_display[[input$select_mk_pg]]
              aggregated_data <- aggregate_pg_contrib(dt_pg_select, input$select_contrib_pg)

              output$dt_pg <- renderDataTable(server = FALSE, {
                datatable(
                  aggregated_data,
                  selection = list(mode = "single", target = "row"),
                  options = list(iDisplayLength = 10, ordering = FALSE, autowidth = TRUE, searching = FALSE),
                  rownames = FALSE
                ) %>%
                  formatSignif(columns = c("Probability"), digits = 3)
              })
            }

            rv_tmp$update <- FALSE
          }
        })

        ## Parameters and diagnostics

        output$select_noc_diag <- renderUI({
          cand_noc <- sort(unique(c(hnc_hp, hnc_hd)))
          selectInput(session$ns("select_noc_diag"), label = "Number of contributors", choices = cand_noc, selected = cand_noc[1], multiple = FALSE, selectize = FALSE)
        })

        output$button_traceplot <- renderUI({
          actionButton(session$ns("act_traceplot_lp"), label = "Trace plot")
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

        output$traceplot_lp <- renderPlot({
          noc_diag <- input$select_noc_diag
          if(noc_diag == hnc_hp){
            traceplot <- traceplot_hp
          }else{
            traceplot <- traceplot_hd
          }
          traceplot
        })

        observeEvent(input$select_noc_diag, {
          noc_diag <- input$select_noc_diag
          if(noc_diag == hnc_hp){
            dt_diag <- dt_diag_hp
          }else{
            dt_diag <- dt_diag_hd
          }

          output$dt_diag <- renderDataTable({
            datatable(
              dt_diag,
              selection = "none",
              options = list(dom = "t", iDisplayLength = nrow(dt_diag), ordering = FALSE, autowidth = TRUE,
                             columnDefs = list(list(width = "175px", targets = 0),
                                               list(width = "90px", targets = 1),
                                               list(width = "100px", targets = 2),
                                               list(width = "125px", targets = 3),
                                               list(width = "50px", targets = 4))),
              rownames = FALSE
            ) %>%
              formatSignif(columns = c("Estimated value", "Rhat"), digits = 4) %>%
              formatStyle(names(dt_diag), lineHeight = "80%")
          })
        })

        output$dt_lr <- renderDataTable({
          datatable(
            dt_lr,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_lr), ordering = FALSE, autowidth = TRUE,
                           columnDefs = list(list(targets = c(1, 2), visible = FALSE))),
            rownames = FALSE
          ) %>%
            formatSignif(columns = c("LR"), digits = 3) %>%
            formatStyle(names(dt_lr), lineHeight = "80%")
        })
      }
    }
  )
}
