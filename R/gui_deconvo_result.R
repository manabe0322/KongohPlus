#' deconvo_result_ui
#'
#' @description The function to create the ui module for results of deconvolution
deconvo_result_ui <- function(id){
  ns <- NS(id)

  tabPanel("Result",
           br(),
           fluidRow(
             column(2,
                    wellPanel(
                      h4(textOutput(ns("title_warnings"))),
                      htmlOutput(ns("warning_message"))
                    ),
                    br(),
                    uiOutput(ns("act_lr"))
             ),
             column(4,
                    wellPanel(
                      h4(textOutput(ns("title_like"))),
                      dataTableOutput(ns("dt_log10like")),
                      br()
                    )
             ),
             column(6,
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
                      uiOutput(ns("diag")),
                      br()
                    )
             )
           )
  )
}

#' deconvo_result_server
#'
#' @description The function to create the server module for results of deconvolution
deconvo_result_server <- function(id, session_top, rv_deconvo){
  moduleServer(
    id,
    function(input, output, session){
      data_deconvo <- rv_deconvo$data_deconvo
      if(!is.null(data_deconvo)){
        markers_use <- data_deconvo$markers_use
        message_sr <- data_deconvo$message_sr
        hnc_all <- data_deconvo$hnc_all
        data_hnc_list <- data_deconvo$data_hnc_list

        rv_deconvo_result <- reactiveValues()
        rv_deconvo_result$dt_diag <- NULL
        rv_deconvo_result$error_message <- NULL
        rv_deconvo_result$traceplot <- NULL

        output$title_warnings <- renderText({"Warnings"})
        output$title_like <- renderText({"Likelihood"})
        output$title_pg <- renderText({"Probabilistic genotyping"})
        output$title_par_diag <- renderText({"Parameters and diagnostics"})

        output$warning_message <- renderUI({HTML(message_sr)})

        output$act_lr <- renderUI({actionButton(session$ns("act_lr"), label = "Likelihood ratio", class = "btn btn-primary btn-lg")})

        observeEvent(input$act_lr, {
          updateNavbarPage(session = session_top, inputId = "navbar", selected = "Likelihood ratio")
        })

        rowCallback <- c(
          "function(row, data){",
          "  for(var i=0; i<data.length; i++){",
          "    if(data[i] === null){",
          "      $('td:eq('+i+')', row).html('NA')",
          "        .css({'color': 'rgb(151,151,151)', 'font-style': 'italic'});",
          "    }",
          "  }",
          "}"
        )

        output$dt_log10like <- renderDataTable(server = FALSE, {
          n_hnc <- length(hnc_all)
          mat_tmp <- matrix(NA, length(markers_use) + 1, n_hnc)
          for(i in 1:n_hnc){
            data_hnc <- data_hnc_list[[i]]
            if(is.null(data_hnc$error_message)){
              mat_tmp[, i] <- round(data_hnc$dt_log10like$Log10_Likelihood, 2)
            }
          }
          dt_log10like_display <- data.table(Marker = c(markers_use, "Total"))
          dt_log10like_display <- cbind(dt_log10like_display, as.data.table(mat_tmp))
          names(dt_log10like_display) <- c("Marker", paste0("NoC = ", hnc_all))

          datatable(
            dt_log10like_display,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_log10like_display), ordering = FALSE, autowidth = TRUE, rowCallback = JS(rowCallback)),
            rownames = FALSE,
            escape = FALSE
          ) %>%
            formatStyle(names(dt_log10like_display), lineHeight = "80%")
        })

        ## Parameters and diagnostics

        output$select_noc_diag <- renderUI({
          selectInput(session$ns("select_noc_diag"), label = "Number of contributors", choices = hnc_all, selected = hnc_all[1], multiple = FALSE, selectize = FALSE)
        })

        observeEvent(input$select_noc_diag, {
          index_hnc <- which(hnc_all == input$select_noc_diag)
          data_hnc <- data_hnc_list[[index_hnc]]
          error_message <- data_hnc$error_message
          if(is.null(data_hnc$error_message)){
            rv_deconvo_result$dt_diag <- data_hnc$dt_diag
            rv_deconvo_result$error_message <- NULL
            rv_deconvo_result$traceplot <- data_hnc$traceplot
            enable("act_traceplot_lp")
          }else{
            rv_deconvo_result$dt_diag <- NULL
            rv_deconvo_result$error_message <- error_message
            rv_deconvo_result$traceplot <- NULL
            disable("act_traceplot_lp")
          }
        })

        output$diag <- renderUI({
          dt_diag <- rv_deconvo_result$dt_diag
          if(is.null(dt_diag)){
            HTML(rv_deconvo_result$error_message)
          }else{
            dataTableOutput(session$ns("dt_diag"))
          }
        })

        output$dt_diag <- renderDataTable({
          dt_diag <- rv_deconvo_result$dt_diag

          if(!is.null(dt_diag)){
            datatable(
              dt_diag,
              selection = "none",
              options = list(dom = "t", iDisplayLength = nrow(dt_diag), ordering = FALSE, autowidth = TRUE,
                             columnDefs = list(list(width = "200px", targets = 0),
                                               list(width = "125px", targets = 1),
                                               list(width = "50px", targets = 2))),
              rownames = FALSE
            ) %>%
              formatSignif(columns = c("Estimated value", "Rhat"), digits = 4) %>%
              formatStyle(names(dt_diag), lineHeight = "80%")
          }
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
          rv_deconvo_result$traceplot
        })
      }
    }
  )
}
