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
                      h4(textOutput(ns("title_param"))),
                      uiOutput(ns("select_noc_param")),
                      uiOutput(ns("param")),
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
        rv_deconvo_result$dt_param <- NULL
        rv_deconvo_result$error_message <- NULL

        output$title_warnings <- renderText({"Warnings"})
        output$title_like <- renderText({"Likelihood"})
        output$title_pg <- renderText({"Probabilistic genotyping"})
        output$title_param <- renderText({"Parameters"})

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
          dt_log10like_display <- data.table(Marker = c(markers_use, "Overall"))
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

        output$select_noc_param <- renderUI({
          selectInput(session$ns("select_noc_param"), label = "Number of contributors", choices = hnc_all, selected = hnc_all[1], multiple = FALSE, selectize = FALSE)
        })

        observeEvent(input$select_noc_param, {
          index_hnc <- which(hnc_all == input$select_noc_param)
          data_hnc <- data_hnc_list[[index_hnc]]
          error_message <- data_hnc$error_message
          if(is.null(data_hnc$error_message)){
            rv_deconvo_result$dt_param <- data_hnc$dt_param
            rv_deconvo_result$error_message <- NULL
          }else{
            rv_deconvo_result$dt_param <- NULL
            rv_deconvo_result$error_message <- error_message
          }
        })

        output$param <- renderUI({
          dt_param <- rv_deconvo_result$dt_param
          if(is.null(dt_param)){
            HTML(rv_deconvo_result$error_message)
          }else{
            dataTableOutput(session$ns("dt_param"))
          }
        })

        output$dt_param <- renderDataTable({
          dt_param <- rv_deconvo_result$dt_param

          if(!is.null(dt_param)){
            datatable(
              dt_param,
              selection = "none",
              options = list(dom = "t", iDisplayLength = nrow(dt_param), ordering = FALSE, autowidth = TRUE,
                             columnDefs = list(list(width = "150px", targets = 0),
                                               list(width = "70px", targets = 1),
                                               list(width = "50px", targets = 2),
                                               list(width = "70px", targets = 3))),
              rownames = FALSE
            ) %>%
              formatSignif(columns = c("Estimated value"), digits = 4) %>%
              formatStyle(names(dt_param), lineHeight = "80%")
          }
        })
      }
    }
  )
}
