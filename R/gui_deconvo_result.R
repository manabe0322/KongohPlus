#' deconvo_result_ui
#'
#' @description The function to create the ui module for results of deconvolution
#' @param id The ui ID
deconvo_result_ui <- function(id){
  ns <- NS(id)

  tabPanel("Result",
           tags$style("
             .edit-modal .modal-dialog {
               width: 90% !important;   /* Expand up to 90% of the screen width */
               max-width: 2000px;       /* Set maximum width */
             }
             .edit-modal .modal-body {
               max-height: 80vh;        /* Limit height to 80% of the viewport */
               overflow-y: auto;        /* Enable vertical scrolling if needed */
             }
           "),
           br(),
           fluidRow(
             column(2,
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
                      fluidRow(
                        column(
                          3,
                          uiOutput(ns("select_noc_param"))
                        ),
                        column(
                          5,
                          br(),
                          uiOutput(ns("button_ave_ph"))
                        ),
                        column(
                          4,
                          br(),
                          uiOutput(ns("button_all_mle"))
                        )
                      ),
                      uiOutput(ns("param")),
                      br()
                    ),
                    wellPanel(
                      h4(textOutput(ns("title_warnings"))),
                      htmlOutput(ns("warning_message")),
                      br()
                    )
             )
           )
  )
}

#' deconvo_result_server
#'
#' @description The function to create the server module for results of deconvolution
#' @param id The server ID
#' @param session_top The top session
#' @param rv_deconvo The reactive values derived from the setting server for deconvolution
deconvo_result_server <- function(id, session_top, rv_deconvo){
  moduleServer(
    id,
    function(input, output, session){
      rv_deconvo_result <- reactiveValues()
      rv_deconvo_result$index_mle_best <- NULL
      rv_deconvo_result$index_converged <- NULL
      rv_deconvo_result$dt_all_mle <- NULL
      rv_deconvo_result$dt_param <- NULL
      rv_deconvo_result$error_message <- NULL
      rv_deconvo_result$markers_use <- NULL
      rv_deconvo_result$message_sr <- NULL
      rv_deconvo_result$hnc_all <- NULL
      rv_deconvo_result$data_hnc_list <- NULL
      rv_deconvo_result$min_size <- NULL
      rv_deconvo_result$max_size <- NULL
      rv_deconvo_result$base_size <- NULL

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

      observe({
        req(rv_deconvo$data_deconvo)

        data_deconvo <- rv_deconvo$data_deconvo
        dt_anameth <- data_deconvo$dt_anameth

        rv_deconvo_result$index_mle_best <- NULL
        rv_deconvo_result$index_converged <- NULL
        rv_deconvo_result$dt_all_mle <- NULL
        rv_deconvo_result$dt_param <- NULL
        rv_deconvo_result$error_message <- NULL
        rv_deconvo_result$markers_use <- data_deconvo$markers_use
        rv_deconvo_result$message_sr <- data_deconvo$message_sr
        rv_deconvo_result$hnc_all <- data_deconvo$hnc_all
        rv_deconvo_result$data_hnc_list <- data_deconvo$data_hnc_list
        rv_deconvo_result$min_size <- data_deconvo$min_size
        rv_deconvo_result$max_size <- data_deconvo$max_size
        rv_deconvo_result$base_size <- dt_anameth$Value[dt_anameth$Subcategory == "Base size"]
      })

      output$title_warnings <- renderText({
        req(rv_deconvo$data_deconvo)
        "Warnings"
      })

      output$title_like <- renderText({
        req(rv_deconvo$data_deconvo)
        "Likelihood"
      })

      output$title_pg <- renderText({
        req(rv_deconvo$data_deconvo)
        "Probabilistic genotyping"
      })

      output$title_param <- renderText({
        req(rv_deconvo$data_deconvo)
        "Parameters"
      })

      output$warning_message <- renderUI({
        data_hnc_list <- rv_deconvo_result$data_hnc_list
        message_sr <- rv_deconvo_result$message_sr

        warning_message <- ""

        if(!is.null(data_hnc_list)){
          for(i in 1:length(data_hnc_list)){
            data_hnc <- data_hnc_list[[i]]
            tmp <- data_hnc$warning_message
            if(!is.null(tmp)){
              warning_message <- paste0(warning_message, tmp)
            }
          }
        }

        warning_message <- paste0(warning_message, message_sr)
        HTML(warning_message)
      })

      output$act_lr <- renderUI({
        req(rv_deconvo$data_deconvo)
        actionButton(session$ns("act_lr"), label = "Likelihood ratio", class = "btn btn-primary btn-lg")
      })

      observeEvent(input$act_lr, {
        req(rv_deconvo$data_deconvo)
        updateNavbarPage(session = session_top, inputId = "navbar", selected = "Likelihood ratio")
      })

      output$dt_log10like <- renderDataTable(server = FALSE, {
        req(rv_deconvo_result$hnc_all)
        req(rv_deconvo_result$data_hnc_list)
        req(rv_deconvo_result$markers_use)

        hnc_all <- rv_deconvo_result$hnc_all
        n_hnc <- length(hnc_all)
        data_hnc_list <- rv_deconvo_result$data_hnc_list
        markers_use <- rv_deconvo_result$markers_use
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
        req(rv_deconvo_result$hnc_all)
        hnc_all <- rv_deconvo_result$hnc_all
        selectInput(session$ns("select_noc_param"), label = "Number of contributors", choices = hnc_all, selected = hnc_all[1], multiple = FALSE, selectize = FALSE)
      })

      output$button_ave_ph <- renderUI({
        req(rv_deconvo$data_deconvo)
        actionButton(session$ns("act_ave_ph"), label = "Illustrate the average peak height")
      })

      observeEvent(input$act_ave_ph, {
        showModal(modalDialog(
          title = "Average peak height",
          plotOutput(session$ns("plot_ave_ph"), height = "350px"),
          easyClose = TRUE,
          footer = NULL
        ))
      })

      output$button_all_mle <- renderUI({
        req(rv_deconvo$data_deconvo)
        actionButton(session$ns("act_all_mle"), label = "View all MLE attempts")
      })

      observeEvent(input$act_all_mle, {
        showModal(modalDialog(
          title = "All MLE attempts",

          div(
            span(style = "color:#BFEFFF;", "■"), " : Selected maximum likelihood",
            span(style = "color:#ADD8E6;", "■"), " : Close to the maximum likelihood"
          ),

          dataTableOutput(session$ns("dt_all_mle")),
          easyClose = TRUE,
          footer = NULL
        ) %>% tagAppendAttributes(class = "edit-modal")
        )
      })

      output$dt_all_mle <- renderDataTable({
        req(rv_deconvo_result$index_mle_best)
        req(rv_deconvo_result$index_converged)
        req(rv_deconvo_result$dt_all_mle)

        index_mle_best <- rv_deconvo_result$index_mle_best
        index_converged <- rv_deconvo_result$index_converged
        dt_all_mle <- rv_deconvo_result$dt_all_mle

        name_best <- paste0("MLE ", index_mle_best)
        name_converged <- setdiff(paste0("MLE ", index_converged), name_best)

        datatable(
          dt_all_mle,
          selection = "none",
          options = list(dom = "t", iDisplayLength = nrow(dt_all_mle), ordering = FALSE, autowidth = TRUE),
          rownames = FALSE
        ) %>%
          formatSignif(columns = paste0("MLE ", 1:(ncol(dt_all_mle) - 1)), digits = 6) %>%
          formatStyle(names(dt_all_mle), lineHeight = "80%") %>%
          formatStyle(name_best, backgroundColor = "#BFEFFF") %>%
          formatStyle(name_converged, backgroundColor = "#ADD8E6")
      })

      output$plot_ave_ph <- renderPlot({
        req(rv_deconvo_result$dt_param)
        req(rv_deconvo_result$min_size)
        req(rv_deconvo_result$max_size)
        req(rv_deconvo_result$base_size)

        create_fig_eph(rv_deconvo_result$dt_param, rv_deconvo_result$min_size, rv_deconvo_result$max_size, rv_deconvo_result$base_size)
      })

      observeEvent(input$select_noc_param, {
        req(rv_deconvo_result$hnc_all)
        req(rv_deconvo_result$data_hnc_list)

        hnc_all <- rv_deconvo_result$hnc_all
        data_hnc_list <- rv_deconvo_result$data_hnc_list

        select_noc_param <- as.integer(input$select_noc_param)
        index_hnc <- which(hnc_all == select_noc_param)

        if(length(index_hnc) != 0){
          data_hnc <- data_hnc_list[[index_hnc]]
          error_message <- data_hnc$error_message
          if(is.null(data_hnc$error_message)){
            rv_deconvo_result$index_mle_best <- data_hnc$index_mle_best
            rv_deconvo_result$index_converged <- data_hnc$index_converged
            rv_deconvo_result$dt_all_mle <- data_hnc$dt_all_mle
            rv_deconvo_result$dt_param <- data_hnc$dt_param
            rv_deconvo_result$error_message <- NULL
          }else{
            rv_deconvo_result$index_mle_best <- NULL
            rv_deconvo_result$index_converged <- NULL
            rv_deconvo_result$dt_all_mle <- NULL
            rv_deconvo_result$dt_param <- NULL
            rv_deconvo_result$error_message <- error_message
          }
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
        req(rv_deconvo_result$dt_param)
        dt_param <- rv_deconvo_result$dt_param

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
      })
    }
  )
}
