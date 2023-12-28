est_par_result_ui <- function(id){
  ns <- NS(id)

  tabPanel("Result",
           br(),
           sidebarLayout(
             sidebarPanel(width = 2,
                          disabled(actionButton(ns("act_traceplot"), label = "Trace plot")),
                          br(),
                          br(),
                          br(),
                          disabled(actionButton(ns("act_dataplot"), label = "Data plot")),
                          br(),
                          br(),
                          br(),
                          disabled(actionButton(ns("act_alcor"), label = "Corrected allele")),
                          br(),
                          br(),
                          br(),
                          disabled(actionButton(ns("act_aic_model"), label = "AIC")),
                          br(),
                          br(),
                          br(),
                          disabled(actionButton(ns("act_save_par"), label = "Save", class = "btn btn-primary btn-lg"))
             ),
             mainPanel(width = 10,
                       tabsetPanel(id = ns("tabset_est_par_result"),
                                   type = "tabs",
                                   tabPanel("AE",
                                            br(),
                                            dataTableOutput(ns("dt_result_ae")),
                                            br()
                                   ),
                                   tabPanel("BSR",
                                            br(),
                                            dataTableOutput(ns("dt_result_bsr")),
                                            br()
                                   ),
                                   tabPanel("FSR",
                                            br(),
                                            dataTableOutput(ns("dt_result_fsr")),
                                            br()
                                   ),
                                   tabPanel("DSR",
                                            br(),
                                            dataTableOutput(ns("dt_result_dsr")),
                                            br()
                                   ),
                                   tabPanel("M2SR",
                                            br(),
                                            dataTableOutput(ns("dt_result_m2sr")),
                                            br()
                                   ),
                                   tabPanel("Gamma",
                                            br(),
                                            dataTableOutput(ns("dt_result_gamma")),
                                            br()
                                   )
                       )
             )
           )
  )
}

est_par_result_server <- function(id, path_pack, rv_output){
  moduleServer(
    id,
    function(input, output, session){
      results_est_par <- rv_output$results_est_par

      if(!is.null(results_est_par)){
        # Extract result lists
        markers_par <- results_est_par$markers_par
        ae_result_list <- results_est_par$ae_result_list
        bsr_result_list <- results_est_par$bsr_result_list
        fsr_result_list <- results_est_par$fsr_result_list
        dsr_result_list <- results_est_par$dsr_result_list
        m2sr_result_list <- results_est_par$m2sr_result_list
        gamma_result_list <- results_est_par$gamma_result_list

        # Extract the combined data.table for corrected alleles
        dt_alcor_combined <- results_est_par$dt_alcor_combined

        # Extract the results of estimated parameters
        dt_par <- results_est_par$dt_par

        # Extract traceplot
        traceplot_ae <- ae_result_list$traceplot
        traceplot_bsr <- bsr_result_list$traceplot
        traceplot_fsr <- fsr_result_list$traceplot
        traceplot_dsr <- dsr_result_list$traceplot
        traceplot_m2sr <- m2sr_result_list$traceplot
        traceplot_gamma <- gamma_result_list$traceplot

        # Extract plot data
        plotdata_ae <- ae_result_list$plotdata
        plotdata_bsr <- bsr_result_list$plotdata
        plotdata_fsr <- fsr_result_list$plotdata
        plotdata_dsr <- dsr_result_list$plotdata
        plotdata_m2sr <- m2sr_result_list$plotdata

        # Extract information on the corrected alleles
        dt_alcor_bsr <- bsr_result_list$dt_alcor
        dt_alcor_fsr <- fsr_result_list$dt_alcor
        dt_alcor_dsr <- dsr_result_list$dt_alcor

        # Extract information on AIC values
        dt_aic_bsr <- bsr_result_list$dt_aic
        dt_aic_fsr <- fsr_result_list$dt_aic
        dt_aic_dsr <- dsr_result_list$dt_aic

        # Extract the table of results
        dt_result_ae <- ae_result_list$dt_result
        dt_result_bsr <- bsr_result_list$dt_result
        dt_result_fsr <- fsr_result_list$dt_result
        dt_result_dsr <- dsr_result_list$dt_result
        dt_result_m2sr <- m2sr_result_list$dt_result
        dt_result_gamma <- gamma_result_list$dt_result

        #### Change the status of each button ####

        enable("act_traceplot")

        enable("act_dataplot")

        # Change the status of the 'Corrected allele' button
        observeEvent(input$tabset_est_par_result, {
          if(is.element(input$tabset_est_par_result, c("AE", "M2SR", "Gamma"))){
            disable("act_alcor")
          }else{
            enable("act_alcor")
          }
        })

        # Change the status of the 'AIC' button
        observeEvent(input$tabset_est_par_result, {
          if(is.element(input$tabset_est_par_result, c("AE", "M2SR", "Gamma"))){
            disable("act_aic_model")
          }else{
            enable("act_aic_model")
          }
        })

        enable("act_save_par")

        #### Trace plot ####

        # Define reactive values for displaying the trace plot
        rv_traceplot <- reactiveValues()
        rv_traceplot$ae <- NULL
        rv_traceplot$bsr <- NULL
        rv_traceplot$fsr <- NULL
        rv_traceplot$dsr <- NULL
        rv_traceplot$m2sr <- NULL

        # Show the pop-up window to display the trace plot
        observeEvent(input$act_traceplot, {
          tabset_est_par_result <- input$tabset_est_par_result

          # AE
          if(tabset_est_par_result == "AE"){
            row_index_select <- input$dt_result_ae_rows_selected
            if(isTruthy(row_index_select)){
              rv_traceplot$ae <- traceplot_ae[[row_index_select]]

              showModal(modalDialog(
                title = "Trace plot",
                plotOutput(session$ns("traceplot_ae")),
                footer = tagList(
                  modalButton("Close")
                ),
                size = "l"
              ))

            }else{
              showModal(modalDialog(title = "Error", "Select a row!", easyClose = TRUE, footer = NULL))
            }
          }

          # BSR
          else if(tabset_est_par_result == "BSR"){
            row_index_select <- input$dt_result_bsr_rows_selected
            if(isTruthy(row_index_select)){
              if(dt_result_bsr$Model[row_index_select] != "Not consider"){
                rv_traceplot$bsr <- traceplot_bsr[[row_index_select]]

                showModal(modalDialog(
                  title = "Trace plot",
                  plotOutput(session$ns("traceplot_bsr")),
                  footer = tagList(
                    modalButton("Close")
                  ),
                  size = "l"
                ))
              }else{
                showModal(modalDialog(title = "Information", "The back stutter ratio is not considered in the selected marker!", easyClose = TRUE, footer = NULL))
              }
            }else{
              showModal(modalDialog(title = "Error", "Select a row!", easyClose = TRUE, footer = NULL))
            }
          }

          # FSR
          else if(tabset_est_par_result == "FSR"){
            row_index_select <- input$dt_result_fsr_rows_selected
            if(isTruthy(row_index_select)){
              if(dt_result_fsr$Model[row_index_select] != "Not consider"){
                rv_traceplot$fsr <- traceplot_fsr[[row_index_select]]

                showModal(modalDialog(
                  title = "Trace plot",
                  plotOutput(session$ns("traceplot_fsr")),
                  footer = tagList(
                    modalButton("Close")
                  ),
                  size = "l"
                ))
              }else{
                showModal(modalDialog(title = "Information", "The forward stutter ratio is not considered in the selected marker!", easyClose = TRUE, footer = NULL))
              }
            }else{
              showModal(modalDialog(title = "Error", "Select a row!", easyClose = TRUE, footer = NULL))
            }
          }

          # DSR
          else if(tabset_est_par_result == "DSR"){
            row_index_select <- input$dt_result_dsr_rows_selected
            if(isTruthy(row_index_select)){
              if(dt_result_dsr$Model[row_index_select] != "Not consider"){
                rv_traceplot$dsr <- traceplot_dsr[[row_index_select]]

                showModal(modalDialog(
                  title = "Trace plot",
                  plotOutput(session$ns("traceplot_dsr")),
                  footer = tagList(
                    modalButton("Close")
                  ),
                  size = "l"
                ))
              }else{
                showModal(modalDialog(title = "Information", "The double-back stutter ratio is not considered in the selected marker!", easyClose = TRUE, footer = NULL))
              }
            }else{
              showModal(modalDialog(title = "Error", "Select a row!", easyClose = TRUE, footer = NULL))
            }
          }

          # M2SR
          else if(tabset_est_par_result == "M2SR"){
            row_index_select <- input$dt_result_m2sr_rows_selected
            if(isTruthy(row_index_select)){
              if(dt_result_m2sr$Model[row_index_select] != "Not consider"){
                rv_traceplot$m2sr <- traceplot_m2sr[[row_index_select]]

                showModal(modalDialog(
                  title = "Trace plot",
                  plotOutput(session$ns("traceplot_m2sr")),
                  footer = tagList(
                    modalButton("Close")
                  ),
                  size = "l"
                ))
              }else{
                showModal(modalDialog(title = "Information", "The minus 2-nt stutter ratio is not considered in the selected marker!", easyClose = TRUE, footer = NULL))
              }
            }else{
              showModal(modalDialog(title = "Error", "Select a row!", easyClose = TRUE, footer = NULL))
            }
          }

          # Gamma
          else if(tabset_est_par_result == "Gamma"){
            showModal(modalDialog(
              title = "Trace plot",
              plotOutput(session$ns("traceplot_gamma")),
              footer = tagList(
                modalButton("Close")
              ),
              size = "l"
            ))
          }
        })

        # Display the trace plot for AE
        output$traceplot_ae <- renderPlot({
          rv_traceplot$ae
        })

        # Display the trace plot for BSR
        output$traceplot_bsr <- renderPlot({
          rv_traceplot$bsr
        })

        # Display the trace plot for FSR
        output$traceplot_fsr <- renderPlot({
          rv_traceplot$fsr
        })

        # Display the trace plot for DSR
        output$traceplot_dsr <- renderPlot({
          rv_traceplot$dsr
        })

        # Display the trace plot for M2SR
        output$traceplot_m2sr <- renderPlot({
          rv_traceplot$m2sr
        })

        # Display the trace plot for the gamma parameter
        output$traceplot_gamma <- renderPlot({
          traceplot_gamma
        })

        #### Data plot ####

        # Define reactive values for displaying the data plot
        rv_dataplot <- reactiveValues()
        rv_dataplot$bsr_xval <- NULL
        rv_dataplot$bsr_yval <- NULL
        rv_dataplot$fsr_xval <- NULL
        rv_dataplot$fsr_yval <- NULL
        rv_dataplot$dsr_xval <- NULL
        rv_dataplot$dsr_yval <- NULL
        rv_dataplot$m2sr_xval <- NULL
        rv_dataplot$m2sr_yval <- NULL

        # Show the pop-up window to display the data plot
        observeEvent(input$act_dataplot, {
          tabset_est_par_result <- input$tabset_est_par_result

          # AE
          if(tabset_est_par_result == "AE"){
            showModal(modalDialog(
              title = "Data plot",
              plotOutput(session$ns("dataplot_ae")),
              footer = tagList(
                modalButton("Close")
              ),
              size = "l"
            ))
          }

          # BSR
          else if(tabset_est_par_result == "BSR"){
            row_index_select <- input$dt_result_bsr_rows_selected
            if(isTruthy(row_index_select)){
              adopt_model_bsr <- dt_result_bsr$Model[row_index_select]
              if(adopt_model_bsr != "Not consider"){
                bsr_data_selected <- plotdata_bsr[[row_index_select]]
                if(is.element(adopt_model_bsr, c("LUS", "Multi-sequence"))){
                  rv_dataplot$bsr_xval <- bsr_data_selected$sr_ca
                }else{
                  rv_dataplot$bsr_xval <- bsr_data_selected$sr_al
                }
                rv_dataplot$bsr_yval <- bsr_data_selected$sr

                showModal(modalDialog(
                  title = "Data plot",
                  plotOutput(session$ns("dataplot_bsr")),
                  footer = tagList(
                    modalButton("Close")
                  ),
                  size = "l"
                ))
              }else{
                showModal(modalDialog(title = "Information", "The back stutter ratio is not considered in the selected marker!", easyClose = TRUE, footer = NULL))
              }
            }else{
              showModal(modalDialog(title = "Error", "Select a row!", easyClose = TRUE, footer = NULL))
            }
          }

          # FSR
          else if(tabset_est_par_result == "FSR"){
            row_index_select <- input$dt_result_fsr_rows_selected
            if(isTruthy(row_index_select)){
              adopt_model_fsr <- dt_result_fsr$Model[row_index_select]
              if(adopt_model_fsr != "Not consider"){
                fsr_data_selected <- plotdata_fsr[[row_index_select]]
                if(is.element(adopt_model_fsr, c("LUS", "Multi-sequence"))){
                  rv_dataplot$fsr_xval <- fsr_data_selected$sr_ca
                }else{
                  rv_dataplot$fsr_xval <- fsr_data_selected$sr_al
                }
                rv_dataplot$fsr_yval <- fsr_data_selected$sr

                showModal(modalDialog(
                  title = "Data plot",
                  plotOutput(session$ns("dataplot_fsr")),
                  footer = tagList(
                    modalButton("Close")
                  ),
                  size = "l"
                ))
              }else{
                showModal(modalDialog(title = "Information", "The back stutter ratio is not considered in the selected marker!", easyClose = TRUE, footer = NULL))
              }
            }else{
              showModal(modalDialog(title = "Error", "Select a row!", easyClose = TRUE, footer = NULL))
            }
          }

          # DSR
          else if(tabset_est_par_result == "DSR"){
            row_index_select <- input$dt_result_dsr_rows_selected
            if(isTruthy(row_index_select)){
              adopt_model_dsr <- dt_result_dsr$Model[row_index_select]
              if(adopt_model_dsr != "Not consider"){
                dsr_data_selected <- plotdata_dsr[[row_index_select]]
                if(is.element(adopt_model_dsr, c("LUS", "Multi-sequence"))){
                  rv_dataplot$dsr_xval <- dsr_data_selected$sr_ca
                }else{
                  rv_dataplot$dsr_xval <- dsr_data_selected$sr_al
                }
                rv_dataplot$dsr_yval <- dsr_data_selected$sr

                showModal(modalDialog(
                  title = "Data plot",
                  plotOutput(session$ns("dataplot_dsr")),
                  footer = tagList(
                    modalButton("Close")
                  ),
                  size = "l"
                ))
              }else{
                showModal(modalDialog(title = "Information", "The back stutter ratio is not considered in the selected marker!", easyClose = TRUE, footer = NULL))
              }
            }else{
              showModal(modalDialog(title = "Error", "Select a row!", easyClose = TRUE, footer = NULL))
            }
          }

          # M2SR
          else if(tabset_est_par_result == "M2SR"){
            row_index_select <- input$dt_result_m2sr_rows_selected
            if(isTruthy(row_index_select)){
              adopt_model_m2sr <- dt_result_m2sr$Model[row_index_select]
              if(adopt_model_m2sr != "Not consider"){
                m2sr_data_selected <- plotdata_m2sr[[row_index_select]]
                rv_dataplot$m2sr_xval <- m2sr_data_selected$sr_al
                rv_dataplot$m2sr_yval <- m2sr_data_selected$sr

                showModal(modalDialog(
                  title = "Data plot",
                  plotOutput(session$ns("dataplot_m2sr")),
                  footer = tagList(
                    modalButton("Close")
                  ),
                  size = "l"
                ))
              }else{
                showModal(modalDialog(title = "Information", "The back stutter ratio is not considered in the selected marker!", easyClose = TRUE, footer = NULL))
              }
            }else{
              showModal(modalDialog(title = "Error", "Select a row!", easyClose = TRUE, footer = NULL))
            }
          }
        })

        # Display the data plot for AE
        output$dataplot_ae <- renderPlot({
          boxplot(plotdata_ae)
        })

        # Display the data plot for BSR
        output$dataplot_bsr <- renderPlot({
          plot(rv_dataplot$bsr_xval, rv_dataplot$bsr_yval)
        })

        # Display the data plot for FSR
        output$dataplot_fsr <- renderPlot({
          plot(rv_dataplot$fsr_xval, rv_dataplot$fsr_yval)
        })

        # Display the data plot for DSR
        output$dataplot_dsr <- renderPlot({
          plot(rv_dataplot$dsr_xval, rv_dataplot$dsr_yval)
        })

        # Display the data plot for M2SR
        output$dataplot_m2sr <- renderPlot({
          plot(rv_dataplot$m2sr_xval, rv_dataplot$m2sr_yval)
        })

        #### Show allele repeat corrections ####

        rv_dt_alcor <- reactiveVal(NULL)

        observeEvent(input$act_alcor, {
          bool_make_dt <- FALSE
          tabset_est_par_result <- input$tabset_est_par_result

          # BSR
          if(tabset_est_par_result == "BSR"){
            row_index_select <- input$dt_result_bsr_rows_selected
            if(isTruthy(row_index_select)){
              adopt_model_bsr_all_mk <- dt_result_bsr[, Model]
              adopt_model_bsr <- adopt_model_bsr_all_mk[row_index_select]
              if(is.element(adopt_model_bsr, c("LUS", "Multi-sequence"))){
                bool_make_dt <- TRUE
                mk <- markers_par[row_index_select]
                setkey(dt_alcor_bsr, Marker)
                dt_alcor_bsr_mk <- dt_alcor_bsr[mk]
                alleles <- dt_alcor_bsr_mk[, Allele]
                alcors <- dt_alcor_bsr_mk[, Corrected]
              }else{
                showModal(modalDialog(title = "Information", "Information on corrected alleles are shown only when the adopted model is 'LUS' or 'Multi-sequence'.", easyClose = TRUE, footer = NULL))
              }
            }else{
              showModal(modalDialog(title = "Error", "Select a row!", easyClose = TRUE, footer = NULL))
            }
          }

          # FSR
          else if(tabset_est_par_result == "FSR"){
            row_index_select <- input$dt_result_fsr_rows_selected
            if(isTruthy(row_index_select)){
              adopt_model_fsr_all_mk <- dt_result_fsr[, Model]
              adopt_model_fsr <- adopt_model_fsr_all_mk[row_index_select]
              if(is.element(adopt_model_fsr, c("LUS", "Multi-sequence"))){
                bool_make_dt <- TRUE
                mk <- markers_par[row_index_select]
                setkey(dt_alcor_fsr, Marker)
                dt_alcor_fsr_mk <- dt_alcor_fsr[mk]
                alleles <- dt_alcor_fsr_mk[, Allele]
                alcors <- dt_alcor_fsr_mk[, Corrected]
              }else{
                showModal(modalDialog(title = "Information", "Information on corrected alleles are shown only when the adopted model is 'LUS' or 'Multi-sequence'.", easyClose = TRUE, footer = NULL))
              }
            }else{
              showModal(modalDialog(title = "Error", "Select a row!", easyClose = TRUE, footer = NULL))
            }
          }

          # DSR
          else if(tabset_est_par_result == "DSR"){
            row_index_select <- input$dt_result_dsr_rows_selected
            if(isTruthy(row_index_select)){
              adopt_model_dsr_all_mk <- dt_result_dsr[, Model]
              adopt_model_dsr <- adopt_model_dsr_all_mk[row_index_select]
              if(is.element(adopt_model_dsr, c("LUS", "Multi-sequence"))){
                bool_make_dt <- TRUE
                mk <- markers_par[row_index_select]
                setkey(dt_alcor_dsr, Marker)
                dt_alcor_dsr_mk <- dt_alcor_dsr[mk]
                alleles <- dt_alcor_dsr_mk[, Allele]
                alcors <- dt_alcor_dsr_mk[, Corrected]
              }else{
                showModal(modalDialog(title = "Information", "Information on corrected alleles are shown only when the adopted model is 'LUS' or 'Multi-sequence'.", easyClose = TRUE, footer = NULL))
              }
            }else{
              showModal(modalDialog(title = "Error", "Select a row!", easyClose = TRUE, footer = NULL))
            }
          }

          if(bool_make_dt){
            dt_alcor <- data.table(Allele = alleles,
                                   Corrected = signif(alcors, 3))
            rv_dt_alcor(dt_alcor)

            showModal(modalDialog(
              title = "Corrected allele",
              dataTableOutput(session$ns("dt_alcor")),
              footer = tagList(
                modalButton("Close")
              )
            ))
          }else{
            rv_dt_alcor(NULL)
          }
        })

        output$dt_alcor <- renderDataTable(server = FALSE, {
          dt_alcor <- rv_dt_alcor()
          datatable(
            dt_alcor,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_alcor), ordering = FALSE, autowidth = TRUE),
            rownames = FALSE
          ) %>%
            formatRound(columns = c("Corrected"), digits = 3)
        })

        #### Show AIC values and parameters ####

        # Define reactive values for displaying AIC values and parameters of each candidate model
        rv_dt_aic <- reactiveValues()
        rv_dt_aic$bsr <- NULL
        rv_dt_aic$fsr <- NULL
        rv_dt_aic$dsr <- NULL

        # Show the pop-up window to show AIC values and parameters of candidate models
        observeEvent(input$act_aic_model, {
          tabset_est_par_result <- input$tabset_est_par_result

          # BSR
          if(tabset_est_par_result == "BSR"){
            row_index_select <- input$dt_result_bsr_rows_selected
            if(isTruthy(row_index_select)){
              adopt_model_bsr <- dt_result_bsr$Model[row_index_select]
              mk <- markers_par[row_index_select]
              markers_aic <- dt_aic_bsr[, Marker]
              if(is.element(mk, markers_aic)){
                setkey(dt_aic_bsr, Marker)
                dt_aic_bsr_mk <- dt_aic_bsr[mk]
                dt_aic_bsr_mk[, Marker := NULL]
                rv_dt_aic$bsr <- dt_aic_bsr_mk

                showModal(modalDialog(
                  title = "AIC values and parameters of candidate models",
                  dataTableOutput(session$ns("dt_aic_bsr")),
                  br(),
                  h5(div(paste0("Adopted model: ", adopt_model_bsr), style = "color:#555555;font-weight:bold;")),
                  footer = tagList(
                    modalButton("Close")
                  ),
                  size = "l"
                ))
              }else if(adopt_model_bsr == "Not consider"){
                showModal(modalDialog(title = "Information", "The back stutter ratio is not considered in the selected marker.", easyClose = TRUE, footer = NULL))
              }else{
                showModal(modalDialog(title = "Information", "Information on AIC values are shown only when 'Best' is selected in the setting tab.", easyClose = TRUE, footer = NULL))
              }
            }else{
              showModal(modalDialog(title = "Error", "Select a row!", easyClose = TRUE, footer = NULL))
            }
          }

          # FSR
          else if(tabset_est_par_result == "FSR"){
            row_index_select <- input$dt_result_fsr_rows_selected
            if(isTruthy(row_index_select)){
              adopt_model_fsr <- dt_result_fsr$Model[row_index_select]
              mk <- markers_par[row_index_select]
              markers_aic <- dt_aic_fsr[, Marker]
              if(is.element(mk, markers_aic)){
                setkey(dt_aic_fsr, Marker)
                dt_aic_fsr_mk <- dt_aic_fsr[mk]
                dt_aic_fsr_mk[, Marker := NULL]
                rv_dt_aic$fsr <- dt_aic_fsr_mk

                showModal(modalDialog(
                  title = "AIC values and parameters of candidate models",
                  dataTableOutput(session$ns("dt_aic_fsr")),
                  br(),
                  h5(div(paste0("Adopted model: ", adopt_model_fsr), style = "color:#555555;font-weight:bold;")),
                  footer = tagList(
                    modalButton("Close")
                  ),
                  size = "l"
                ))
              }else if(adopt_model_fsr == "Not consider"){
                showModal(modalDialog(title = "Information", "The forward stutter ratio is not considered in the selected marker.", easyClose = TRUE, footer = NULL))
              }else{
                showModal(modalDialog(title = "Information", "Information on AIC values are shown only when 'Best' is selected in the setting tab.", easyClose = TRUE, footer = NULL))
              }
            }else{
              showModal(modalDialog(title = "Error", "Select a row!", easyClose = TRUE, footer = NULL))
            }
          }

          # DSR
          else if(tabset_est_par_result == "DSR"){
            row_index_select <- input$dt_result_dsr_rows_selected
            if(isTruthy(row_index_select)){
              adopt_model_dsr <- dt_result_dsr$Model[row_index_select]
              mk <- markers_par[row_index_select]
              markers_aic <- dt_aic_dsr[, Marker]
              if(is.element(mk, markers_aic)){
                setkey(dt_aic_dsr, Marker)
                dt_aic_dsr_mk <- dt_aic_dsr[mk]
                dt_aic_dsr_mk[, Marker := NULL]
                rv_dt_aic$dsr <- dt_aic_dsr_mk

                showModal(modalDialog(
                  title = "AIC values and parameters of candidate models",
                  dataTableOutput(session$ns("dt_aic_dsr")),
                  br(),
                  h5(div(paste0("Adopted model: ", adopt_model_dsr), style = "color:#555555;font-weight:bold;")),
                  footer = tagList(
                    modalButton("Close")
                  ),
                  size = "l"
                ))
              }else if(adopt_model_dsr == "Not consider"){
                showModal(modalDialog(title = "Information", "The double-back stutter ratio is not considered in the selected marker.", easyClose = TRUE, footer = NULL))
              }else{
                showModal(modalDialog(title = "Information", "Information on AIC values are shown only when 'Best' is selected in the setting tab.", easyClose = TRUE, footer = NULL))
              }
            }else{
              showModal(modalDialog(title = "Error", "Select a row!", easyClose = TRUE, footer = NULL))
            }
          }
        })

        # Display information on AIC values and parameters of each back stutter model
        output$dt_aic_bsr <- renderDataTable(server = FALSE, {
          dt_aic_bsr <- rv_dt_aic$bsr
          datatable(
            dt_aic_bsr,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_aic_bsr), ordering = FALSE, autowidth = TRUE),
            rownames = FALSE
          ) %>%
            formatRound(columns = c("AIC", "Slope", "Intercept", "Rhat_Slope", "Rhat_Intercept", "Rhat_Sigma", "Rhat_Likelihood"), digits = 3) %>%
            formatRound(columns = c("Sigma", "MinUS"), digits = 0)
        })

        # Display information on AIC values and parameters of each forward stutter model
        output$dt_aic_fsr <- renderDataTable(server = FALSE, {
          dt_aic_fsr <- rv_dt_aic$fsr
          datatable(
            dt_aic_fsr,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_aic_fsr), ordering = FALSE, autowidth = TRUE),
            rownames = FALSE
          ) %>%
            formatRound(columns = c("AIC", "Slope", "Intercept", "Rhat_Slope", "Rhat_Intercept", "Rhat_Sigma", "Rhat_Likelihood"), digits = 3) %>%
            formatRound(columns = c("Sigma", "MinUS"), digits = 0)
        })

        # Display information on AIC values and parameters of each double-back stutter model
        output$dt_aic_dsr <- renderDataTable(server = FALSE, {
          dt_aic_dsr <- rv_dt_aic$dsr
          datatable(
            dt_aic_dsr,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_aic_dsr), ordering = FALSE, autowidth = TRUE),
            rownames = FALSE
          ) %>%
            formatRound(columns = c("AIC", "Slope", "Intercept", "Rhat_Slope", "Rhat_Intercept", "Rhat_Sigma", "Rhat_Likelihood"), digits = 3) %>%
            formatRound(columns = c("Sigma", "MinUS"), digits = 0)
        })

        #### Save parameters ####

        rv_par_names <- reactiveVal(NULL)

        observeEvent(input$act_save_par, {
          par_names <- list.files(paste0(path_pack, "/extdata/parameter"))
          rv_par_names(gsub(".{4}$", "", par_names))

          showModal(modalDialog(
            title = "Save parameters",
            textInput(session$ns("proj_par"), label = "Enter the project name"),
            footer = tagList(
              disabled(actionButton(session$ns("act_download_par"), "Save")),
              modalButton("Close")
            )
          ))
        })

        observeEvent(input$proj_par, {
          proj_par <- input$proj_par
          if(nchar(proj_par) == 0){
            showFeedbackDanger(inputId = "proj_par", text = "Required")
            disable("act_download_par")
          }else if(is.element(proj_par, rv_par_names())){
            showFeedbackDanger(inputId = "proj_par", text = "This name is already registered.")
            disable("act_download_par")
          }else{
            hideFeedback("proj_par")
            enable("act_download_par")
          }
        })

        observeEvent(input$act_download_par, {
          proj_par <- input$proj_par
          write.csv(dt_par, paste0(path_pack, "/extdata/parameter/", proj_par, ".csv"), row.names = FALSE)
          write.csv(dt_alcor_combined, paste0(path_pack, "/extdata/corrected_allele/alcor_", proj_par, ".csv"), row.names = FALSE)
          removeModal()
        })

        #### Result table ####

        # Display the result for AE
        output$dt_result_ae <- renderDataTable(server = FALSE, {
          datatable(
            dt_result_ae,
            selection = list(mode = "single", target = "row"),
            options = list(dom = "t", iDisplayLength = nrow(dt_result_ae), ordering = FALSE, autowidth = TRUE),
            rownames = FALSE
          ) %>%
            formatRound(columns = c("Myu", "Rhat_Myu", "Rhat_Sigma", "Rhat_Likelihood"), digits = 3) %>%
            formatRound(columns = c("Sigma"), digits = 0)
        })

        # Display the result for BSR
        output$dt_result_bsr <- renderDataTable(server = FALSE, {
          datatable(
            dt_result_bsr,
            selection = list(mode = "single", target = "row"),
            options = list(dom = "t", iDisplayLength = nrow(dt_result_bsr), ordering = FALSE, autowidth = TRUE),
            rownames = FALSE
          ) %>%
            formatRound(columns = c("Myu", "Slope", "Intercept", "Rhat_Myu", "Rhat_Slope", "Rhat_Intercept", "Rhat_Sigma", "Rhat_Likelihood"), digits = 3) %>%
            formatRound(columns = c("Sigma", "MinUS"), digits = 0)
        })

        # Display the result for FSR
        output$dt_result_fsr <- renderDataTable(server = FALSE, {
          datatable(
            dt_result_fsr,
            selection = list(mode = "single", target = "row"),
            options = list(dom = "t", iDisplayLength = nrow(dt_result_fsr), ordering = FALSE, autowidth = TRUE),
            rownames = FALSE
          ) %>%
            formatRound(columns = c("Myu", "Slope", "Intercept", "Rhat_Myu", "Rhat_Slope", "Rhat_Intercept", "Rhat_Sigma", "Rhat_Likelihood"), digits = 3) %>%
            formatRound(columns = c("Sigma", "MinUS"), digits = 0)
        })

        # Display the result for DSR
        output$dt_result_dsr <- renderDataTable(server = FALSE, {
          datatable(
            dt_result_dsr,
            selection = list(mode = "single", target = "row"),
            options = list(dom = "t", iDisplayLength = nrow(dt_result_dsr), ordering = FALSE, autowidth = TRUE),
            rownames = FALSE
          ) %>%
            formatRound(columns = c("Myu", "Slope", "Intercept", "Rhat_Myu", "Rhat_Slope", "Rhat_Intercept", "Rhat_Sigma", "Rhat_Likelihood"), digits = 3) %>%
            formatRound(columns = c("Sigma", "MinUS"), digits = 0)
        })

        # Display the result for M2SR
        output$dt_result_m2sr <- renderDataTable(server = FALSE, {
          datatable(
            dt_result_m2sr,
            selection = list(mode = "single", target = "row"),
            options = list(dom = "t", iDisplayLength = nrow(dt_result_m2sr), ordering = FALSE, autowidth = TRUE),
            rownames = FALSE
          ) %>%
            formatRound(columns = c("Myu", "Rhat_Myu", "Rhat_Sigma", "Rhat_Likelihood"), digits = 3) %>%
            formatRound(columns = c("Sigma"), digits = 0)
        })

        # Display the result for the gamma parameter
        output$dt_result_gamma <- renderDataTable(server = FALSE, {
          datatable(
            dt_result_gamma,
            selection = list(mode = "single", target = "row"),
            options = list(dom = "t", iDisplayLength = nrow(dt_result_gamma), ordering = FALSE, autowidth = TRUE),
            rownames = FALSE
          ) %>%
            formatRound(columns = c("Beta", "Rhat_Beta", "Rhat_Likelihood"), digits = 3)
        })
      }
    }
  )
}
