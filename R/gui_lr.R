#' lr_ui
#'
#' @description The function to create the ui module for the likelihood ratio
lr_ui <- function(id){
  ns <- NS(id)

  tabPanel("Likelihood ratio",
           useWaiter(),
           fluidRow(
             column(3,
                    wellPanel(
                      uiOutput(ns("hp")),
                      uiOutput(ns("noc_hp")),
                      uiOutput(ns("hd")),
                      uiOutput(ns("noc_hd")),
                      br(),
                      br(),
                      uiOutput(ns("act_analyze_lr"))
                    )
             ),
             column(3,
                    wellPanel(
                      h4(textOutput(ns("title_lr"))),
                      dataTableOutput(ns("dt_lr")),
                      br()
                    )
             ),
             column(6,
                    wellPanel(
                      h4(textOutput(ns("title_lr_calc_record"))),
                      dataTableOutput(ns("dt_lr_calc_record"))
                    )
             )
           )
  )
}

#' create_contrib_order
#'
#' @description The function to create a vector of the ordered contributor
#' @param hyp_id_best A vector of the contributors' IDs for the best hypothesis
#' @param refs A vector of the references' names
create_contrib_order <- function(hyp_id_best, refs){
  hnc <- length(hyp_id_best)
  contrib_order <- rep("U", hnc)
  for(i in 1:hnc){
    id1 <- hyp_id_best[i]
    if(id1 != 0){
      contrib_order[i] <- refs[id1]
    }
  }
  return(contrib_order)
}

#' lr_server
#'
#' @description The function to create the server module for the likelihood ratio
lr_server <- function(id, rv_input, rv_deconvo){
  moduleServer(
    id,
    function(input, output, session){
      rv_lr <- reactiveValues()
      rv_lr$dt_ref <- NULL
      rv_lr$fn_ref <- NULL
      rv_lr$refs <- NULL
      rv_lr$lr_list <- NULL
      rv_lr$dt_lr_calc_record <- NULL
      rv_lr$dt_lr_display <- NULL

      observe({
        req(rv_input)
        req(rv_lr)

        if(!is.null(rv_input$dt_ref)){
          rv_lr$dt_ref <- dt_ref <- rv_input$dt_ref
          rv_lr$fn_ref <- rv_input$fn_ref
          if(!is.null(dt_ref)){
            rv_lr$refs <- unique(dt_ref[, SampleName])
          }
        }
      })

      observe({
        req(rv_deconvo)
        data_deconvo <- rv_deconvo$data_deconvo

        if(!is.null(data_deconvo)){
          hnc_all <- data_deconvo$hnc_all

          output$hp <- renderUI({
            checkboxGroupInput(session$ns("hp"), label = h4("Prosecutor hypothesis (Hp)"), choices = rv_lr$refs)
          })

          output$noc_hp <- renderUI({
            selectInput(session$ns("noc_hp"), label = "Number of contributors", choices = hnc_all, selected = hnc_all[1], selectize = FALSE)
          })

          output$hd <- renderUI({
            checkboxGroupInput(session$ns("hd"), label = h4("Defense hypothesis (Hd)"), choices = rv_lr$refs)
          })

          output$noc_hd <- renderUI({
            selectInput(session$ns("noc_hd"), label = "Number of contributors", choices = hnc_all, selected = hnc_all[1], selectize = FALSE)
          })

          output$"act_analyze_lr" <- renderUI({
            actionButton(session$ns("act_analyze_lr"), label = "Calculate", class = "btn btn-primary btn-lg")
          })

          observeEvent(input$act_analyze_lr, {
            waiter_show(html = spin_3k(), color = "white")

            # Fix setting
            hp <- input$hp
            noc_hp <- as.numeric(input$noc_hp)
            hd <- input$hd
            noc_hd <- as.numeric(input$noc_hd)
            refs <- rv_lr$refs
            dt_ref <- rv_lr$dt_ref

            dt_csp <- data_deconvo$dt_csp
            dt_af <- data_deconvo$dt_af
            dt_kit <- data_deconvo$dt_kit
            dt_kit_auto <- extract_kit_auto(dt_kit)
            markers_auto <- unique(dt_kit_auto[, Marker])
            error_message <- check_main_file(dt_csp, dt_ref, dt_af, markers_auto, method = "lr")

            if(error_message != ""){
              waiter_hide()
              showModal(modalDialog(title = "Error", HTML(error_message), easyClose = TRUE, footer = NULL))
            }else{
              # Extract data from data_deconvo
              markers_use <- data_deconvo$markers_use
              peaks_list <- data_deconvo$peaks_list
              pop_al_list <- data_deconvo$pop_al_list
              pop_prob_list <- data_deconvo$pop_prob_list
              theta <- data_deconvo$theta
              data_hnc_list <- data_deconvo$data_hnc_list

              # Extract data for Hp
              data_hp <- data_hnc_list[[which(hnc_all == noc_hp)]]
              error_message_hp <- data_hp$error_message
              dt_pg_list_hp <- data_hp$dt_pg_list

              # Extract data for Hd
              data_hd <- data_hnc_list[[which(hnc_all == noc_hd)]]
              error_message_hd <- data_hd$error_message
              dt_pg_list_hd <- data_hd$dt_pg_list

              if(length(hp) > noc_hp || length(hd) > noc_hd){
                waiter_hide()
                showModal(modalDialog(title = "Error", "The number of contributors is smaller than the number of known individuals who have been checked as contributors.", easyClose = TRUE, footer = NULL))
              }else if(!is.null(error_message_hp)){
                waiter_hide()
                showModal(modalDialog(title = "Error", HTML(error_message_hp), easyClose = TRUE, footer = NULL))
              }else if(!is.null(error_message_hd)){
                waiter_hide()
                showModal(modalDialog(title = "Error", HTML(error_message_hd), easyClose = TRUE, footer = NULL))
              }else{
                tmp <- calculate_lr(noc_hp, hp, noc_hd, hd, refs,
                                    markers_use, peaks_list, pop_al_list, pop_prob_list, theta, dt_pg_list_hp, dt_pg_list_hd, dt_ref)

                lr_data <- tmp$lr_data
                rv_lr$dt_lr_display <- lr_data$dt_lr

                lr_list <- rv_lr$lr_list
                if(is.null(lr_list)){
                  lr_list <- list()
                  lr_list[[1]] <- lr_data
                }else{
                  lr_list[[length(lr_list) + 1]] <- lr_data
                }
                rv_lr$lr_list <- lr_list

                rv_lr$dt_lr_calc_record <- rbind(rv_lr$dt_lr_calc_record, tmp$dt_lr_calc_record_1)

                waiter_hide()
              }
            }
          })

          output$title_lr <- renderText({
            dt_lr_calc_record <- rv_lr$dt_lr_calc_record
            if(is.null(dt_lr_calc_record)){
              ""
            }else{
              "Likelihood ratio"
            }
          })

          output$dt_lr <- renderDataTable({
            dt_lr <- rv_lr$dt_lr_display
            if(!is.null(dt_lr)){
              datatable(
                dt_lr,
                selection = "none",
                options = list(dom = "t", iDisplayLength = nrow(dt_lr), ordering = FALSE, autowidth = TRUE,
                               columnDefs = list(list(targets = c(1, 2), visible = FALSE))),
                rownames = FALSE
              ) %>%
                formatSignif(columns = c("LR"), digits = 3) %>%
                formatStyle(names(dt_lr), lineHeight = "80%")
            }
          })

          output$title_lr_calc_record <- renderText({
            dt_lr_calc_record <- rv_lr$dt_lr_calc_record
            if(is.null(dt_lr_calc_record)){
              ""
            }else{
              "Calculation record"
            }
          })

          output$dt_lr_calc_record <- renderDataTable({
            dt_lr_calc_record <- rv_lr$dt_lr_calc_record
            if(!is.null(dt_lr_calc_record)){
              datatable(
                dt_lr_calc_record,
                selection = list(mode = "single", target = "row"),
                options = list(dom = "t", iDisplayLength = nrow(dt_lr_calc_record), ordering = FALSE, autowidth = TRUE),
                rownames = FALSE
              ) %>%
                formatSignif(columns = c("LR"), digits = 3) %>%
                formatStyle(names(dt_lr_calc_record), lineHeight = "80%")
            }
          })

          observeEvent(input$dt_lr_calc_record_rows_selected, {
            index_select <- input$dt_lr_calc_record_rows_selected

            if(!is.null(index_select)){
              lr_tmp <- rv_lr$lr_list[[index_select]]
              rv_lr$dt_lr_display <- lr_tmp$dt_lr
            }
          })
        }
      })

      return(rv_lr)
    }
  )
}
