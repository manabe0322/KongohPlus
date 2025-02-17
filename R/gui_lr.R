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
                      checkboxInput(ns("sub_source"), label = "Sub-source LR", value = TRUE),
                      br(),
                      br(),
                      actionButton(ns("act_analyze_lr"), label = "Analyze", class = "btn btn-primary btn-lg")
                    )
             ),
             column(3,
                    wellPanel(
                      h4(textOutput(ns("title_lr"))),
                      dataTableOutput(ns("dt_lr")),
                      br()
                    )
             ),
             column(6)
           )
  )
}

#' lr_server
#'
#' @description The function to create the server module for the likelihood ratio
lr_server <- function(id, rv_file, rv_deconvo){
  moduleServer(
    id,
    function(input, output, session){
      rv_lr_setting <- reactiveValues()
      rv_lr_setting$dt_ref <- NULL
      rv_lr_setting$fn_ref <- NULL
      rv_lr_setting$refs <- NULL

      rv_lr <- reactiveValues()
      rv_lr$dt_ref <- NULL
      rv_lr$fn_ref <- NULL
      rv_lr$refs <- NULL
      rv_lr$dt_lr <- NULL

      observe({
        req(rv_file)
        rv_lr_setting$dt_ref <- dt_ref <- rv_file$dt_ref
        rv_lr_setting$fn_ref <- rv_file$fn_ref
        if(!is.null(dt_ref)){
          rv_lr_setting$refs <- unique(dt_ref[, SampleName])
        }
      })

      observe({
        req(rv_deconvo)
        data_deconvo <- rv_deconvo$data_deconvo

        if(!is.null(data_deconvo)){
          hnc_all <- data_deconvo$hnc_all

          output$hp <- renderUI({
            checkboxGroupInput(session$ns("hp"), label = h4("Prosecutor hypothesis"), choices = rv_lr_setting$refs)
          })

          output$noc_hp <- renderUI({
            selectInput(session$ns("noc_hp"), label = "Number of contributors", choices = hnc_all, selected = hnc_all[1], selectize = FALSE)
          })

          output$hd <- renderUI({
            checkboxGroupInput(session$ns("hd"), label = h4("Defense hypothesis"), choices = rv_lr_setting$refs)
          })

          output$noc_hd <- renderUI({
            selectInput(session$ns("noc_hd"), label = "Number of contributors", choices = hnc_all, selected = hnc_all[1], selectize = FALSE)
          })

          observeEvent(input$act_analyze_lr, {
            waiter_show(html = spin_3k(), color = "white")

            # Fix setting
            hp <- input$hp
            noc_hp <- as.numeric(input$noc_hp)
            hd <- input$hd
            noc_hd <- as.numeric(input$noc_hd)
            refs <- rv_lr_setting$refs
            dt_ref <- rv_lr_setting$dt_ref
            fn_ref <- rv_lr_setting$fn_ref
            sub_source <- input$sub_source

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
              dt_anameth <- data_deconvo$dt_anameth
              theta <- as.numeric(dt_anameth$Value[dt_anameth$Subcategory == "Theta"])
              data_hnc_list <- data_deconvo$data_hnc_list

              # Extract data for Hp
              data_hp <- data_hnc_list[[which(hnc_all == noc_hp)]]
              error_message_hp <- data_hp$error_message
              gtc_use_hp <- data_hp$gtc_use
              lpd_gtc_use_1st_hp <- data_hp$lpd_gtc_use_1st
              lpd_gtc_use_2nd_hp <- data_hp$lpd_gtc_use_2nd

              # Extract data for Hd
              data_hd <- data_hnc_list[[which(hnc_all == noc_hd)]]
              error_message_hd <- data_hd$error_message
              gtc_use_hd <- data_hd$gtc_use
              lpd_gtc_use_1st_hd <- data_hd$lpd_gtc_use_1st
              lpd_gtc_use_2nd_hd <- data_hd$lpd_gtc_use_2nd

              if(is.null(error_message_hp) && is.null(error_message_hd)){
                n_ref <- length(refs)
                if(n_ref > noc_hp || n_ref > noc_hd){
                  n_row <- length(markers_use) + 1
                  dt_lr <- data.table(Marker = c(markers_use, "Total"),
                                      Log10_Likelihood_Hp = rep(NA, n_row),
                                      Log10_Likelihood_Hd = rep(NA, n_row),
                                      LR = rep(NA, n_row))
                }else{
                  log10like_hp <- assign_log10like(noc_hp, markers_use, peaks_list, pop_al_list, pop_prob_list, theta, gtc_use_hp, lpd_gtc_use_1st_hp, lpd_gtc_use_2nd_hp, dt_ref, hp, refs, sub_source)
                  log10like_hd <- assign_log10like(noc_hd, markers_use, peaks_list, pop_al_list, pop_prob_list, theta, gtc_use_hd, lpd_gtc_use_1st_hd, lpd_gtc_use_2nd_hd, dt_ref, hd, refs, sub_source)
                  dt_lr <- data.table(Marker = c(markers_use, "Total"),
                                      Log10_Likelihood_Hp = log10like_hp,
                                      Log10_Likelihood_Hd = log10like_hd,
                                      LR = 10^(log10like_hp - log10like_hd))
                }
              }else{
                n_row <- length(markers_use) + 1
                dt_lr <- data.table(Marker = c(markers_use, "Total"),
                                    Log10_Likelihood_Hp = rep(NA, n_row),
                                    Log10_Likelihood_Hd = rep(NA, n_row),
                                    LR = rep(NA, n_row))
              }
              names(dt_lr) <- c("Marker", "log10 (f(E | Hp))", "log10 (f(E | Hd))", "LR")

              rv_lr$dt_ref <- dt_ref
              rv_lr$fn_ref <- fn_ref
              rv_lr$refs <- refs
              rv_lr$hp <- hp
              rv_lr$hd <- hd
              rv_lr$dt_lr <- dt_lr

              waiter_hide()
            }
          })

          output$title_lr <- renderText({
            dt_lr <- rv_lr$dt_lr
            if(is.null(dt_lr)){
              ""
            }else{
              "Likelihood ratio"
            }
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

          output$dt_lr <- renderDataTable({
            dt_lr <- rv_lr$dt_lr
            if(!is.null(dt_lr)){
              datatable(
                dt_lr,
                selection = "none",
                options = list(dom = "t", iDisplayLength = nrow(dt_lr), ordering = FALSE, autowidth = TRUE, rowCallback = JS(rowCallback),
                               columnDefs = list(list(targets = c(1, 2), visible = FALSE))),
                rownames = FALSE
              ) %>%
                formatSignif(columns = c("LR"), digits = 3) %>%
                formatStyle(names(dt_lr), lineHeight = "80%")
            }
          })
        }
      })

      return(rv_lr)
    }
  )
}
