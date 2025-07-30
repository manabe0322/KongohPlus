#' hd_true_importance_ui
#'
#' @description The function to create the ui module for the Hd true test
hd_true_importance_ui <- function(id){
  ns <- NS(id)
  hd_true_importance_ui_sub(ns("hd_true_importance_sub"))
}

#' hd_true_importance_ui_sub
#'
#' @description The function to create the sub-ui module for the Hd true test
hd_true_importance_ui_sub <- function(id){
  ns <- NS(id)

  tabPanel("Hd true test (importance sampling)",
           useWaiter(),
           titlePanel("Hd true test (importance sampling)"),
           fluidRow(
             column(3,
                    wellPanel(
                      h4("Crime stain profile"),
                      textOutput(ns("fn_csp")),
                      br(),
                      h4("Allele frequency"),
                      textOutput(ns("fn_af"))
                    ),
                    wellPanel(
                      uiOutput(ns("noc")),
                      br(),
                      numericInput(ns("n_simu"), label = h4("Number of simulations"), value = 100),
                      br(),
                      br(),
                      actionButton(ns("act_analyze_hd_true"), label = "Analyze", class = "btn btn-primary btn-lg"),
                      br(),
                      textOutput(ns("act_reason"))
                    )
             ),
             column(9,
                    textOutput(ns("expect_lr_hd")),
                    br(),
                    plotOutput(ns("hist_lr_importance"), height = "600px")
             )
           )
  )
}

#' hd_true_importance_server
#'
#' @description The function to create the server module for the Hd true test
hd_true_importance_server <- function(id, path_pack, rv_deconvo){
  moduleServer(
    id,
    function(input, output, session){
      rv_hd_true_importance <- reactiveValues()
      rv_hd_true_importance$lr_importance <- NULL
      rv_hd_true_importance$bias <- NULL
      rv_hd_true_importance$expect_lr_hd <- NULL

      observe({
        req(rv_deconvo)
        req(rv_hd_true_importance)
        rv_hd_true_importance <- hd_true_importance_server_sub("hd_true_importance_sub", path_pack, rv_deconvo, rv_hd_true_importance)
      })

      return(rv_hd_true_importance)
    }
  )
}

#' hd_true_importance_server_sub
#'
#' @description The function to create the sub-server module for the Hd true test
hd_true_importance_server_sub <- function(id, path_pack, rv_deconvo, rv_hd_true_importance){
  moduleServer(
    id,
    function(input, output, session){
      data_deconvo <- rv_deconvo$data_deconvo

      if(is.null(data_deconvo)){
        disable("act_analyze_hd_true")
      }else{
        hnc_all <- data_deconvo$hnc_all
        rv_hd_true_importance$theta <- data_deconvo$theta

        if(rv_hd_true_importance$theta == 0){
          enable("act_analyze_hd_true")
        }else{
          disable("act_analyze_hd_true")
        }

        output$act_reason <- renderText({
          if(rv_hd_true_importance$theta == 0){
            ""
          }else{
            "This feature is available only when theta is set to 0."
          }
        })

        output$fn_csp <- renderText({paste0(data_deconvo$fn_csp)})
        output$fn_af <- renderText({paste0(data_deconvo$fn_af)})

        output$noc <- renderUI({
          selectInput(session$ns("noc"), label = "Number of contributors", choices = hnc_all, selected = hnc_all[1], selectize = FALSE)
        })

        observeEvent(input$act_analyze_hd_true, {
          waiter_show(html = spin_3k(), color = "white")

          hnc <- as.numeric(input$noc)
          n_simu <- input$n_simu
          tmp <- perform_hd_true(hnc, n_simu, data_deconvo, importance_sampling = TRUE)
          rv_hd_true_importance$lr_importance <- lr_importance <- tmp$lr_rand
          rv_hd_true_importance$bias <- bias <- tmp$bias

          if(any(is.na(lr_importance))){
            data_hnc_list <- data_deconvo$data_hnc_list
            data_h1 <- data_hnc_list[[which(hnc_all == hnc)]]
            error_message_h1 <- data_h1$error_message
            waiter_hide()
            showModal(modalDialog(title = "Error", HTML(error_message_h1), easyClose = TRUE, footer = NULL))
          }else{
            rv_hd_true_importance$expect_lr_hd <- sum(lr_importance * bias) / n_simu
            waiter_hide()
          }
        })

        output$expect_lr_hd <- renderText({
          expect_lr_hd <- rv_hd_true_importance$expect_lr_hd
          if(!is.null(expect_lr_hd)){
            paste0(expect_lr_hd)
          }
        })

        output$hist_lr_importance <- renderPlot({
          lr_importance <- rv_hd_true_importance$lr_importance
          if(!is.null(lr_importance)){
            log10_lr_importance <- log10(lr_importance)
            log10_lr_importance[log10_lr_importance < -100] <- -100
            hist(log10_lr_importance)
          }
        })
      }

      return(rv_hd_true_importance)
    }
  )
}
