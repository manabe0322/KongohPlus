#' hd_true_ui
#'
#' @description The function to create the ui module for the Hd true test
hd_true_ui <- function(id){
  ns <- NS(id)
  hd_true_ui_sub(ns("hd_true_sub"))
}

#' hd_true_ui_sub
#'
#' @description The function to create the sub-ui module for the Hd true test
hd_true_ui_sub <- function(id){
  ns <- NS(id)

  tabPanel("Hd true test",
           useWaiter(),
           titlePanel("Hd true test"),
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
                      checkboxInput(ns("sub_source"), label = "Sub-source LR", value = TRUE),
                      br(),
                      br(),
                      actionButton(ns("act_analyze_hd_true"), label = "Analyze", class = "btn btn-primary btn-lg")
                    )
             ),
             column(9,
                    plotOutput(ns("hist_lr_rand"), height = "600px")
             )
           )
  )
}

#' hd_true_server
#'
#' @description The function to create the server module for the Hd true test
hd_true_server <- function(id, path_pack, rv_deconvo){
  moduleServer(
    id,
    function(input, output, session){
      rv_hd_true <- reactiveValues()
      rv_hd_true$lr_rand <- NULL

      observe({
        req(rv_deconvo)
        req(rv_hd_true)
        rv_hd_true <- hd_true_server_sub("hd_true_sub", path_pack, rv_deconvo, rv_hd_true)
      })

      return(rv_hd_true)
    }
  )
}

#' hd_true_server_sub
#'
#' @description The function to create the sub-server module for the Hd true test
hd_true_server_sub <- function(id, path_pack, rv_deconvo, rv_hd_true){
  moduleServer(
    id,
    function(input, output, session){
      data_deconvo <- rv_deconvo$data_deconvo

      if(is.null(data_deconvo)){
        disable("act_analyze_hd_true")
      }else{
        hnc_all <- data_deconvo$hnc_all

        enable("act_analyze_hd_true")

        output$fn_csp <- renderText({paste0(data_deconvo$fn_csp)})
        output$fn_af <- renderText({paste0(data_deconvo$fn_af)})

        output$noc <- renderUI({
          selectInput(session$ns("noc"), label = "Number of contributors", choices = hnc_all, selected = hnc_all[1], selectize = FALSE)
        })

        observeEvent(input$act_analyze_hd_true, {
          waiter_show(html = spin_3k(), color = "white")

          hnc <- as.numeric(input$noc)
          n_simu <- input$n_simu
          sub_source <- input$sub_source

          # Extract data from data_deconvo
          markers_use <- data_deconvo$markers_use
          peaks_list <- data_deconvo$peaks_list
          pop_al_list <- data_deconvo$pop_al_list
          pop_prob_list <- data_deconvo$pop_prob_list
          dt_anameth <- data_deconvo$dt_anameth
          theta <- dt_anameth$Value[dt_anameth$Subcategory == "Theta"]
          data_hnc_list <- data_deconvo$data_hnc_list

          data_h1 <- data_hnc_list[[which(hnc_all == hnc)]]
          error_message_h1 <- data_h1$error_message
          dt_pg_list_h1 <- data_h1$dt_pg_list

          dt_af <- data_deconvo$dt_af
          tmp <- set_al_prob_randgt(markers_use, dt_af)
          pop_al_list_randgt <- tmp$pop_al_list
          pop_prob_list_randgt <- tmp$pop_prob_list
          randgt <- t(generate_randgt(n_simu, pop_al_list_randgt, pop_prob_list_randgt))

          lr_rand <- rep(NA, n_simu)
          for(i in 1:n_simu){
            ref1 <- randgt[, c(2 * i - 1, 2 * i)]
            ref1 <- t(apply(ref1, 1, sort))
            refs <- paste0("Ref", i)
            dt_ref_rand <- data.table(SampleName = rep(refs, length(markers_use)),
                                      Marker = markers_use,
                                      Allele1 = ref1[, 1],
                                      Allele2 = ref1[, 2])

            log10like_hp <- assign_log10like(hnc, markers_use, peaks_list, pop_al_list, pop_prob_list, theta, dt_pg_list_h1, dt_ref_rand, refs, refs, sub_source)
            log10like_hd <- assign_log10like(hnc, markers_use, peaks_list, pop_al_list, pop_prob_list, theta, dt_pg_list_h1, dt_ref_rand, NULL, refs, sub_source)

            lr_rand[i] <- 10^(sum(log10like_hp) - sum(log10like_hd))
          }
          rv_hd_true$lr_rand <- lr_rand

          waiter_hide()
        })

        output$hist_lr_rand <- renderPlot({
          lr_rand <- rv_hd_true$lr_rand
          if(!is.null(lr_rand)){
            log10_lr_rand <- log10(lr_rand)
            log10_lr_rand[log10_lr_rand < -100] <- -100
            hist(log10_lr_rand)
          }
        })
      }

      return(rv_hd_true)
    }
  )
}
