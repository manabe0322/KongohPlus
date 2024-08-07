#' hd_true_ui
#'
#' @description The function to create the ui module for the Hd true test
hd_true_ui <- function(id){
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
                      numericInput(ns("n_simu"), label = h4("Number of simulations"), value = 100),
                      br(),
                      actionButton(ns("act_analyze_hd_true"), label = "Analyze", class = "btn btn-primary btn-lg"),
                      span(textOutput(ns("error_message")), style = "color:red; font-weight:bold;")
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

      rv_hd_true_setting <- reactiveValues()
      rv_hd_true_setting$data_deconvo <- NULL

      observe({
        req(rv_deconvo)
        rv_hd_true_setting$data_deconvo <- rv_deconvo$data_deconvo
      })

      if(is.null(rv_hd_true_setting$data_deconvo)){
        disable("act_analyze_hd_true")
        output$error_message <- renderText({paste0("Perform deconvolution before starting Hd true test!")})
      }else{
        data_deconvo <- rv_hd_true_setting$data_deconvo
        enable("act_analyze_hd_true")
        output$error_message <- renderText({paste0("")})

        output$fn_csp <- renderText({paste0(data_deconvo$fn_csp)})
        output$fn_af <- renderText({paste0(data_deconvo$fn_af)})

        observeEvent(input$act_analyze_hd_true, {
          waiter_show(html = spin_3k(), color = "white")

          n_simu <- input$n_simu

          # Extract data from data_deconvo
          hnc <- data_deconvo$hnc
          markers_use <- data_deconvo$markers_use
          peaks_list <- data_deconvo$peaks_list
          pop_al_list <- data_deconvo$pop_al_list
          pop_prob_list <- data_deconvo$pop_prob_list
          theta <- data_deconvo$theta
          gtc_use <- data_deconvo$gtc_use
          lpd_gtc_use_1st <- data_deconvo$lpd_gtc_use_1st
          lpd_gtc_use_2nd <- data_deconvo$lpd_gtc_use_2nd

          dt_af <- data_deconvo$dt_af
          dt_anameth <- fread(paste0(path_pack, "/extdata/anameth/", data_deconvo$anameth, ".csv"))
          method_al_prob <- dt_anameth$Value[dt_anameth$Category == "Allele_probability" & dt_anameth$Subcategory == "Method"]
          tmp <- set_al_prob_randgt(markers_use, dt_af, method_al_prob)
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
            data_hp <- assign_log10like_randgt(hnc, markers_use, peaks_list, pop_al_list, pop_prob_list, theta, gtc_use, lpd_gtc_use_1st, lpd_gtc_use_2nd, dt_ref_rand, refs, refs)
            data_hd <- assign_log10like_randgt(hnc, markers_use, peaks_list, pop_al_list, pop_prob_list, theta, gtc_use, lpd_gtc_use_1st, lpd_gtc_use_2nd, dt_ref_rand, NULL, refs)
            lr_rand[i] <- 10^(data_hp$combined_log10like_adopt - data_hd$combined_log10like_adopt)
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
