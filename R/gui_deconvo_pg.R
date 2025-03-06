#' deconvo_pg_ui
#'
#' @description The function to create the ui module for probabilistic genotyping
deconvo_pg_ui <- function(id){
  ns <- NS(id)

  tabPanel("Probabilistic genotyping",
           br(),
           fluidRow(
             column(6,
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
                      uiOutput(ns("pg"))
                    )
             ),
             column(6,
                    plotOutput(ns("ph_graph"))
             )
           )
  )
}

#' aggregate_pg_contrib
#'
#' @description The function to change the format of the pg table from per combination to per contributor
#' @param dt_pg The data.table of the probabilistic genotyping in one marker
#' @param contrib The name of the selected contributor
aggregate_pg_contrib <- function(dt_pg, contrib){
  aggregated_data <- dt_pg %>%
    group_by_at(which(names(dt_pg) == contrib)) %>%
    summarise(Probability = sum(Probability))
  aggregated_data <- aggregated_data %>%
    arrange(desc(Probability))
  return(aggregated_data)
}

#' deconvo_pg_server
#'
#' @description The function to create the server module for probabilistic genotyping
deconvo_pg_server <- function(id, rv_deconvo){
  moduleServer(
    id,
    function(input, output, session){
      data_deconvo <- rv_deconvo$data_deconvo
      if(!is.null(data_deconvo)){
        markers_use <- data_deconvo$markers_use
        hnc_all <- data_deconvo$hnc_all
        data_hnc_list <- data_deconvo$data_hnc_list
        dt_anameth <- data_deconvo$dt_anameth
        list_dt_average_sr <- data_deconvo$list_dt_average_sr
        dt_kit <- data_deconvo$dt_kit
        peaks_list <- data_deconvo$peaks_list
        heights_list <- data_deconvo$heights_list
        mean_size <- data_deconvo$mean_size

        rv_deconvo_pg <- reactiveValues()
        rv_deconvo_pg$cand_contrib <- paste0("U", 1:hnc_all[1])
        rv_deconvo_pg$dt_pg_per_comb <- NULL
        rv_deconvo_pg$dt_pg_per_contrib <- NULL
        rv_deconvo_pg$error_message <- NULL

        output$select_meth_pg <- renderUI({
          selectInput(session$ns("select_meth_pg"), label = "Method", choices = c("Per combination", "Per contributor"), selected = "Per combination", selectize = FALSE)
        })

        output$select_noc_pg <- renderUI({
          selectInput(session$ns("select_noc_pg"), label = "Number of contributors", choices = hnc_all, selected = hnc_all[1], selectize = FALSE)
        })

        output$select_mk_pg <- renderUI({
          selectInput(session$ns("select_mk_pg"), label = "Marker", choices = markers_use, selected = markers_use[1], selectize = FALSE)
        })

        output$select_contrib_pg <- renderUI({
          selectInput(session$ns("select_contrib_pg"), label = "Contributor", choices = rv_deconvo_pg$cand_contrib, selected = rv_deconvo_pg$cand_contrib[1], selectize = FALSE)
        })

        react_input_pg <- reactive({
          list(input$select_meth_pg, input$select_noc_pg, input$select_mk_pg, input$select_contrib_pg)
        })

        observeEvent(react_input_pg(), ignoreInit = TRUE, {
          hnc <- as.numeric(input$select_noc_pg)
          rv_deconvo_pg$cand_contrib <- paste0("U", 1:hnc)

          index_hnc <- which(hnc_all == hnc)
          data_hnc <- data_hnc_list[[index_hnc]]
          rv_deconvo_pg$error_message <- data_hnc$error_message
          if(is.null(rv_deconvo_pg$error_message)){
            dt_pg_list <- data_hnc$dt_pg_list

            if(input$select_meth_pg == "Per combination"){
              rv_deconvo_pg$dt_pg_per_comb <- dt_pg_list[[input$select_mk_pg]]
              rv_deconvo_pg$dt_pg_per_contrib <- NULL
            }else{
              rv_deconvo_pg$dt_pg_per_comb <- NULL
              rv_deconvo_pg$dt_pg_per_contrib <- aggregate_pg_contrib(dt_pg_list[[input$select_mk_pg]], input$select_contrib_pg)
            }
          }else{
            rv_deconvo_pg$dt_pg_per_comb <- NULL
            rv_deconvo_pg$dt_pg_per_contrib <- NULL
          }
        })

        output$pg <- renderUI({
          dt_pg_per_comb <- rv_deconvo_pg$dt_pg_per_comb
          dt_pg_per_contrib <- rv_deconvo_pg$dt_pg_per_contrib
          error_message <- rv_deconvo_pg$error_message
          if(is.null(error_message)){
            if(is.null(dt_pg_per_contrib)){
              dataTableOutput(session$ns("dt_pg_per_comb"))
            }else{
              dataTableOutput(session$ns("dt_pg_per_contrib"))
            }
          }else{
            HTML(error_message)
          }
        })

        output$dt_pg_per_comb <- renderDataTable(server = FALSE, {
          dt_pg_per_comb <- rv_deconvo_pg$dt_pg_per_comb
          hnc <- as.numeric(input$select_noc_pg)
          if(!is.null(dt_pg_per_comb)){
            datatable(
              dt_pg_per_comb,
              selection = list(mode = "single", target = "row"),
              options = list(iDisplayLength = 10, ordering = FALSE, autowidth = TRUE, searching = FALSE,
                             columnDefs = list(list(targets = c(hnc, hnc + 1, hnc + 2), visible = FALSE))),
              rownames = FALSE
            ) %>%
              formatSignif(columns = c("Probability"), digits = 3) %>%
              formatStyle(names(dt_pg_per_comb), lineHeight = "80%")
          }
        })

        output$dt_pg_per_contrib <- renderDataTable(server = FALSE, {
          dt_pg_per_contrib <- rv_deconvo_pg$dt_pg_per_contrib
          if(!is.null(dt_pg_per_contrib)){
            datatable(
              dt_pg_per_contrib,
              selection = list(mode = "single", target = "row"),
              options = list(iDisplayLength = 10, ordering = FALSE, autowidth = TRUE, searching = FALSE),
              rownames = FALSE
            ) %>%
              formatSignif(columns = c("Probability"), digits = 3) %>%
              formatStyle(names(dt_pg_per_contrib), lineHeight = "80%")
          }
        })

        observeEvent(input$dt_pg_per_comb_rows_selected, {
          index_select <- input$dt_pg_per_comb_rows_selected

          if(is.null(index_select)){

          }else{
            hnc <- input$select_noc_pg
            mk <- input$select_mk_pg
            index_hnc <- which(hnc_all == hnc)
            data_hnc <- data_hnc_list[[index_hnc]]
            heights_all_rep <- heights_list[[mk]]

            dt_pg_list <- data_hnc$dt_pg_list
            dt_pg <- dt_pg_list[[mk]]
            gt_select_list <- vector(mode = "list", length = hnc)
            for(i in 1:hnc){
              tmp <- strsplit(dt_pg[[paste0("U", i)]][index_select], ", ")[[1]]
              tmp[tmp == "Q"] <- "99"
              gt_select_list[[i]] <- as.numeric(tmp)
            }

            dt_diag <- data_hnc$dt_diag
            param_names <- dt_diag[["Parameter"]]
            est_vals <- dt_diag[["Estimated value"]]
            tphs <- est_vals[param_names == "Total peak height"]
            mp <- est_vals[param_names == "Mixture proportion"]
            deg <- est_vals[param_names == "Degradation"]
            ae_use <- est_vals[param_names == "Amplification efficiency"]
            markers_use <- dt_diag$Marker[param_names == "Amplification efficiency"]
            ae <- ae_use[markers_use == mk]

            dt_average_sr <- list_dt_average_sr[[mk]]

            dt_size <- dt_kit[Marker == mk]
            dt_size <- dt_size[, c("Allele", "Size")]

            # Calculate expected peak heights under the selected genotype combination
            n_rep <- length(tphs)
            eph_list <- vector(mode = "list", length = n_rep)
            for(i in 1:n_rep){
              eph_list[[i]] <- calc_eph(gt_select_list, tphs[i], mp, deg, ae, dt_average_sr, dt_size, mean_size)
            }

            # Plot (1st replicate)
            peaks_eph <- as.numeric(names(eph_list[[1]]))
            peaks_oph <- as.numeric(peaks_list[[mk]])
            peaks_all <- sort(unique(c(peaks_oph, peaks_eph)))
            peaks_all[peaks_all == 99] <- max(setdiff(peaks_all, 99)) + 3
            x_min <- min(peaks_all)
            x_max <- max(peaks_all)
            peaks_eph[peaks_eph == 99] <- x_max
            peaks_oph[peaks_oph == 99] <- x_max

            ephs <- ophs <- rep(0, length(peaks_all))
            ephs[match(peaks_eph, peaks_all)] <- eph_list[[1]]
            ophs[match(peaks_oph, peaks_all)] <- heights_all_rep[[1]]

            df_eph <- data.frame(peak = peaks_all,
                                 height = ephs,
                                 Label = "Expected")

            df_oph <- data.frame(peak = peaks_all,
                                 height = ophs,
                                 Label = "Observed")

            df_combined <- rbind(df_eph, df_oph)

            output$ph_graph <- renderPlot({
              ggplot(df_combined, aes(x = peak, y = height, fill = Label)) +
                geom_col(position = "dodge", width = 0.2) +
                labs(title = "", x = "Allele", y = "Peak height (RFU)") +
                theme_minimal()
            })
          }
        })

        observeEvent(input$dt_pg_per_contrib_rows_selected, {
          index_select <- input$dt_pg_per_contrib_rows_selected

          if(is.null(index_select)){

          }else{

          }
        })
      }
    }
  )
}
