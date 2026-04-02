#' deconvo_pg_ui
#'
#' @description The function to create the ui module for probabilistic genotyping
#' @param id The ui ID
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
                    fluidRow(
                      column(6, uiOutput(ns("slider_x"))),
                      column(6, uiOutput(ns("slider_y")))
                    ),
                    uiOutput(ns("ph_graph_ui"))
             )
           )
  )
}

#' deconvo_pg_server
#'
#' @description The function to create the server module for probabilistic genotyping
#' @param id The server ID
#' @param rv_deconvo The reactive values derived from the deconvolution server
deconvo_pg_server <- function(id, rv_deconvo){
  moduleServer(
    id,
    function(input, output, session){
      rv_deconvo_pg <- reactiveValues()
      rv_deconvo_pg$markers_use <- NULL
      rv_deconvo_pg$hnc_all <- NULL
      rv_deconvo_pg$data_hnc_list <- NULL
      rv_deconvo_pg$dt_anameth <- NULL
      rv_deconvo_pg$list_dt_average_sr <- NULL
      rv_deconvo_pg$dt_kit <- NULL
      rv_deconvo_pg$peaks_list <- NULL
      rv_deconvo_pg$heights_list <- NULL
      rv_deconvo_pg$ae_use <- NULL
      rv_deconvo_pg$n_rep <- NULL
      rv_deconvo_pg$cand_contrib <- NULL
      rv_deconvo_pg$dt_pg_per_comb <- NULL
      rv_deconvo_pg$dt_pg_per_contrib <- NULL
      rv_deconvo_pg$error_message <- NULL
      rv_deconvo_pg$df_combined_list <- NULL
      rv_deconvo_pg$peaks_plot <- NULL
      rv_deconvo_pg$peaks_label <- NULL
      rv_deconvo_pg$x_min <- 0
      rv_deconvo_pg$x_max <- 50
      rv_deconvo_pg$y_max <- 50000
      rv_deconvo_pg$plot_mk <- NULL
      rv_deconvo_pg$at <- 0

      observe({
        req(rv_deconvo$data_deconvo)

        data_deconvo <- rv_deconvo$data_deconvo
        if(!is.null(data_deconvo)){
          rv_deconvo_pg$markers_use <- data_deconvo$markers_use
          rv_deconvo_pg$hnc_all <- data_deconvo$hnc_all
          rv_deconvo_pg$data_hnc_list <- data_deconvo$data_hnc_list
          rv_deconvo_pg$dt_anameth <- data_deconvo$dt_anameth
          rv_deconvo_pg$list_dt_average_sr <- data_deconvo$list_dt_average_sr
          rv_deconvo_pg$dt_kit <- data_deconvo$dt_kit
          rv_deconvo_pg$peaks_list <- data_deconvo$peaks_list
          rv_deconvo_pg$heights_list <- data_deconvo$heights_list
          rv_deconvo_pg$ae_use <- data_deconvo$ae_use
          rv_deconvo_pg$n_rep <- length(rv_deconvo_pg$heights_list[[1]])
          rv_deconvo_pg$cand_contrib <- paste0("U", 1:(rv_deconvo_pg$hnc_all[1]))
        }
      })

      output$select_meth_pg <- renderUI({
        req(rv_deconvo_pg$hnc_all)
        selectInput(session$ns("select_meth_pg"), label = "Method", choices = c("Per combination", "Per contributor"), selected = "Per combination", selectize = FALSE)
      })

      output$select_noc_pg <- renderUI({
        req(rv_deconvo_pg$hnc_all)
        hnc_all <- rv_deconvo_pg$hnc_all
        selectInput(session$ns("select_noc_pg"), label = "Number of contributors", choices = hnc_all, selected = hnc_all[1], selectize = FALSE)
      })

      output$select_mk_pg <- renderUI({
        req(rv_deconvo_pg$markers_use)
        markers_use <- rv_deconvo_pg$markers_use
        selectInput(session$ns("select_mk_pg"), label = "Marker", choices = markers_use, selected = markers_use[1], selectize = FALSE)
      })

      output$select_contrib_pg <- renderUI({
        req(rv_deconvo_pg$cand_contrib)
        cand_contrib <- rv_deconvo_pg$cand_contrib
        selectInput(session$ns("select_contrib_pg"), label = "Contributor", choices = cand_contrib, selected = cand_contrib[1], selectize = FALSE)
      })

      react_input_pg <- reactive({
        list(input$select_meth_pg, input$select_noc_pg, input$select_mk_pg, input$select_contrib_pg)
      })

      observeEvent(react_input_pg(), ignoreInit = TRUE, {
        req(rv_deconvo_pg$hnc_all)
        req(rv_deconvo_pg$data_hnc_list)
        hnc_all <- rv_deconvo_pg$hnc_all
        data_hnc_list <- rv_deconvo_pg$data_hnc_list

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

      output$dt_pg_per_comb <- renderDataTable(server = TRUE, {
        req(rv_deconvo_pg$dt_pg_per_comb)
        dt_pg_per_comb <- rv_deconvo_pg$dt_pg_per_comb
        label_dt_pg <- names(dt_pg_per_comb)
        hnc <- as.numeric(input$select_noc_pg)

        datatable(
          dt_pg_per_comb,
          selection = list(mode = "single", target = "row"),
          options = list(iDisplayLength = 10, ordering = FALSE, autowidth = TRUE, searching = FALSE,
                         columnDefs = list(list(targets = c(grep("Allele", label_dt_pg), which(is.element(label_dt_pg, "Logprob_gtc")), grep("Loglike", label_dt_pg)) - 1, visible = FALSE))),
          rownames = FALSE
        ) %>%
          formatSignif(columns = c("Probability"), digits = 3) %>%
          formatStyle(names(dt_pg_per_comb), lineHeight = "80%")
      })

      output$dt_pg_per_contrib <- renderDataTable(server = TRUE, {
        req(rv_deconvo_pg$dt_pg_per_contrib)
        dt_pg_per_contrib <- rv_deconvo_pg$dt_pg_per_contrib

        datatable(
          dt_pg_per_contrib,
          selection = list(mode = "single", target = "row"),
          options = list(iDisplayLength = 10, ordering = FALSE, autowidth = TRUE, searching = FALSE),
          rownames = FALSE
        ) %>%
          formatSignif(columns = c("Probability"), digits = 3) %>%
          formatStyle(names(dt_pg_per_contrib), lineHeight = "80%")
      })

      observeEvent(input$dt_pg_per_comb_rows_selected, {
        req(rv_deconvo_pg$markers_use)
        req(rv_deconvo_pg$hnc_all)
        req(rv_deconvo_pg$data_hnc_list)
        req(rv_deconvo_pg$list_dt_average_sr)
        req(rv_deconvo_pg$dt_kit)
        req(rv_deconvo_pg$peaks_list)
        req(rv_deconvo_pg$heights_list)
        req(rv_deconvo_pg$ae_use)
        req(rv_deconvo_pg$dt_anameth)
        req(rv_deconvo_pg$n_rep)
        req(rv_deconvo_pg$dt_pg_per_comb)

        markers_use <- rv_deconvo_pg$markers_use
        hnc_all <- rv_deconvo_pg$hnc_all
        data_hnc_list <- rv_deconvo_pg$data_hnc_list
        list_dt_average_sr <- rv_deconvo_pg$list_dt_average_sr
        dt_kit <- rv_deconvo_pg$dt_kit
        peaks_list <- rv_deconvo_pg$peaks_list
        heights_list <- rv_deconvo_pg$heights_list
        ae_use <- rv_deconvo_pg$ae_use
        dt_anameth <- rv_deconvo_pg$dt_anameth
        n_rep <- rv_deconvo_pg$n_rep
        dt_pg_per_comb <- rv_deconvo_pg$dt_pg_per_comb

        at_all_mk <- dt_anameth$Value[dt_anameth$Category == "AT"]
        markers_at <- dt_anameth$Subcategory[dt_anameth$Category == "AT"]
        base_size <- dt_anameth$Value[dt_anameth$Subcategory == "Base size"]

        index_select <- input$dt_pg_per_comb_rows_selected

        if(!is.null(index_select)){
          hnc <- input$select_noc_pg
          mk <- input$select_mk_pg
          index_hnc <- which(hnc_all == hnc)
          data_hnc <- data_hnc_list[[index_hnc]]
          heights_all_rep <- heights_list[[mk]]
          at <- at_all_mk[markers_at == mk]

          gt_select_list <- vector(mode = "list", length = hnc)
          for(i in 1:hnc){
            tmp <- strsplit(dt_pg_per_comb[[paste0("U", i)]][index_select], ", ")[[1]]
            tmp[tmp == "Q"] <- "99"
            gt_select_list[[i]] <- as.numeric(tmp)
          }

          dt_param <- data_hnc$dt_param
          param_names <- dt_param[["Parameter"]]
          est_vals <- dt_param[["Estimated value"]]
          temp_total <- est_vals[param_names == "Template amount"]
          mp <- est_vals[param_names == "Mixture proportion"]
          temps <- temp_total * mp
          deg <- est_vals[param_names == "Degradation"]
          rep_ae <- est_vals[param_names == "Replicate amplification efficiency"]
          ae <- ae_use[markers_use == mk]

          dt_average_sr <- list_dt_average_sr[[mk]]

          dt_size <- dt_kit[Marker == mk]
          rep_unit_len <- dt_size[["Repeat Unit Length"]][1]
          dt_size <- dt_size[, c("Allele", "Size")]

          df_combined_list <- vector(mode = "list", length = n_rep)
          for(i in 1:n_rep){
            # Calculate expected peak heights under the selected genotype combination
            eph <- calc_eph(gt_select_list, temps, deg, rep_ae[i], ae, dt_average_sr, dt_size, base_size)
            peaks_eph <- as.numeric(names(eph))

            # Remove 99 from oph data
            peaks_oph <- as.numeric(peaks_list[[mk]])
            ophs_one_rep <- heights_all_rep[[i]]
            ophs_one_rep <- ophs_one_rep[peaks_oph != 99]
            peaks_oph <- peaks_oph[peaks_oph != 99]

            if(i == 1){
              peaks_all <- sort(unique(c(peaks_oph, peaks_eph)))
              peaks_plot <- create_peaks_plot(peaks_all, rep_unit_len)
              peaks_label <- peaks_all
              peaks_label[peaks_label == 991] <- "Q (U1)"
              peaks_label[peaks_label == 992] <- "Q (U2)"
              peaks_label[peaks_label == 993] <- "Q (U3)"
              peaks_label[peaks_label == 994] <- "Q (U4)"
            }

            eph_graph <- oph_graph <- rep(0, length(peaks_all))
            eph_graph[match(peaks_eph, peaks_all)] <- eph
            oph_graph[match(peaks_oph, peaks_all)] <- ophs_one_rep

            df_eph <- data.frame(peak = peaks_plot,
                                 height = eph_graph,
                                 Label = "Expected")

            df_oph <- data.frame(peak = peaks_plot,
                                 height = oph_graph,
                                 Label = "Observed")

            df_combined_list[[i]] <- rbind(df_eph, df_oph)
          }

          rv_deconvo_pg$peaks_plot <- peaks_plot
          rv_deconvo_pg$peaks_label <- peaks_label
          rv_deconvo_pg$df_combined_list <- df_combined_list
          rv_deconvo_pg$x_min <- floor(min(peaks_plot) - 1)
          rv_deconvo_pg$x_max <- ceiling(max(peaks_plot) + 1)
          rv_deconvo_pg$y_max <- max(c(eph_graph, oph_graph))
          rv_deconvo_pg$plot_mk <- mk
          rv_deconvo_pg$at <- at
        }
      })

      observeEvent(input$dt_pg_per_contrib_rows_selected, {
        req(rv_deconvo_pg$markers_use)
        req(rv_deconvo_pg$hnc_all)
        req(rv_deconvo_pg$data_hnc_list)
        req(rv_deconvo_pg$list_dt_average_sr)
        req(rv_deconvo_pg$dt_kit)
        req(rv_deconvo_pg$peaks_list)
        req(rv_deconvo_pg$heights_list)
        req(rv_deconvo_pg$ae_use)
        req(rv_deconvo_pg$dt_anameth)
        req(rv_deconvo_pg$n_rep)
        req(rv_deconvo_pg$dt_pg_per_contrib)

        markers_use <- rv_deconvo_pg$markers_use
        hnc_all <- rv_deconvo_pg$hnc_all
        data_hnc_list <- rv_deconvo_pg$data_hnc_list
        list_dt_average_sr <- rv_deconvo_pg$list_dt_average_sr
        dt_kit <- rv_deconvo_pg$dt_kit
        peaks_list <- rv_deconvo_pg$peaks_list
        heights_list <- rv_deconvo_pg$heights_list
        ae_use <- rv_deconvo_pg$ae_use
        dt_anameth <- rv_deconvo_pg$dt_anameth
        n_rep <- rv_deconvo_pg$n_rep
        dt_pg_per_contrib <- rv_deconvo_pg$dt_pg_per_contrib

        at_all_mk <- dt_anameth$Value[dt_anameth$Category == "AT"]
        markers_at <- dt_anameth$Subcategory[dt_anameth$Category == "AT"]
        base_size <- dt_anameth$Value[dt_anameth$Subcategory == "Base size"]

        index_select <- input$dt_pg_per_contrib_rows_selected

        if(!is.null(index_select)){
          hnc <- input$select_noc_pg
          mk <- input$select_mk_pg
          contrib <- input$select_contrib_pg
          index_hnc <- which(hnc_all == hnc)
          data_hnc <- data_hnc_list[[index_hnc]]
          heights_all_rep <- heights_list[[mk]]
          at <- at_all_mk[markers_at == mk]

          gt_select_list <- vector(mode = "list", length = 1)
          tmp <- strsplit(dt_pg_per_contrib[[contrib]][index_select], ", ")[[1]]
          tmp[tmp == "Q"] <- "99"
          gt_select_list[[1]] <- as.numeric(tmp)

          index_contrib <- as.integer(gsub("U", "", contrib))

          dt_param <- data_hnc$dt_param
          param_names <- dt_param[["Parameter"]]
          est_vals <- dt_param[["Estimated value"]]
          temp_total <- est_vals[param_names == "Template amount"]
          mp <- est_vals[param_names == "Mixture proportion"]
          temps <- temp_total * mp
          deg <- est_vals[param_names == "Degradation"]
          rep_ae <- est_vals[param_names == "Replicate amplification efficiency"]
          ae <- ae_use[markers_use == mk]

          temp_c1 <- temps[index_contrib]
          deg_c1 <- deg[index_contrib]

          dt_average_sr <- list_dt_average_sr[[mk]]

          dt_size <- dt_kit[Marker == mk]
          rep_unit_len <- dt_size[["Repeat Unit Length"]][1]
          dt_size <- dt_size[, c("Allele", "Size")]

          df_combined_list <- vector(mode = "list", length = n_rep)
          for(i in 1:n_rep){
            # Calculate expected peak heights under the selected genotype combination
            eph <- calc_eph(gt_select_list, temp_c1, deg_c1, rep_ae[i], ae, dt_average_sr, dt_size, base_size)
            peaks_eph <- as.numeric(names(eph))

            # Remove 99 from oph data
            peaks_oph <- as.numeric(peaks_list[[mk]])
            ophs_one_rep <- heights_all_rep[[i]]
            ophs_one_rep <- ophs_one_rep[peaks_oph != 99]
            peaks_oph <- peaks_oph[peaks_oph != 99]

            if(i == 1){
              peaks_all <- sort(unique(c(peaks_oph, peaks_eph)))
              peaks_plot <- create_peaks_plot(peaks_all, rep_unit_len)
              peaks_label <- peaks_all
              peaks_label[peaks_label == 991] <- paste0("Q (", contrib, ")")
            }

            eph_graph <- oph_graph <- rep(0, length(peaks_all))
            eph_graph[match(peaks_eph, peaks_all)] <- eph
            oph_graph[match(peaks_oph, peaks_all)] <- ophs_one_rep

            df_eph <- data.frame(peak = peaks_plot,
                                 height = eph_graph,
                                 Label = "Expected")

            df_oph <- data.frame(peak = peaks_plot,
                                 height = oph_graph,
                                 Label = "Observed")

            df_combined_list[[i]] <- rbind(df_eph, df_oph)
          }

          rv_deconvo_pg$peaks_plot <- peaks_plot
          rv_deconvo_pg$peaks_label <- peaks_label
          rv_deconvo_pg$df_combined_list <- df_combined_list
          rv_deconvo_pg$x_min <- floor(min(peaks_plot) - 1)
          rv_deconvo_pg$x_max <- ceiling(max(peaks_plot) + 1)
          rv_deconvo_pg$y_max <- max(c(eph_graph, oph_graph))
          rv_deconvo_pg$plot_mk <- mk
          rv_deconvo_pg$at <- at
        }
      })

      output$slider_x <- renderUI({
        x_min <- rv_deconvo_pg$x_min
        x_max <- rv_deconvo_pg$x_max

        sliderInput(session$ns("x_range"),
                    "Range of the x-axis",
                    min = x_min,
                    max = x_max,
                    value = c(x_min, x_max),
                    step = 1)
      })

      output$slider_y <- renderUI({
        y_max <- rv_deconvo_pg$y_max

        sliderInput(session$ns("y_range"),
                    "Range of the y-axis (log10 scale)",
                    min = 0,
                    max = ceiling(log10(y_max)),
                    value = round(log10(y_max), 2),
                    step = 0.01)
      })

      output$ph_graph_ui <- renderUI({
        req(rv_deconvo_pg$n_rep)
        n_rep <- rv_deconvo_pg$n_rep

        select_meth_pg <- input$select_meth_pg
        if(!is.null(select_meth_pg)){
          plot_height <- paste0(350 * n_rep, "px")
          plotOutput(session$ns("ph_graph"), width = "100%", height = plot_height)
        }
      })

      output$ph_graph <- renderPlot({
        req(rv_deconvo_pg$n_rep)
        req(rv_deconvo_pg$peaks_plot)
        req(rv_deconvo_pg$peaks_label)
        req(rv_deconvo_pg$df_combined_list)
        req(rv_deconvo_pg$plot_mk)
        req(rv_deconvo_pg$at)

        n_rep <- rv_deconvo_pg$n_rep
        peaks_plot <- rv_deconvo_pg$peaks_plot
        peaks_label <- rv_deconvo_pg$peaks_label
        df_combined_list <- rv_deconvo_pg$df_combined_list
        mk <- rv_deconvo_pg$plot_mk
        at <- rv_deconvo_pg$at

        ggplot_list <- vector(mode = "list", length = n_rep)
        options(warn = -1)
        for(i in 1:n_rep){
          df_combined <- df_combined_list[[i]]

          ggplot_list[[i]] <- ggplot(df_combined, aes(x = peak, y = height, fill = Label)) +
            geom_col(position = "dodge", width = 0.1) +
            geom_hline(yintercept = at, linetype = "dashed", linewidth = 0.5) +
            scale_x_continuous(
              breaks = peaks_plot,
              labels = peaks_label,
              limits = c(input$x_range[1], input$x_range[2])
            ) +
            coord_cartesian(ylim = c(0, 10^input$y_range)) +
            labs(title = paste0(mk, ", Replicate ", i), x = "Allele", y = "Peak height (RFU)", fill = NULL) +
            theme_minimal() +
            theme(
              plot.title = element_text(size = 16),
              legend.text = element_text(size = 14),
              axis.title.x = element_text(size = 14, margin = margin(t = 20)),
              axis.title.y = element_text(size = 14, margin = margin(r = 20)),
              axis.text.x = element_text(size = 14),
              axis.text.y = element_text(size = 14)
            )
        }
        grid.arrange(grobs = ggplot_list[1:n_rep], ncol = 1, nrow = n_rep)
        options(warn = 0)
      })
    }
  )
}
