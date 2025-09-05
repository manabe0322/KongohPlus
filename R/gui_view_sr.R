#' view_sr_ui
#'
#' @description The function to create the ui module for viewing stutter ratios
view_sr_ui <- function(id){
  ns <- NS(id)

  tabPanel(
    "View stutter ratios",
    titlePanel("View stutter ratios"),
    br(),
    fluidRow(
      column(3,
             wellPanel(
               h4("Stutter data"),
               fluidRow(
                 column(3, br(), shinyFilesButton(ns("load_stutter_data"), label = "Load", title = "Select a file", multiple = FALSE)),
                 column(9, uiOutput(ns("fn_stutter_data")))
               ),
               uiOutput(ns("select_mk")),
               uiOutput(ns("select_st_type"))
             )
      ),
      column(4,
             wellPanel(
               h4("Average stutter ratio"),
               dataTableOutput(ns("dt_average_sr_display"))
             )
      ),
      column(5,
             wellPanel(
               h4("Experimental data"),
               uiOutput(ns("ui_plot_sr")),
               br(),
               dataTableOutput(ns("dt_sr_data_display"))
             )
      )
    )
  )
}

#' view_sr_server
#'
#' @description The function to create the server module for viewing stutter ratios
view_sr_server <- function(id, path_pack){
  moduleServer(
    id,
    function(input, output, session){
      rv_view_sr <- reactiveValues()
      rv_view_sr$fn_stutter_data <- NULL
      rv_view_sr$markers_auto <- character(0)
      rv_view_sr$list_dt_average_sr <- NULL
      rv_view_sr$dt_average_sr <- NULL
      rv_view_sr$dt_setting_est_sr <- NULL
      rv_view_sr$dt_sr_data <- NULL
      rv_view_sr$name_kit_stutter <- NULL
      rv_view_sr$dt_st_type <- NULL
      rv_view_sr$meth_all_mk <- NULL
      rv_view_sr$meth <- NULL
      rv_view_sr$dt_st_mk <- NULL
      rv_view_sr$dt_kit_mk <- NULL
      rv_view_sr$dt_sr_data_display <- NULL
      rv_view_sr$df_select_sr <- NULL

      shinyFileChoose(input, "load_stutter_data", roots = c(documents = paste0(path_pack, "/extdata/stutter")), session = session, filetypes = c("RData"))

      observeEvent(input$load_stutter_data, {
        req(input$load_stutter_data)
        path_stutter_data <- parseFilePaths(roots = c(documents = paste0(path_pack, "/extdata/stutter")), input$load_stutter_data)$datapath

        if(length(path_stutter_data) > 0){
          rv_view_sr$fn_stutter_data <- basename(path_stutter_data)

          load(path_stutter_data)

          rv_view_sr$markers_auto <- names(data_est_sr$list_dt_average_sr)
          rv_view_sr$list_dt_average_sr <- data_est_sr$list_dt_average_sr
          rv_view_sr$dt_setting_est_sr <- data_est_sr$dt_setting_est_sr
          rv_view_sr$dt_sr_data <- data_est_sr$dt_sr_data
          rv_view_sr$name_kit_stutter <- data_est_sr$name_kit_stutter
        }
      })

      output$fn_stutter_data <- renderUI({
        disabled(textInput(session$ns("fn_stutter_data"), label = "", value = rv_view_sr$fn_stutter_data))
      })

      output$select_mk <- renderUI({
        markers_auto <- rv_view_sr$markers_auto

        if(length(markers_auto) == 0){
          selectInput(session$ns("select_mk"), label = h4("Marker"), choices = markers_auto, selected = NULL, multiple = FALSE, selectize = FALSE)
        }else{
          selectInput(session$ns("select_mk"), label = h4("Marker"), choices = markers_auto, selected = markers_auto[1], multiple = FALSE, selectize = FALSE)
        }
      })

      output$select_st_type <- renderUI({
        fn_stutter_data <- rv_view_sr$fn_stutter_data

        if(is.null(fn_stutter_data)){
          selectInput(session$ns("select_st_type"), label = h4("Stutter type"), choices = character(0), selected = NULL, multiple = FALSE, selectize = FALSE)
        }else{
          selectInput(session$ns("select_st_type"), label = h4("Stutter type"), choices = c("Back stutter", "Forward stutter", "Double-back stutter", "Minus 2-nt stutter"), selected = "Back stutter", multiple = FALSE, selectize = FALSE)
        }
      })

      observeEvent(list(input$select_mk, input$select_st_type), {
        select_mk <- input$select_mk
        select_st_type <- input$select_st_type

        rv_view_sr$dt_average_sr <- NULL
        rv_view_sr$dt_st_mk <- NULL
        rv_view_sr$dt_st_type <- NULL
        rv_view_sr$meth_all_mk <- NULL
        rv_view_sr$meth <- NULL
        rv_view_sr$dt_st_mk <- NULL
        rv_view_sr$dt_kit_mk <- NULL
        rv_view_sr$df_select_sr <- NULL

        if(!is.null(select_mk) && !is.null(select_st_type)){
          if(!is.null(rv_view_sr$list_dt_average_sr)){
            rv_view_sr$dt_average_sr <- rv_view_sr$list_dt_average_sr[[select_mk]]
          }

          dt_setting_est_sr <- rv_view_sr$dt_setting_est_sr
          dt_sr_data <- rv_view_sr$dt_sr_data
          name_kit_stutter <- rv_view_sr$name_kit_stutter

          if(!is.null(dt_setting_est_sr) && !is.null(dt_sr_data) && !is.null(name_kit_stutter)){

            if(select_st_type == "Back stutter"){
              dt_st_type <- dt_sr_data[dt_sr_data$Bool_BS]
              meth_all_mk <- dt_setting_est_sr[["BSR model"]]
              meth <- meth_all_mk[dt_setting_est_sr$Marker == select_mk]
            }else if(select_st_type == "Forward stutter"){
              dt_st_type <- dt_sr_data[dt_sr_data$Bool_FS]
              meth_all_mk <- dt_setting_est_sr[["FSR model"]]
              meth <- meth_all_mk[dt_setting_est_sr$Marker == select_mk]
            }else if(select_st_type == "Double-back stutter"){
              dt_st_type <- dt_sr_data[dt_sr_data$Bool_DS]
              meth_all_mk <- dt_setting_est_sr[["DSR model"]]
              meth <- meth_all_mk[dt_setting_est_sr$Marker == select_mk]
            }else if(select_st_type == "Minus 2-nt stutter"){
              dt_st_type <- dt_sr_data[dt_sr_data$Bool_M2S]
              meth_all_mk <- dt_setting_est_sr[["M2SR model"]]
              meth <- meth_all_mk[dt_setting_est_sr$Marker == select_mk]
            }else{
              dt_st_type <- dt_sr_data[rep(FALSE, nrow(dt_sr_data))]
              meth_all_mk <- NULL
              meth <- NULL
            }
            rv_view_sr$dt_st_type <- dt_st_type
            rv_view_sr$meth_all_mk <- meth_all_mk
            rv_view_sr$meth <- meth

            if(nrow(dt_st_type) != 0){
              dt_st_mk <- dt_st_type[Marker == select_mk]
              dt_st_mk <- dt_st_mk[Allele != "Multiple"]
              dt_st_mk <- dt_st_mk[!dt_st_mk$Bool_SharedAllele]
              dt_st_mk <- dt_st_mk[!dt_st_mk$Bool_AlleleDropout]
              dt_st_mk <- dt_st_mk[dt_st_mk$Bool_InfoAllele]
              rv_view_sr$dt_st_mk <- dt_st_mk
            }

            dt_kit <- fread(paste0(path_pack, "/extdata/kit/", name_kit_stutter, ".csv"))
            dt_kit <- change_data_type_kit(dt_kit)
            rv_view_sr$dt_kit_mk <- dt_kit[Marker == select_mk]
          }
        }
      })


      output$dt_average_sr_display <- renderDataTable({
        dt_average_sr <- rv_view_sr$dt_average_sr

        if(!is.null(dt_average_sr)){
          dt_average_sr_display <- dt_average_sr[, !c("Corrected allele (BSR)", "Corrected allele (FSR)", "Corrected allele (DSR)", "Corrected allele (M2SR)"), with = FALSE]

          datatable(
            dt_average_sr_display,
            selection = "none",
            options = list(dom = "t", iDisplayLength = nrow(dt_average_sr_display), ordering = FALSE, autoWidth = FALSE, scrollY = "600px", fixedHeader = TRUE),
            rownames = FALSE
          ) %>%
            formatSignif(columns = c("Average BSR", "Average FSR", "Average DSR", "Average M2SR"), digits = 4) %>%
            formatStyle(names(dt_average_sr_display), lineHeight = "80%")
        }
      })

      output$ui_plot_sr <- renderUI({
        dt_st_type <- rv_view_sr$dt_st_type
        meth_all_mk <- rv_view_sr$meth_all_mk
        meth <- rv_view_sr$meth
        dt_st_mk <- rv_view_sr$dt_st_mk
        dt_kit_mk <- rv_view_sr$dt_kit_mk
        markers_auto <- rv_view_sr$markers_auto

        if(!is.null(dt_st_type) && !is.null(meth_all_mk) && !is.null(meth) && !is.null(dt_st_mk) && !is.null(dt_kit_mk) && length(markers_auto) != 0){
          plotOutput(session$ns("plot_sr"), width = "100%", height = "300px")
        }
      })

      output$plot_sr <- renderPlot({
        dt_st_type <- rv_view_sr$dt_st_type
        meth_all_mk <- rv_view_sr$meth_all_mk
        meth <- rv_view_sr$meth
        dt_st_mk <- rv_view_sr$dt_st_mk
        dt_kit_mk <- rv_view_sr$dt_kit_mk
        markers_auto <- rv_view_sr$markers_auto

        if(!is.null(dt_st_type) && !is.null(meth_all_mk) && !is.null(meth) && !is.null(dt_st_mk) && !is.null(dt_kit_mk) && length(markers_auto) != 0){
          select_mk <- input$select_mk
          select_st_type <- input$select_st_type
          df_select_sr <- rv_view_sr$df_select_sr

          x_al <- numeric(0)
          y_sr <- numeric(0)

          if(nrow(dt_st_mk) != 0){
            # Extract (corrected) allele repeat numbers
            if(meth == "Linear (Corrected Allele)"){
              alcor <- dt_kit_mk[["Corrected Allele"]]
              x_al <- alcor[match(dt_st_mk$Allele, dt_kit_mk[["Allele"]])]
            }else{
              x_al <- as.numeric(dt_st_mk$Allele)
            }

            y_sr <- dt_st_mk$SR
          }

          x_al_multi_other_mk <- numeric(0)
          y_sr_multi_other_mk <- numeric(0)

          if(meth == "Uniform (Multiple loci)"){

            for(i in 1:length(markers_auto)){
              if(markers_auto[i] != select_mk && meth_all_mk[i] == "Uniform (Multiple loci)"){
                dt_st_other_mk <- dt_st_type[Marker == markers_auto[i]]
                dt_st_other_mk <- dt_st_other_mk[Allele != "Multiple"]
                dt_st_other_mk <- dt_st_other_mk[!dt_st_other_mk$Bool_SharedAllele]
                dt_st_other_mk <- dt_st_other_mk[!dt_st_other_mk$Bool_AlleleDropout]
                dt_st_other_mk <- dt_st_other_mk[dt_st_other_mk$Bool_InfoAllele]

                x_al_multi_other_mk <- c(x_al_multi_other_mk, as.numeric(dt_st_other_mk$Allele))
                y_sr_multi_other_mk <- c(y_sr_multi_other_mk, dt_st_other_mk$SR)
              }
            }
          }

          # Label for x-axis
          if(meth == "Linear (Corrected Allele)"){
            x_lab <- "Corrected Allele"
          }else{
            x_lab <- "Allele"
          }

          # Label for y-axis
          if(select_st_type == "Back stutter"){
            ylab = "BSR"
          }else if(select_st_type == "Forward stutter"){
            ylab = "FSR"
          }else if(select_st_type == "Double-back stutter"){
            ylab = "DSR"
          }else if(select_st_type == "Minus 2-nt stutter"){
            ylab = "M2SR"
          }else{
            ylab = "SR"
          }

          df_target_mk <- data.frame(Allele = x_al, SR = y_sr)
          color_values <- c("Selected markers" = "#00afcc")

          if(length(x_al_multi_other_mk) == 0){
            base_plot <- ggplot() +
              geom_point(data = df_target_mk, aes(Allele, SR, color = "Selected markers"))

          }else{
            df_other_mk <- data.frame(Allele = x_al_multi_other_mk, SR = y_sr_multi_other_mk)

            base_plot <- ggplot() +
              geom_point(data = df_other_mk, aes(Allele, SR, color = "Other markers")) +
              geom_point(data = df_target_mk, aes(Allele, SR, color = "Selected markers"))

            color_values <- c(color_values, "Other markers" = "#c0c0c0")
          }

          if(!is.null(df_select_sr)){
            base_plot <- base_plot + geom_point(data = df_select_sr, aes(Allele, SR, color = "Selected data"))

            color_values <- c(color_values, "Selected data" = "#e9bc00")
          }

          base_plot +
            labs(x = x_lab, y = ylab, color = NULL) +
            scale_color_manual(values = color_values) +
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
      })

      output$dt_sr_data_display <- renderDataTable({
        dt_st_mk <- rv_view_sr$dt_st_mk

        if(!is.null(dt_st_mk)){
          rv_view_sr$dt_sr_data_display <- dt_sr_data_display <- dt_st_mk[, c("SampleName", "Allele", "SR", "Peak height (stutter)", "Peak height (allele)"), with = FALSE]

          datatable(
            dt_sr_data_display,
            selection = list(mode = "single", target = "row"),
            extensions = "FixedColumns",
            options = list(dom = "t", iDisplayLength = nrow(dt_sr_data_display), autoWidth = FALSE, scrollX = TRUE, scrollY = "300px",
                           fixedHeader = TRUE, fixedColumns = list(leftColumns = 1)),
            rownames = FALSE
          ) %>%
            formatSignif(columns = c("SR"), digits = 4) %>%
            formatStyle(names(dt_sr_data_display), lineHeight = "80%")
        }
      })

      observeEvent(input$dt_sr_data_display_rows_selected, {
        index_select <- input$dt_sr_data_display_rows_selected
        if(length(index_select) == 1){
          dt_sr_data_display <- rv_view_sr$dt_sr_data_display
          al <- dt_sr_data_display[["Allele"]][index_select]
          sr <- dt_sr_data_display[["SR"]][index_select]

          meth <- rv_view_sr$meth
          if(meth == "Linear (Corrected Allele)"){
            dt_kit_mk <- rv_view_sr$dt_kit_mk
            alcor <- dt_kit_mk[["Corrected Allele"]][dt_kit_mk[["Allele"]] == al]
          }else{
            alcor <- as.numeric(al)
          }

          rv_view_sr$df_select_sr <- data.frame(Allele = alcor, SR = sr)
        }else{
          rv_view_sr$df_select_sr <- NULL
        }
      })
    }
  )
}
