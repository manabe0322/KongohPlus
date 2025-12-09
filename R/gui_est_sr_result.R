#' est_sr_result_ui
#'
#' @description The function to create the ui module to create the result tab for estimating average stutter ratios
est_sr_result_ui <- function(id){
  ns <- NS(id)

  tabPanel("Result",
           br(),
           fluidRow(
             column(3,
                    wellPanel(
                      uiOutput(ns("select_mk")),
                      uiOutput(ns("select_st_type"))
                    ),
                    br(),
                    actionButton(ns("act_saveas_est_sr_1"), label = "Save as", class = "btn btn-primary btn-lg")
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

#' est_sr_setting_server
#'
#' @description The function to create the server module to create the result tab for estimating average stutter ratios
est_sr_result_server <- function(id, path_pack, rv_est_sr){
  moduleServer(
    id,
    function(input, output, session){
      rv_est_sr_result <- reactiveValues()
      rv_est_sr_result$markers_auto <- NULL
      rv_est_sr_result$list_dt_average_sr <- NULL
      rv_est_sr_result$dt_setting_est_sr <- NULL
      rv_est_sr_result$dt_sr_data <- NULL
      rv_est_sr_result$name_kit_stutter <- NULL
      rv_est_sr_result$dt_average_sr <- NULL
      rv_est_sr_result$dt_st_type <- NULL
      rv_est_sr_result$meth_all_mk <- NULL
      rv_est_sr_result$meth <- NULL
      rv_est_sr_result$dt_st_mk <- NULL
      rv_est_sr_result$dt_kit_mk <- NULL
      rv_est_sr_result$dt_sr_data_display <- NULL
      rv_est_sr_result$df_select_sr <- NULL

      observe({
        req(rv_est_sr)
        data_est_sr <- rv_est_sr$data_est_sr

        if(!is.null(data_est_sr)){
          rv_est_sr_result$markers_auto <- names(data_est_sr$list_dt_average_sr)
          rv_est_sr_result$list_dt_average_sr <- data_est_sr$list_dt_average_sr
          rv_est_sr_result$dt_setting_est_sr <- data_est_sr$dt_setting_est_sr
          rv_est_sr_result$dt_sr_data <- data_est_sr$dt_sr_data
          rv_est_sr_result$name_kit_stutter <- data_est_sr$name_kit_stutter
          rv_est_sr_result$dt_average_sr <- NULL
          rv_est_sr_result$dt_st_type <- NULL
          rv_est_sr_result$meth_all_mk <- NULL
          rv_est_sr_result$meth <- NULL
          rv_est_sr_result$dt_st_mk <- NULL
          rv_est_sr_result$dt_kit_mk <- NULL
          rv_est_sr_result$dt_sr_data_display <- NULL
          rv_est_sr_result$df_select_sr <- NULL
          enable("act_saveas_est_sr_1")
        }else{
          rv_est_sr_result$markers_auto <- NULL
          rv_est_sr_result$list_dt_average_sr <- NULL
          rv_est_sr_result$dt_setting_est_sr <- NULL
          rv_est_sr_result$dt_sr_data <- NULL
          rv_est_sr_result$name_kit_stutter <- NULL
          rv_est_sr_result$dt_average_sr <- NULL
          rv_est_sr_result$dt_st_type <- NULL
          rv_est_sr_result$meth_all_mk <- NULL
          rv_est_sr_result$meth <- NULL
          rv_est_sr_result$dt_st_mk <- NULL
          rv_est_sr_result$dt_kit_mk <- NULL
          rv_est_sr_result$dt_sr_data_display <- NULL
          rv_est_sr_result$df_select_sr <- NULL
          disable("act_saveas_est_sr_1")
        }
      })

      output$select_mk <- renderUI({
        req(rv_est_sr_result$markers_auto)
        markers_auto <- rv_est_sr_result$markers_auto
        selectInput(session$ns("select_mk"), label = h4("Marker"), choices = markers_auto, selected = markers_auto[1], multiple = FALSE, selectize = FALSE)
      })

      output$select_st_type <- renderUI({
        req(rv_est_sr_result$markers_auto)
        selectInput(session$ns("select_st_type"), label = h4("Stutter type"), choices = c("Back stutter", "Forward stutter", "Double-back stutter", "Minus 2-nt stutter"), selected = "Back stutter", multiple = FALSE, selectize = FALSE)
      })

      observeEvent(input$act_saveas_est_sr_1, {
        data_est_sr <- rv_est_sr$data_est_sr

        if(!is.null(data_est_sr)){
          showModal(modalDialog(
            title = "Save stutter data",
            textInput(session$ns("fn_est_sr"), label = "Enter the file name"),
            footer = tagList(
              disabled(actionButton(session$ns("act_saveas_est_sr_2"), "Save")),
              modalButton("Close")
            )
          ))
        }

        else{
          showModal(modalDialog(title = "Error", "No stutter data.", easyClose = TRUE, footer = NULL))
        }
      })

      observeEvent(input$fn_est_sr, {
        fn_est_sr <- input$fn_est_sr
        sr_file_names <- list.files(paste0(path_pack, "/extdata/stutter"))
        if(is.element(paste0(fn_est_sr, ".RData"), sr_file_names)){
          showFeedbackDanger(inputId = "fn_est_sr", text = "A file with the same name already exists.")
          disable("act_saveas_est_sr_2")
        }else if(nchar(fn_est_sr) == 0){
          hideFeedback("fn_est_sr")
          disable("act_saveas_est_sr_2")
        }else{
          hideFeedback("fn_est_sr")
          enable("act_saveas_est_sr_2")
        }
      })

      observeEvent(input$act_saveas_est_sr_2, {
        req(rv_est_sr)
        data_est_sr <- rv_est_sr$data_est_sr
        fn_est_sr <- input$fn_est_sr

        if(!is.null(data_est_sr)){
          save(list = c("data_est_sr"), file = paste0(path_pack, "/extdata/stutter/", fn_est_sr, ".RData"))
          removeModal()
        }

        else{
          removeModal()
          showModal(modalDialog(title = "Error", "No stutter data.", easyClose = TRUE, footer = NULL))
        }
      })

      observeEvent(list(input$select_mk, input$select_st_type), {
        req(rv_est_sr_result$list_dt_average_sr)
        req(rv_est_sr_result$dt_setting_est_sr)
        req(rv_est_sr_result$dt_sr_data)
        req(rv_est_sr_result$name_kit_stutter)
        list_dt_average_sr <- rv_est_sr_result$list_dt_average_sr
        dt_setting_est_sr <- rv_est_sr_result$dt_setting_est_sr
        dt_sr_data <- rv_est_sr_result$dt_sr_data
        name_kit_stutter <- rv_est_sr_result$name_kit_stutter

        rv_est_sr_result$dt_average_sr <- NULL
        rv_est_sr_result$dt_st_type <- NULL
        rv_est_sr_result$meth_all_mk <- NULL
        rv_est_sr_result$meth <- NULL
        rv_est_sr_result$dt_st_mk <- NULL
        rv_est_sr_result$dt_kit_mk <- NULL
        rv_est_sr_result$df_select_sr <- NULL

        select_mk <- input$select_mk
        select_st_type <- input$select_st_type

        rv_est_sr_result$dt_average_sr <- list_dt_average_sr[[select_mk]]

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
        rv_est_sr_result$dt_st_type <- dt_st_type
        rv_est_sr_result$meth_all_mk <- meth_all_mk
        rv_est_sr_result$meth <- meth

        if(nrow(dt_st_type) != 0){
          dt_st_mk <- dt_st_type[Marker == select_mk]
          dt_st_mk <- dt_st_mk[Allele != "Multiple"]
          dt_st_mk <- dt_st_mk[!dt_st_mk$Bool_SharedAllele]
          dt_st_mk <- dt_st_mk[!dt_st_mk$Bool_AlleleDropout]
          dt_st_mk <- dt_st_mk[dt_st_mk$Bool_InfoAllele]
          rv_est_sr_result$dt_st_mk <- dt_st_mk
        }

        dt_kit <- fread(paste0(path_pack, "/extdata/kit/", name_kit_stutter, ".csv"))
        dt_kit <- try(change_data_type_kit(dt_kit))
        if(inherits(dt_kit, "try-error")){
          rv_est_sr_result$dt_kit_mk <- NULL
          showModal(modalDialog(title = "Error", "The file format of the selected kit is invalid.", easyClose = TRUE, footer = NULL))
        }else{
          rv_est_sr_result$dt_kit_mk <- dt_kit[Marker == select_mk]
        }
      })

      output$dt_average_sr_display <- renderDataTable({
        req(rv_est_sr_result$dt_average_sr)
        dt_average_sr <- rv_est_sr_result$dt_average_sr

        dt_average_sr_display <- dt_average_sr[, !c("Corrected allele (BSR)", "Corrected allele (FSR)", "Corrected allele (DSR)", "Corrected allele (M2SR)"), with = FALSE]

        datatable(
          dt_average_sr_display,
          selection = "none",
          options = list(dom = "t", iDisplayLength = nrow(dt_average_sr_display), ordering = FALSE, autoWidth = TRUE, scrollY = "600px", fixedHeader = TRUE),
          rownames = FALSE
        ) %>%
          formatSignif(columns = c("Average BSR", "Average FSR", "Average DSR", "Average M2SR"), digits = 4) %>%
          formatStyle(names(dt_average_sr_display), lineHeight = "80%")
      })

      output$ui_plot_sr <- renderUI({
        req(rv_est_sr_result$dt_st_type)
        req(rv_est_sr_result$meth_all_mk)
        req(rv_est_sr_result$meth)
        req(rv_est_sr_result$dt_st_mk)
        req(rv_est_sr_result$dt_kit_mk)
        req(rv_est_sr_result$markers_auto)

        plotOutput(session$ns("plot_sr"), width = "100%", height = "300px")
      })

      output$plot_sr <- renderPlot({
        req(rv_est_sr_result$dt_st_type)
        req(rv_est_sr_result$meth_all_mk)
        req(rv_est_sr_result$meth)
        req(rv_est_sr_result$dt_st_mk)
        req(rv_est_sr_result$dt_kit_mk)
        req(rv_est_sr_result$markers_auto)
        dt_st_type <- rv_est_sr_result$dt_st_type
        meth_all_mk <- rv_est_sr_result$meth_all_mk
        meth <- rv_est_sr_result$meth
        dt_st_mk <- rv_est_sr_result$dt_st_mk
        dt_kit_mk <- rv_est_sr_result$dt_kit_mk
        markers_auto <- rv_est_sr_result$markers_auto
        df_select_sr <- rv_est_sr_result$df_select_sr  # It can be either null or not.

        select_mk <- input$select_mk
        select_st_type <- input$select_st_type

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
      })

      output$dt_sr_data_display <- renderDataTable({
        req(rv_est_sr_result$dt_st_mk)
        dt_st_mk <- rv_est_sr_result$dt_st_mk

        dt_sr_data_display <- dt_st_mk[, c("SampleName", "Allele", "SR", "Peak height (stutter)", "Peak height (allele)"), with = FALSE]
        dt_sr_data_display[["SampleName"]] <- factor(dt_sr_data_display[["SampleName"]], levels = unique(dt_sr_data_display[["SampleName"]]), labels = unique(dt_sr_data_display[["SampleName"]]))
        tmp_alleles <- as.character(sort(as.numeric(unique(dt_sr_data_display[["Allele"]]))))
        dt_sr_data_display[["Allele"]] <- factor(dt_sr_data_display[["Allele"]], levels = tmp_alleles, labels = tmp_alleles)
        rv_est_sr_result$dt_sr_data_display <- dt_sr_data_display

        datatable(
          dt_sr_data_display,
          filter = "top",
          selection = list(mode = "single", target = "row"),
          options = list(
            iDisplayLength = 10,
            autoWidth = FALSE,
            scrollX = TRUE,
            columnDefs = list(
              list(
                targets = 0,
                width = '250px'
              ),
              list(
                targets = 1,
                width = '75px'
              ),
              list(
                targets = c(2, 3, 4),
                searchable = FALSE
              )
            )
          ),
          rownames = FALSE
        ) %>%
          formatSignif(columns = c("SR"), digits = 4) %>%
          formatStyle(names(dt_sr_data_display), lineHeight = "80%")
      })

      observeEvent(input$dt_sr_data_display_rows_selected, {
        req(rv_est_sr_result$dt_kit_mk)
        req(rv_est_sr_result$dt_sr_data_display)
        req(rv_est_sr_result$meth)
        dt_kit_mk <- rv_est_sr_result$dt_kit_mk
        dt_sr_data_display <- rv_est_sr_result$dt_sr_data_display
        meth <- rv_est_sr_result$meth
        index_select <- input$dt_sr_data_display_rows_selected

        if(length(index_select) == 1){
          al <- as.character(dt_sr_data_display[["Allele"]][index_select])
          sr <- dt_sr_data_display[["SR"]][index_select]

          if(meth == "Linear (Corrected Allele)"){
            alcor <- dt_kit_mk[["Corrected Allele"]][dt_kit_mk[["Allele"]] == al]
          }else{
            alcor <- as.numeric(al)
          }

          rv_est_sr_result$df_select_sr <- data.frame(Allele = alcor, SR = sr)
        }else{
          rv_est_sr_result$df_select_sr <- NULL
        }
      })
    }
  )
}
