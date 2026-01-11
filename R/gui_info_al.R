#' info_al_ui
#'
#' @description The function to create the ui module for information on the alleles
#' @param id The ui ID
info_al_ui <- function(id){
  ns <- NS(id)

  tabPanel("2. Register alleles",
           titlePanel("2. Register alleles"),
           br(),
           sidebarLayout(
             sidebarPanel(width = 3,
                          h4("Kit"),
                          fluidRow(
                            column(3, br(), shinyFilesButton(ns("load_kit"), label = "Load", title = "Kit", multiple = FALSE)),
                            column(9, uiOutput(ns("name_kit")))
                          ),
                          uiOutput(ns("select_marker")),
                          br(),
                          h4("Add an allele"),
                          textInput(ns("allele_add"), label = "Enter an allele", value = ""),
                          actionButton(ns("act_add_allele"), label = "Add"),
                          br(),
                          br(),
                          h4("Delete an allele"),
                          uiOutput(ns("allele_delete")),
                          actionButton(ns("act_delete_allele"), label = "Delete"),
                          br(),
                          br(),
                          br(),
                          disabled(actionButton(ns("act_save_allele"), label = "Save", class = "btn btn-primary btn-lg"))
             ),
             mainPanel(width = 9,
                       rHandsontableOutput(ns("dt_info_al"))
             )
           )
  )
}

#' info_al_server
#'
#' @description The function to create the server module for information on the alleles
#' @param id The server ID
#' @param path_pack The path of the package KongohPlus
info_al_server <- function(id, path_pack){
  moduleServer(
    id,
    function(input, output, session){
      rv_info_al <- reactiveValues()
      rv_info_al$dt_kit <- NULL
      rv_info_al$name_kit <- ""
      rv_info_al$list_dt_info_al <- NULL
      rv_info_al$markers <- NULL
      rv_info_al$rep_uni_len_all_mk <- NULL
      rv_info_al$chr_all_mk <- NULL
      rv_info_al$select_marker <- NULL
      rv_info_al$select_chr <- "Autosomal"

      shinyFileChoose(input, "load_kit", roots = c(documents = paste0(path_pack, "/extdata/kit")), session = session, filetypes = c("csv"))

      observeEvent(input$load_kit, {
        req(input$load_kit)
        path_kit <- parseFilePaths(roots = c(documents = paste0(path_pack, "/extdata/kit")), input$load_kit)$datapath
        if(length(path_kit) > 0){
          dt_kit <- fread(path_kit)
          dt_kit <- try(change_data_type_kit(dt_kit))
          if(inherits(dt_kit, "try-error")){
            showModal(modalDialog(title = "Error", "The file format of the selected kit is invalid. Please update a file in the correct format.", easyClose = TRUE, footer = NULL))
          }else{
            setindex(dt_kit, Marker)
            markers <- unique(dt_kit$Marker)
            n_mk <- length(markers)
            list_dt_info_al <- vector(mode = "list", length = n_mk)
            names(list_dt_info_al) <- markers
            rep_uni_len_all_mk <- rep(NA, n_mk)
            names(rep_uni_len_all_mk) <- markers
            chr_all_mk <- rep("", n_mk)
            names(chr_all_mk) <- markers
            for(i in 1:n_mk){
              dt_kit_mk <- dt_kit[Marker == markers[i]]
              alleles <- dt_kit_mk[["Allele"]]
              if(length(alleles) == 1 && is.na(alleles)){
                list_dt_info_al[[i]] <- data.table("Allele" = character(0),
                                                   "Corrected Allele" = numeric(0),
                                                   "Size" = numeric(0))
              }else{
                list_dt_info_al[[i]] <- data.table("Allele" = dt_kit_mk[["Allele"]],
                                                   "Corrected Allele" = dt_kit_mk[["Corrected Allele"]],
                                                   "Size" = dt_kit_mk[["Size"]])
              }
              rep_uni_len_all_mk[i] <- unique(dt_kit_mk[["Repeat Unit Length"]])[1]
              chr_all_mk[i] <- unique(dt_kit_mk[["Chromosome"]])[1]
            }
            rv_info_al$dt_kit <- dt_kit
            rv_info_al$name_kit <- str_replace(basename(path_kit), ".{4}$", "")
            rv_info_al$list_dt_info_al <- list_dt_info_al
            rv_info_al$markers <- markers
            rv_info_al$rep_uni_len_all_mk <- rep_uni_len_all_mk
            rv_info_al$chr_all_mk <- chr_all_mk

            enable("act_save_allele")
          }
        }
      })

      output$name_kit <- renderUI({
        disabled(textInput(session$ns("name_kit"), label = "", value = rv_info_al$name_kit))
      })

      output$select_marker <- renderUI({
        markers <- rv_info_al$markers
        selectInput(session$ns("select_marker"), label = h4("Marker"), choices = markers, selected = markers[1], multiple = FALSE, selectize = FALSE)
      })

      observeEvent(input$select_marker, {
        # Keep data for the marker selected one before
        dt_info_al <- hot_to_r(input$dt_info_al)
        if(!is.null(dt_info_al)){
          rv_info_al$list_dt_info_al[[rv_info_al$select_marker]] <- dt_info_al
        }

        # Update information on the newly selected marker
        select_marker <- input$select_marker
        rv_info_al$select_chr <- rv_info_al$chr_all_mk[select_marker]
        rv_info_al$select_marker <- select_marker
      })

      observeEvent(input$allele_add, {
        if(rv_info_al$name_kit != ""){
          allele_add <- input$allele_add
          select_marker <- input$select_marker
          select_chr <- rv_info_al$select_chr
          list_dt_info_al <- rv_info_al$list_dt_info_al
          dt_info_al <- list_dt_info_al[[select_marker]]
          alleles <- dt_info_al[["Allele"]]

          if(is.element(allele_add, alleles)){
            showFeedbackDanger(inputId = "allele_add", text = "The allele is already registered.")
            disable("act_add_allele")
          }else if(select_chr == "Autosomal"){
            options(warn = -1)
            allele_add <- as.numeric(allele_add)
            options(warn = 0)
            if(is.na(allele_add) || allele_add < 0){
              showFeedbackDanger(inputId = "allele_add", text = "Only a positive number is allowed.")
              disable("act_add_allele")
            }else{
              hideFeedback("allele_add")
              enable("act_add_allele")
            }
          }else{
            if(allele_add == ""){
              showFeedbackDanger(inputId = "allele_add", text = "Enter the allele name.")
              disable("act_add_allele")
            }else{
              hideFeedback("allele_add")
              enable("act_add_allele")
            }
          }
        }else{
          hideFeedback("allele_add")
          disable("act_add_allele")
        }
      })

      observeEvent(input$act_add_allele, {
        dt_info_al <- hot_to_r(input$dt_info_al)
        dt_info_al_add <- data.table("Allele" = input$allele_add,
                                     "Corrected Allele" = NA,
                                     "Size" = NA)
        dt_info_al <- rbind(dt_info_al, dt_info_al_add)

        if(rv_info_al$select_chr == "Autosomal"){
          dt_info_al <- dt_info_al[order(as.numeric(dt_info_al[["Allele"]]))]
        }else{
          dt_info_al <- dt_info_al[order(Allele)]
        }

        rv_info_al$list_dt_info_al[[input$select_marker]] <- dt_info_al
      })

      output$allele_delete <- renderUI({
        select_marker <- input$select_marker
        if(!is.null(select_marker)){
          list_dt_info_al <- rv_info_al$list_dt_info_al
          dt_info_al <- list_dt_info_al[[select_marker]]
          alleles <- dt_info_al[["Allele"]]
          if(length(alleles) == 0){
            alleles <- NULL
            select_al <- NULL
          }else{
            select_al <- alleles[1]
          }
        }else{
          alleles <- NULL
          select_al <- NULL
        }
        selectInput(session$ns("allele_delete"), label = "Select an allele", choices = alleles, selected = select_al, multiple = FALSE, selectize = FALSE)
      })

      observeEvent(input$act_delete_allele, {
        allele_delete <- input$allele_delete
        if(!is.null(allele_delete)){
          dt_info_al <- hot_to_r(input$dt_info_al)
          index_delete <- which(dt_info_al[["Allele"]] == allele_delete)
          rv_info_al$list_dt_info_al[[input$select_marker]] <- dt_info_al[-index_delete]
        }
      })

      observeEvent(input$act_save_allele, {
        rv_info_al$list_dt_info_al[[input$select_marker]] <- hot_to_r(input$dt_info_al)

        name_kit <- rv_info_al$name_kit
        list_dt_info_al <- rv_info_al$list_dt_info_al
        markers <- rv_info_al$markers
        rep_uni_len_all_mk <- rv_info_al$rep_uni_len_all_mk
        chr_all_mk <- rv_info_al$chr_all_mk
        n_mk <- length(list_dt_info_al)

        # Check whether all sizes are defined or not
        markers_no_allele <- character(0)
        markers_no_size <- character(0)
        for(i in 1:n_mk){
          mk <- markers[i]
          dt_info_al <- list_dt_info_al[[mk]]
          sizes <- dt_info_al[["Size"]]
          if(length(sizes) == 0){
            markers_no_allele <- c(markers_no_allele, mk)
          }else if(any(is.na(sizes))){
            markers_no_size <- c(markers_no_size, mk)
          }
        }

        error_message <- ""
        if(length(markers_no_allele) != 0){
          error_message <- paste0(error_message, "The allele is not defined at all in the following markers:<br>")
          for(i in 1:length(markers_no_allele)){
            error_message <- paste0(error_message, "  - ", markers_no_allele[i], "<br>")
          }
        }
        if(length(markers_no_size) != 0){
          error_message <- "The sizes of all alleles must be defined. Please check the following markers:<br>"
          for(i in 1:length(markers_no_size)){
            error_message <- paste0(error_message, "  - ", markers_no_size[i], "<br>")
          }
        }

        if(error_message != ""){
          showModal(modalDialog(title = "Error", HTML(error_message), easyClose = TRUE, footer = NULL))
        }else{
          dt_kit_new <- data.table("Marker" = character(0),
                                   "Repeat Unit Length" = integer(0),
                                   "Chromosome" = character(0),
                                   "Allele" = character(0),
                                   "Corrected Allele" = numeric(0),
                                   "Size" = numeric(0))
          for(i in 1:n_mk){
            mk <- markers[i]
            dt_info_al <- list_dt_info_al[[mk]]
            n_al <- nrow(dt_info_al)
            dt_kit_mk_left <- data.table("Marker" = rep(mk, n_al),
                                         "Repeat Unit Length" = rep(rep_uni_len_all_mk[mk], n_al),
                                         "Chromosome" = rep(chr_all_mk[mk], n_al))
            dt_kit_mk <- cbind(dt_kit_mk_left, dt_info_al)
            dt_kit_new <- rbind(dt_kit_new, dt_kit_mk)
          }

          write.csv(dt_kit_new, paste0(path_pack, "/extdata/kit/", name_kit, ".csv"), row.names = FALSE)
          showModal(modalDialog(title = "Information", "The allele information was saved successfully.", easyClose = TRUE, footer = NULL))
        }
      })

      output$dt_info_al <- renderRHandsontable({
        select_marker <- rv_info_al$select_marker
        list_dt_info_al <- rv_info_al$list_dt_info_al

        if(!is.null(select_marker) && !is.null(list_dt_info_al)){
          rhandsontable(list_dt_info_al[[select_marker]], rowHeaders = FALSE, overflow = "visible") %>%
            hot_table(fillHandle = FALSE, stretchH = "all", rowHeights = 30) %>%
            hot_col(c("Allele"), type = "text", readOnly = TRUE) %>%
            hot_col(c("Corrected Allele"), type = "numeric") %>%
            hot_col(c("Size"), type = "numeric")
        }
      })
    }
  )
}
