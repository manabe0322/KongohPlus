#' size_ui
#'
#' @description The function to create the ui module for the allele sizes
size_ui <- function(id){
  ns <- NS(id)

  tabPanel("Allele size",
           titlePanel("Allele size"),
           br(),
           sidebarLayout(
             sidebarPanel(width = 3,
                          uiOutput(ns("select_kit_size")),
                          br(),
                          uiOutput(ns("select_marker_size")),
                          br(),
                          actionButton(ns("act_size_edit"), label = "Edit"),
                          br(),
                          br(),
                          actionButton(ns("act_size_add"), label = "Add"),
                          br(),
                          br(),
                          actionButton(ns("act_size_delete"), label = "Delete")
             ),
             mainPanel(width = 9,
                       dataTableOutput(ns("dt_info_size"))
             )
           )
  )
}

#' create_init_dt_info_size
#'
#' @description The function to create the initial data.table for information on the allele sizes
create_init_dt_info_size <- function(path_pack){
  kit_all <- list.files(paste0(path_pack, "/extdata/kit"))
  if(length(kit_all) != 0){
    dt_kit <- fread(paste0(path_pack, "/extdata/kit/", kit_all[1]))
    dt_kit <- change_data_type_kit(dt_kit)
    marker <- dt_kit$Marker[1]
    dt_info_size <- create_dt_info_size(dt_kit, marker)
  }else{
    dt_info_size <- NULL
  }
  return(dt_info_size)
}

#' create_dt_info_size
#'
#' @description The function to create the data.table for information on the allele sizes
create_dt_info_size <- function(dt_kit, marker){
  index_mk <- which(dt_kit$Marker == marker)
  dt_info_size <- data.table(Marker = dt_kit$Marker[index_mk],
                             Repeat_unit_length = dt_kit$Repeat_unit_length[index_mk],
                             Chr = dt_kit$Chr[index_mk],
                             Allele = dt_kit$Allele[index_mk],
                             Size = dt_kit$Size[index_mk])
  return(dt_info_size)
}

#' load_first_dt_kit
#'
#' @description The function to load the initial dt_kit
#' @param path_pack The package path
load_first_dt_kit <- function(path_pack){
  kit_all <- list.files(paste0(path_pack, "/extdata/kit"))
  if(length(kit_all) != 0){
    dt_kit <- fread(paste0(path_pack, "/extdata/kit/", kit_all[1]))
  }else{
    dt_kit <- NULL
  }
  return(dt_kit)
}

#' size_server
#'
#' @description The function to create the server module for the allele sizes
size_server <- function(id, path_pack){
  moduleServer(
    id,
    function(input, output, session){

      ######################################
      # Define the initial reactive values #
      ######################################

      rv_size <- reactiveValues()
      rv_size$dt_info_size <- create_init_dt_info_size(path_pack)
      rv_size$dt_kit <- load_first_dt_kit(path_pack)

      #################################################################
      # Define selectInput and update information on the allele sizes #
      #################################################################

      output$select_kit_size <- renderUI({
        kit_all <- list.files(paste0(path_pack, "/extdata/kit"))
        kit_all <- gsub(".{4}$", "", kit_all)
        selectInput(session$ns("select_kit_size"), label = h4("Kit"), choices = kit_all, selected = kit_all[1], selectize = FALSE)
      })

      observeEvent(input$select_kit_size, {
        dt_kit <- fread(paste0(path_pack, "/extdata/kit/", input$select_kit_size, ".csv"))
        rv_size$dt_kit <- dt_kit <- change_data_type_kit(dt_kit)
        markers <- unique(dt_kit$Marker)

        output$select_marker_size <- renderUI({
          selectInput(session$ns("select_marker_size"), label = h4("Marker"), choices = markers, selected = markers[1], selectize = FALSE)
        })
      })

      observeEvent(input$select_marker_size, {
        rv_size$dt_info_size <- create_dt_info_size(rv_size$dt_kit, input$select_marker_size)
      })

      #######################
      # Edit an allele size #
      #######################

      observeEvent(input$act_size_edit, {
        pos_select <- input$dt_info_size_rows_selected
        if(isTruthy(pos_select)){
          dt_kit <- rv_size$dt_kit
          dt_info_size <- rv_size$dt_info_size
          kit <- input$select_kit_size
          marker <- input$select_marker_size
          chr <- dt_info_size$Chr[1]
          allele <- dt_info_size$Allele[pos_select]

          showModal(modalDialog(
            sidebarLayout(
              sidebarPanel(width = 3,
                           h4("Kit"),
                           p(kit),
                           br(),
                           h4("Marker"),
                           p(marker),
                           br(),
                           h4("Chromosome"),
                           p(chr),
                           br(),
                           h4("Allele"),
                           p(allele)
              ),
              mainPanel(width = 9,
                        numericInput(session$ns("size_edit"), label = h4("Allele size"), value = dt_info_size$Size[pos_select])
              )
            ),
            footer = tagList(
              actionButton(session$ns("act_edit_size_save"), "Save"),
              modalButton("Cancel")
            ),
            size = "l"
          ))
        }else{
          showModal(modalDialog(title = "Error", "Select a row!", easyClose = TRUE, footer = NULL))
        }
      })

      observeEvent(input$act_edit_size_save, {
        pos_select <- input$dt_info_size_rows_selected
        dt_kit <- rv_size$dt_kit
        dt_info_size <- rv_size$dt_info_size
        kit <- input$select_kit_size
        marker <- input$select_marker_size
        allele <- dt_info_size$Allele[pos_select]
        size_edit <- input$size_edit

        # Update information on allele sizes
        dt_info_size$Size[pos_select] <- size_edit
        rv_size$dt_info_size <- dt_info_size

        # Update dt_kit
        index_target <- intersect(which(dt_kit$Marker == marker), which(dt_kit$Allele == allele))
        dt_kit$Size[index_target] <- size_edit
        rv_size$dt_kit <- dt_kit
        write.csv(dt_kit, paste0(path_pack, "/extdata/kit/", kit, ".csv"), row.names = FALSE)

        removeModal()
      }, ignoreInit = TRUE)

      ######################
      # Add an allele size #
      ######################

      observeEvent(input$act_size_add, {
        dt_info_size <- rv_size$dt_info_size
        kit <- input$select_kit_size
        marker <- input$select_marker_size
        chr <- dt_info_size$Chr[1]

        showModal(modalDialog(
          sidebarLayout(
            sidebarPanel(width = 3,
                         h4("Kit"),
                         p(kit),
                         br(),
                         h4("Marker"),
                         p(marker),
                         br(),
                         h4("Chromosome"),
                         p(chr)
            ),
            mainPanel(width = 9,
                      textInput(session$ns("allele_add"), label = h4("Allele"), value = ""),
                      br(),
                      numericInput(session$ns("size_add"), label = h4("Allele size"), value = NULL)
            )
          ),
          footer = tagList(
            disabled(actionButton(session$ns("act_add_size_save"), "Save")),
            modalButton("Cancel")
          ),
          size = "l"
        ))
      })

      observeEvent(input$allele_add, {
        dt_info_size <- rv_size$dt_info_size
        alleles <- dt_info_size$Allele
        chr <- dt_info_size$Chr[1]
        allele_add <- input$allele_add
        hideFeedback("allele_add")
        options(warn = -1)
        if(allele_add == ""){
          showFeedbackDanger(inputId = "allele_add", text = "Enter the allele name in this field.")
          disable("act_add_size_save")
        }else if(is.element(allele_add, alleles)){
          showFeedbackDanger(inputId = "allele_add", text = "This allele has been already registered.")
          disable("act_add_size_save")
        }else if(chr == "Auto" && is.na(as.numeric(allele_add))){
          showFeedbackDanger(inputId = "allele_add", text = "Only a number is allowed for the autosomal marker.")
          disable("act_add_size_save")
        }else{
          enable("act_add_size_save")
        }
        options(warn = 0)
      })

      observeEvent(input$act_add_size_save, {
        dt_kit <- rv_size$dt_kit
        dt_info_size <- rv_size$dt_info_size
        kit <- input$select_kit_size
        marker <- input$select_marker_size
        chr <- dt_info_size$Chr[1]
        alleles <- dt_info_size$Allele
        allele_add <- input$allele_add
        size_add <- input$size_add

        info_marker <- dt_info_size$Marker
        info_rep <- dt_info_size$Repeat_unit_length
        info_chr <- dt_info_size$Chr
        info_allele <- dt_info_size$Allele
        info_size <- dt_info_size$Size

        if(is.na(info_allele[1]) || is.na(info_size[1])){
          alleles <- allele_add
          sizes <- size_add
        }else{
          alleles <- c(info_allele, allele_add)
          sizes <- c(info_size, size_add)
        }

        n_al <- length(alleles)
        if(chr == "Auto"){
          # Sort alleles
          alleles <- as.numeric(alleles)
          index_sort <- order(alleles)
          alleles <- alleles[index_sort]
          sizes <- sizes[index_sort]

          # Change data type of alleles (numeric to character)
          alleles_char <- rep("", n_al)
          for(i in 1:n_al){
            al_1 <- alleles[i]
            if(al_1 %% 1 == 0){
              alleles_char[i] <- as.character(round(al_1, 0))
            }else{
              alleles_char[i] <- as.character(round(al_1, 1))
            }
          }
        }else{
          index_sort <- order(alleles)
          alleles_char <- alleles[index_sort]
          sizes <- sizes[index_sort]
        }

        # Update information on allele sizes
        if(n_al == 1){
          dt_info_size_add <- data.table(Marker = info_marker,
                                         Repeat_unit_length = info_rep,
                                         Chr = info_chr,
                                         Allele = alleles_char,
                                         Size = sizes)
        }else{
          dt_info_size_add <- data.table(Marker = c(info_marker, info_marker[1]),
                                         Repeat_unit_length = c(info_rep, info_rep[1]),
                                         Chr = c(info_chr, info_chr[1]),
                                         Allele = alleles_char,
                                         Size = sizes)
        }
        rv_size$dt_info_size <- dt_info_size_add

        # Update dt_kit
        index_target <- which(dt_kit$Marker == marker)
        if(index_target[1] == 1){
          dt_kit_other_mk <- dt_kit[dt_kit$Marker != marker, ]
          dt_kit_add <- rbind(dt_info_size_add, dt_kit_other_mk)
        }else if(index_target[length(index_target)] == nrow(dt_kit)){
          dt_kit_other_mk <- dt_kit[dt_kit$Marker != marker, ]
          dt_kit_add <- rbind(dt_kit_other_mk, dt_info_size_add)
        }else{
          dt_kit_other_mk_1 <- dt_kit[1:(index_target[1] - 1), ]
          dt_kit_other_mk_2 <- dt_kit[(index_target[length(index_target)] + 1):nrow(dt_kit), ]
          dt_kit_add <- rbind(dt_kit_other_mk_1, dt_info_size_add, dt_kit_other_mk_2)
        }
        write.csv(dt_kit_add, paste0(path_pack, "/extdata/kit/", kit, ".csv"), row.names = FALSE)
        rv_size$dt_kit <- dt_kit_add

        removeModal()
      }, ignoreInit = TRUE)

      #########################
      # Delete an allele size #
      #########################

      observeEvent(input$act_size_delete, {
        pos_select <- input$dt_info_size_rows_selected
        if(isTruthy(pos_select)){
          showModal(modalDialog(
            title = "Delete an allele size",
            "Click Ok button to delete the selected allele size.",
            footer = tagList(
              actionButton(session$ns("act_delete_size_ok"), "OK"),
              modalButton("Cancel")
            ),
            size = "l"
          ))
        }else{
          showModal(modalDialog(title = "Error", "Select a row!", easyClose = TRUE, footer = NULL))
        }
      })

      observeEvent(input$act_delete_size_ok, {
        pos_select <- input$dt_info_size_rows_selected
        dt_kit <- rv_size$dt_kit
        dt_info_size <- rv_size$dt_info_size
        kit <- input$select_kit_size
        marker <- input$select_marker_size

        info_marker <- dt_info_size$Marker
        info_rep <- dt_info_size$Repeat_unit_length
        info_chr <- dt_info_size$Chr
        info_allele <- dt_info_size$Allele
        info_size <- dt_info_size$Size

        if(length(info_marker) == 1){
          dt_info_size_delete <- data.table(Marker = info_marker,
                                            Repeat_unit_length = info_rep,
                                            Chr = info_chr,
                                            Allele = "",
                                            Size = NA)
        }else{
          dt_info_size_delete <- data.table(Marker = info_marker[- pos_select],
                                            Repeat_unit_length = info_rep[- pos_select],
                                            Chr = info_chr[- pos_select],
                                            Allele = info_allele[- pos_select],
                                            Size = info_size[- pos_select])
        }
        rv_size$dt_info_size <- dt_info_size_delete

        # Update dt_kit
        index_target <- which(dt_kit$Marker == marker)
        if(index_target[1] == 1){
          dt_kit_other_mk <- dt_kit[dt_kit$Marker != marker, ]
          dt_kit_delete <- rbind(dt_info_size_delete, dt_kit_other_mk)
        }else if(index_target[length(index_target)] == nrow(dt_kit)){
          dt_kit_other_mk <- dt_kit[dt_kit$Marker != marker, ]
          dt_kit_delete <- rbind(dt_kit_other_mk, dt_info_size_delete)
        }else{
          dt_kit_other_mk_1 <- dt_kit[1:(index_target[1] - 1), ]
          dt_kit_other_mk_2 <- dt_kit[(index_target[length(index_target)] + 1):nrow(dt_kit), ]
          dt_kit_delete <- rbind(dt_kit_other_mk_1, dt_info_size_delete, dt_kit_other_mk_2)
        }
        write.csv(dt_kit_delete, paste0(path_pack, "/extdata/kit/", kit, ".csv"), row.names = FALSE)
        rv_size$dt_kit <- dt_kit_delete

        removeModal()
      }, ignoreInit = TRUE)

      #######################################
      # Display information on allele sizes #
      #######################################

      output$dt_info_size <- renderDataTable({
        dt_info_size <- rv_size$dt_info_size
        if(!is.null(dt_info_size)){
          datatable(
            dt_info_size,
            colnames = c("Marker", "Length of the repeat unit", "Chromosome", "Allele", "Size"),
            selection = list(mode = "single", target = "row"),
            options = list(dom = "t", iDisplayLength = nrow(dt_info_size), ordering = FALSE, autoWidth = TRUE),
            rownames = FALSE
          )
        }
      })
    }
  )
}
