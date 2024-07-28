#' sr_ui
#'
#' @description The function to create the ui module for the stutter ratio
sr_ui <- function(id){
  ns <- NS(id)

  tabPanel("Average stutter ratio",
           titlePanel("Average stutter ratio"),
           br(),
           sidebarLayout(
             sidebarPanel(width = 2,
                          uiOutput(ns("select_param_sr")),
                          br(),
                          uiOutput(ns("select_marker_sr")),
                          br(),
                          actionButton(ns("act_sr_edit"), "Edit"),
                          br(),
                          br(),
                          actionButton(ns("act_sr_add"), "Add"),
                          br(),
                          br(),
                          actionButton(ns("act_sr_delete"), "Delete")
             ),
             mainPanel(width = 10,
                       dataTableOutput(ns("dt_info_sr"))
             )
           )
  )
}

#' create_init_dt_info_sr
#'
#' @description The function to create the initial data.table for information on the stutter ratio
create_init_dt_info_sr <- function(path_pack){
  param_all <- list.files(paste0(path_pack, "/extdata/param"))
  if(length(param_all) != 0){
    load(paste0(path_pack, "/extdata/param/", param_all[1]))
    marker <- dt_sr$Marker[1]
    dt_info_sr <- create_dt_info_sr(dt_sr, marker)
  }else{
    dt_info_sr <- NULL
  }
  return(dt_info_sr)
}

#' create_dt_info_sr
#'
#' @description The function to create the data.table for information on the stutter ratio
create_dt_info_sr <- function(dt_sr, marker){
  index_mk <- which(dt_sr$Marker == marker)
  dt_info_sr <- data.table(Allele = dt_sr$Allele[index_mk],
                           Average_BSR = dt_sr$Average_BSR[index_mk],
                           Average_FSR = dt_sr$Average_FSR[index_mk],
                           Average_DSR = dt_sr$Average_DSR[index_mk],
                           Average_M2SR = dt_sr$Average_M2SR[index_mk])

  # Change data type
  col_char <- c("Allele")
  dt_info_sr[, (col_char) := lapply(.SD, as.character), .SDcols = col_char]

  return(dt_info_sr)
}

#' sr_server
#'
#' @description The function to create the server module for the stutter ratio
sr_server <- function(id, path_pack){
  moduleServer(
    id,
    function(input, output, session){

      ######################################
      # Define the initial reactive values #
      ######################################

      rv_sr <- reactiveValues()
      rv_sr$dt_info_sr <- create_init_dt_info_sr(path_pack)

      ##################################################################
      # Define selectInput and update information on the stutter ratio #
      ##################################################################

      output$select_param_sr <- renderUI({
        param_all <- list.files(paste0(path_pack, "/extdata/param"))
        param_all <- gsub(".{6}$", "", param_all)
        selectInput(session$ns("select_param_sr"), label = h4("Parameter data"), choices = param_all, selected = param_all[1], selectize = FALSE)
      })

      observeEvent(input$select_param_sr, {
        load(paste0(path_pack, "/extdata/param/", input$select_param_sr, ".Rdata"))
        markers_sr <- unique(dt_sr$Marker)

        output$select_marker_sr <- renderUI({
          selectInput(session$ns("select_marker_sr"), label = h4("Marker"), choices = markers_sr, selected = markers_sr[1], selectize = FALSE)
        })
      })

      observeEvent(input$select_marker_sr, {
        load(paste0(path_pack, "/extdata/param/", input$select_param_sr, ".Rdata"))
        rv_sr$dt_info_sr <- create_dt_info_sr(dt_sr, input$select_marker_sr)
      })

      #######################################
      # Edit the stutter ratio of an allele #
      #######################################

      observeEvent(input$act_sr_edit, {
        pos_select <- input$dt_info_sr_rows_selected

        if(isTruthy(pos_select)){
          dt_info_sr <- rv_sr$dt_info_sr
          param <- input$select_param_sr
          marker <- input$select_marker_sr
          allele <- dt_info_sr$Allele[pos_select]

          showModal(modalDialog(
            sidebarLayout(
              sidebarPanel(width = 3,
                           h4("Parameter"),
                           p(param),
                           br(),
                           h4("Marker"),
                           p(marker),
                           br(),
                           h4("Allele"),
                           p(allele)
              ),
              mainPanel(width = 9,
                        numericInput(session$ns("bsr_edit"), label = h4("Average back stutter ratio"), value = dt_info_sr$Average_BSR[pos_select]),
                        br(),
                        numericInput(session$ns("fsr_edit"), label = h4("Average forward stutter ratio"), value = dt_info_sr$Average_FSR[pos_select]),
                        br(),
                        numericInput(session$ns("dsr_edit"), label = h4("Average double-back stutter ratio"), value = dt_info_sr$Average_DSR[pos_select]),
                        br(),
                        numericInput(session$ns("m2sr_edit"), label = h4("Average minus 2-nt stutter ratio"), value = dt_info_sr$Average_M2SR[pos_select])
              )
            ),
            footer = tagList(
              actionButton(session$ns("act_edit_sr_save"), "Save"),
              modalButton("Cancel")
            ),
            size = "l"
          ))
        }else{
          showModal(modalDialog(title = "Error", "Select a row!", easyClose = TRUE, footer = NULL))
        }
      })

      observeEvent(input$act_edit_sr_save, {
        pos_select <- input$dt_info_sr_rows_selected
        dt_info_sr <- rv_sr$dt_info_sr
        param <- input$select_param_sr
        marker <- input$select_marker_sr
        allele <- dt_info_sr$Allele[pos_select]
        bsr_edit <- input$bsr_edit
        fsr_edit <- input$fsr_edit
        dsr_edit <- input$dsr_edit
        m2sr_edit <- input$m2sr_edit

        # Update information on the stutter ratio
        dt_info_sr$Average_BSR[pos_select] <- bsr_edit
        dt_info_sr$Average_FSR[pos_select] <- fsr_edit
        dt_info_sr$Average_DSR[pos_select] <- dsr_edit
        dt_info_sr$Average_M2SR[pos_select] <- m2sr_edit
        rv_sr$dt_info_sr <- dt_info_sr

        # Update dt_sr
        load(paste0(path_pack, "/extdata/param/", param, ".Rdata"))
        index_target <- intersect(which(dt_sr$Marker == marker), which(dt_sr$Allele == allele))
        dt_sr$Average_BSR[index_target] <- bsr_edit
        dt_sr$Average_FSR[index_target] <- fsr_edit
        dt_sr$Average_DSR[index_target] <- dsr_edit
        dt_sr$Average_M2SR[index_target] <- m2sr_edit
        save(list = c("dt_par", "dt_sr"), file = paste0(path_pack, "/extdata/param/", param, ".Rdata"))

        removeModal()
      }, ignoreInit = TRUE)

      ######################################
      # Add the stutter ratio of an allele #
      ######################################

      observeEvent(input$act_sr_add, {
        dt_info_sr <- rv_sr$dt_info_sr
        param <- input$select_param_sr
        marker <- input$select_marker_sr

        showModal(modalDialog(
          sidebarLayout(
            sidebarPanel(width = 3,
                         h4("Parameter"),
                         p(param),
                         br(),
                         h4("Marker"),
                         p(marker)
            ),
            mainPanel(width = 9,
                      numericInput(session$ns("allele_add"), label = h4("Allele"), value = NULL),
                      br(),
                      numericInput(session$ns("bsr_add"), label = h4("Average back stutter ratio"), value = NULL),
                      br(),
                      numericInput(session$ns("fsr_add"), label = h4("Average forward stutter ratio"), value = NULL),
                      br(),
                      numericInput(session$ns("dsr_add"), label = h4("Average double-back stutter ratio"), value = NULL),
                      br(),
                      numericInput(session$ns("m2sr_add"), label = h4("Average minus 2-nt stutter ratio"), value = NULL)
            )
          ),
          footer = tagList(
            actionButton(session$ns("act_add_sr_save"), "Save"),
            modalButton("Cancel")
          ),
          size = "l"
        ))
      })

      observeEvent(input$act_add_sr_save, {
        dt_info_sr <- rv_sr$dt_info_sr
        param <- input$select_param_sr
        marker <- input$select_marker_sr
        allele_add <- input$allele_add
        bsr_add <- input$bsr_add
        fsr_add <- input$fsr_add
        dsr_add <- input$dsr_add
        m2sr_add <- input$m2sr_add

        info_allele_add <- c(dt_info_sr$Allele, allele_add)
        info_bsr_add <- c(dt_info_sr$Average_BSR, bsr_add)
        info_fsr_add <- c(dt_info_sr$Average_FSR, fsr_add)
        info_dsr_add <- c(dt_info_sr$Average_DSR, dsr_add)
        info_m2sr_add <- c(dt_info_sr$Average_M2SR, m2sr_add)

        # Sort alleles
        alleles <- as.numeric(info_allele_add)
        index_sort <- order(alleles)
        alleles <- alleles[index_sort]
        info_bsr_add <- info_bsr_add[index_sort]
        info_fsr_add <- info_fsr_add[index_sort]
        info_dsr_add <- info_dsr_add[index_sort]
        info_m2sr_add <- info_m2sr_add[index_sort]

        # Change data type of alleles (numeric to character)
        n_al <- length(alleles)
        alleles_char <- rep("", n_al)
        for(i in 1:n_al){
          al_1 <- alleles[i]
          if(al_1 %% 1 == 0){
            alleles_char[i] <- as.character(round(al_1, 0))
          }else{
            alleles_char[i] <- as.character(round(al_1, 1))
          }
        }

        # Update information on the stutter ratio
        dt_info_sr_add <- data.table(Allele = alleles_char,
                                     Average_BSR = info_bsr_add,
                                     Average_FSR = info_fsr_add,
                                     Average_DSR = info_dsr_add,
                                     Average_M2SR = info_m2sr_add)
        rv_sr$dt_info_sr <- dt_info_sr_add

        # Update dt_sr
        load(paste0(path_pack, "/extdata/param/", param, ".Rdata"))
        index_allele_add <- which(alleles == allele_add)
        dt_sr_mk <- dt_sr[dt_sr$Marker == marker, ]
        marker_alcor <- rep(marker, n_al)
        model_bsr_alcor <- rep(dt_sr_mk$Model_BSR[1], n_al)
        model_bsr_alcor[index_allele_add] <- "Manual"
        model_fsr_alcor <- rep(dt_sr_mk$Model_FSR[1], n_al)
        model_fsr_alcor[index_allele_add] <- "Manual"
        model_dsr_alcor <- rep(dt_sr_mk$Model_DSR[1], n_al)
        model_dsr_alcor[index_allele_add] <- "Manual"
        model_m2sr_alcor <- rep(dt_sr_mk$Model_M2SR[1], n_al)
        model_m2sr_alcor[index_allele_add] <- "Manual"

        if(index_allele_add == 1){
          corrected_bsr_alcor <- c(NA, dt_sr_mk$Corrected_BSR[2:n_al])
          corrected_fsr_alcor <- c(NA, dt_sr_mk$Corrected_FSR[2:n_al])
          corrected_dsr_alcor <- c(NA, dt_sr_mk$Corrected_DSR[2:n_al])
        }else if(index_allele_add == n_al){
          corrected_bsr_alcor <- c(dt_sr_mk$Corrected_BSR[1:(n_al - 1)], NA)
          corrected_fsr_alcor <- c(dt_sr_mk$Corrected_FSR[1:(n_al - 1)], NA)
          corrected_dsr_alcor <- c(dt_sr_mk$Corrected_DSR[1:(n_al - 1)], NA)
        }else{
          corrected_bsr_alcor <- c(dt_sr_mk$Corrected_BSR[1:(index_allele_add - 1)], NA, dt_sr_mk$Corrected_BSR[(index_allele_add + 1):n_al])
          corrected_fsr_alcor <- c(dt_sr_mk$Corrected_FSR[1:(index_allele_add - 1)], NA, dt_sr_mk$Corrected_FSR[(index_allele_add + 1):n_al])
          corrected_dsr_alcor <- c(dt_sr_mk$Corrected_DSR[1:(index_allele_add - 1)], NA, dt_sr_mk$Corrected_DSR[(index_allele_add + 1):n_al])
        }

        dt_sr_mk_add <- data.table(Marker = marker_alcor,
                                   Allele = alleles_char,
                                   Model_BSR = model_bsr_alcor,
                                   Corrected_BSR = corrected_bsr_alcor,
                                   Average_BSR = info_bsr_add,
                                   Model_FSR = model_fsr_alcor,
                                   Corrected_FSR = corrected_fsr_alcor,
                                   Average_FSR = info_fsr_add,
                                   Model_DSR = model_dsr_alcor,
                                   Corrected_DSR = corrected_dsr_alcor,
                                   Average_DSR = info_dsr_add,
                                   Model_M2SR = model_m2sr_alcor,
                                   Average_M2SR = info_m2sr_add)

        index_target <- which(dt_sr$Marker == marker)
        if(index_target[1] == 1){
          dt_sr_other_mk <- dt_sr[dt_sr$Marker != marker, ]
          dt_sr <- rbind(dt_sr_mk_add, dt_sr_other_mk)
        }else if(index_target[length(index_target)] == nrow(dt_sr)){
          dt_sr_other_mk <- dt_sr[dt_sr$Marker != marker, ]
          dt_sr <- rbind(dt_sr_other_mk, dt_sr_mk_add)
        }else{
          dt_sr_other_mk_1 <- dt_sr[1:(index_target[1] - 1), ]
          dt_sr_other_mk_2 <- dt_sr[(index_target[length(index_target)] + 1):nrow(dt_sr), ]
          dt_sr <- rbind(dt_sr_other_mk_1, dt_sr_mk_add, dt_sr_other_mk_2)
        }
        save(list = c("dt_par", "dt_sr"), file = paste0(path_pack, "/extdata/param/", param, ".Rdata"))

        removeModal()
      }, ignoreInit = TRUE)

      #########################################
      # Delete the stutter ratio of an allele #
      #########################################

      observeEvent(input$act_sr_delete, {
        pos_select <- input$dt_info_sr_rows_selected

        if(isTruthy(pos_select)){
          showModal(modalDialog(
            title = "Delete the stutter ratio of an allele",
            "Click Ok button to delete the stutter ratio of the selected allele.",
            footer = tagList(
              actionButton(session$ns("act_delete_sr_ok"), "OK"),
              modalButton("Cancel")
            ),
            size = "l"
          ))
        }else{
          showModal(modalDialog(title = "Error", "Select a row!", easyClose = TRUE, footer = NULL))
        }
      })

      observeEvent(input$act_delete_sr_ok, {
        pos_select <- input$dt_info_sr_rows_selected
        dt_info_sr <- rv_sr$dt_info_sr
        param <- input$select_param_sr
        marker <- input$select_marker_sr

        # Update information on the stutter ratio
        rv_sr$dt_info_sr <- dt_info_sr[- pos_select, ]

        # Update dt_sr
        load(paste0(path_pack, "/extdata/param/", param, ".Rdata"))
        dt_sr_mk_delete <- copy(dt_sr)
        dt_sr_mk_delete <- dt_sr_mk_delete[dt_sr_mk_delete$Marker == marker, ]
        dt_sr_mk_delete <- dt_sr_mk_delete[- pos_select, ]

        index_target <- which(dt_sr$Marker == marker)
        if(index_target[1] == 1){
          dt_sr_other_mk <- dt_sr[dt_sr$Marker != marker, ]
          dt_sr <- rbind(dt_sr_mk_delete, dt_sr_other_mk)
        }else if(index_target[length(index_target)] == nrow(dt_sr)){
          dt_sr_other_mk <- dt_sr[dt_sr$Marker != marker, ]
          dt_sr <- rbind(dt_sr_other_mk, dt_sr_mk_delete)
        }else{
          dt_sr_other_mk_1 <- dt_sr[1:(index_target[1] - 1), ]
          dt_sr_other_mk_2 <- dt_sr[(index_target[length(index_target)] + 1):nrow(dt_sr), ]
          dt_sr <- rbind(dt_sr_other_mk_1, dt_sr_mk_delete, dt_sr_other_mk_2)
        }
        save(list = c("dt_par", "dt_sr"), file = paste0(path_pack, "/extdata/param/", param, ".Rdata"))

        removeModal()
      }, ignoreInit = TRUE)

      ############################################
      # Display information on the stutter ratio #
      ############################################

      output$dt_info_sr <- renderDataTable({
        dt_info_sr <- rv_sr$dt_info_sr
        if(!is.null(dt_info_sr)){
          datatable(
            dt_info_sr,
            colnames = c("Allele", "Average back stutter ratio", "Average forward stutter ratio", "Average double-back stutter ratio", "Average minus 2-nt stutter ratio"),
            selection = list(mode = "single", target = "row"),
            options = list(dom = "t", iDisplayLength = nrow(dt_info_sr), ordering = FALSE, autoWidth = TRUE),
            rownames = FALSE
          ) %>%
            formatRound(columns = c("Average_BSR", "Average_FSR", "Average_DSR"), digits = 4)
        }
      })
    }
  )
}
