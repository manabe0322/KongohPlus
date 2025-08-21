#' file_ui
#'
#' @description The function to create the ui module for loading files
file_ui <- function(id){
  ns <- NS(id)

  tabPanel("File",
           sidebarLayout(
             sidebarPanel(width = 3,
                          fileInput(ns("file_csp"), label = h4("Crime stain profile"), accept = ".csv", buttonLabel = "Load"),
                          fileInput(ns("file_ref"), label = h4("Reference profile"), accept = ".csv", buttonLabel = "Load"),
                          fileInput(ns("file_af"), label = h4("Allele frequency"), accept = ".csv", buttonLabel = "Load"),
                          br(),
                          actionButton(ns("act_deconvo"), label = "Deconvolution", class = "btn btn-primary btn-lg")
             ),
             mainPanel(width = 9,
                       tabsetPanel(id = ns("tabset_load_file"),
                                   type = "tabs",
                                   tabPanel("Crime stain profile",
                                            dataTableOutput(ns("dt_csp"))
                                   ),
                                   tabPanel("Reference profile",
                                            dataTableOutput(ns("dt_ref"))
                                   ),
                                   tabPanel("Allele frequency",
                                            dataTableOutput(ns("dt_af"))
                                   )
                       )
             )
           )
  )
}

#' file_server
#'
#' @description The function to create the server module for loading files
file_server <- function(id, session_top, bool_load_proj){
  moduleServer(
    id,
    function(input, output, session){
      rv_input <- reactiveValues()
      rv_input$dt_csp <- NULL
      rv_input$dt_ref <- NULL
      rv_input$dt_af <- NULL
      rv_input$fn_csp <- NULL
      rv_input$fn_ref <- NULL
      rv_input$fn_af <- NULL

      # --- Load crime stain profiles --- #

      observeEvent(input$file_csp, {
        file_input <- input$file_csp
        if(!is.null(file_input)){
          dt_csp <- fread(file_input$datapath)
          dt_csp <- change_data_type_csp(dt_csp)
          rv_input$dt_csp <- dt_csp
          rv_input$fn_csp <- file_input$name
          updateTabsetPanel(inputId = "tabset_load_file", selected = "Crime stain profile")
        }
      })

      # --- Load reference profiles --- #

      observeEvent(input$file_ref, {
        file_input <- input$file_ref
        if(!is.null(file_input)){
          dt_ref <- fread(file_input$datapath)
          dt_ref <- change_data_type_ref(dt_ref)
          rv_input$dt_ref <- dt_ref
          rv_input$fn_ref <- file_input$name
          updateTabsetPanel(inputId = "tabset_load_file", selected = "Reference profile")
        }
      })

      # --- Load allele frequencies --- #

      observeEvent(input$file_af, {
        file_input <- input$file_af
        if(!is.null(file_input)){
          dt_af <- fread(file_input$datapath)
          dt_af <- change_data_type_af(dt_af)
          rv_input$dt_af <- dt_af
          rv_input$fn_af <- file_input$name
          updateTabsetPanel(inputId = "tabset_load_file", selected = "Allele frequency")
        }
      })

      # --- Display crime stain profiles --- #

      output$dt_csp <- renderDataTable({
        dt_csp <- rv_input$dt_csp

        if(!is.null(dt_csp)){
          name_dt_csp <- names(dt_csp)

          mat_alleles <- as.matrix(dt_csp[, grep("Allele", name_dt_csp), with = FALSE])
          mat_heights <- as.matrix(dt_csp[, grep("Height", name_dt_csp), with = FALSE])
          n_row <- nrow(mat_alleles)
          paste_alleles <- rep(NA, n_row)
          paste_heights <- rep(NA, n_row)
          for(i in 1:n_row){
            alleles_mk <- mat_alleles[i, ]
            alleles_mk <- alleles_mk[!is.na(alleles_mk)]
            if(length(alleles_mk) != 0){
              paste_alleles[i] <- paste0(alleles_mk, collapse = ", ")
            }

            heights_mk <- mat_heights[i, ]
            heights_mk <- heights_mk[!is.na(heights_mk)]
            if(length(heights_mk) != 0){
              paste_heights[i] <- paste0(heights_mk, collapse = ", ")
            }
          }

          dt_csp_display <- data.table(SampleName = factor(dt_csp$SampleName, levels = unique(dt_csp$SampleName), labels = unique(dt_csp$SampleName)),
                                       Marker = factor(dt_csp$Marker, levels = unique(dt_csp$Marker), labels = unique(dt_csp$Marker)),
                                       Allele = paste_alleles,
                                       Height = paste_heights)


          index_sn <- which(name_dt_csp == "SampleName") - 1
          index_mk <- which(name_dt_csp == "Marker") - 1
          index_num <- setdiff(0:(ncol(dt_csp) - 1), c(index_sn, index_mk))

          datatable(
            dt_csp_display,
            filter = "top",
            selection = "none",
            extensions = c("FixedColumns", "FixedHeader"),
            options = list(dom = "t", iDisplayLength = nrow(dt_csp_display), ordering = FALSE, autoWidth = FALSE, scrollX = TRUE, scrollY = "calc(100vh - 200px)",
                           fixedHeader = TRUE, fixedColumns = list(leftColumns = 1, rightColumns = 0),
                           columnDefs = list(list(targets = c(2, 3), searchable = FALSE))),
            rownames = FALSE
          ) %>%
            formatStyle(names(dt_csp_display), lineHeight = "80%")
        }
      })

      # --- Display reference profiles --- #

      output$dt_ref <- renderDataTable({
        dt_ref <- rv_input$dt_ref

        if(!is.null(dt_ref)){
          dt_ref_display <- NULL
          samples <- unique(dt_ref$SampleName)
          n_sample <- length(samples)
          for(i in 1:n_sample){
            dt_tmp <- dt_ref[dt_ref$SampleName == samples[i]]

            gts <- rep(NA, nrow(dt_tmp))
            for(j in 1:nrow(dt_tmp)){
              gt1 <- c(dt_tmp$Allele1[j], dt_tmp$Allele2[j])
              gt1 <- gt1[!is.na(gt1)]
              gts[j] <- paste0(gt1, collapse = ", ")
            }

            dt_tmp2 <- data.table(Marker = dt_tmp$Marker,
                                  Gt = gts)
            names(dt_tmp2) <- c("Marker", samples[i])

            if(i == 1){
              dt_ref_display <- dt_tmp2
            }else{
              dt_ref_display <- full_join(dt_ref_display, dt_tmp2, by = c("Marker"))
            }
          }

          dt_ref_display$Marker <- factor(dt_ref_display$Marker, levels = unique(dt_ref_display$Marker), labels = unique(dt_ref_display$Marker))

          datatable(
            dt_ref_display,
            filter = "top",
            selection = "none",
            extensions = c("FixedColumns", "FixedHeader"),
            options = list(dom = "t", iDisplayLength = nrow(dt_ref_display), ordering = FALSE, autoWidth = FALSE, scrollX = TRUE, scrollY = "calc(100vh - 200px)",
                           fixedHeader = TRUE, fixedColumns = list(leftColumns = 1, rightColumns = 0),
                           columnDefs = list(list(targets = 1:(ncol(dt_ref_display) - 1), searchable = FALSE))),
            rownames = FALSE
          ) %>%
            formatStyle(names(dt_ref_display), lineHeight = "80%")
        }
      })

      # --- Display allele frequencies --- #

      output$dt_af <- renderDataTable({
        dt_af <- rv_input$dt_af
        if(!is.null(dt_af)){

          datatable(
            dt_af,
            selection = "none",
            extensions = c("FixedColumns", "FixedHeader"),
            options = list(dom = "t", iDisplayLength = nrow(dt_af), ordering = FALSE, autoWidth = FALSE, scrollX = TRUE, scrollY = "calc(100vh - 150px)",
                           fixedHeader = TRUE, fixedColumns = list(leftColumns = 1, rightColumns = 0)),
            rownames = FALSE
            ) %>%
              formatStyle(names(dt_af), lineHeight = "80%")
        }
      })

      observeEvent(input$act_deconvo, {
        updateNavbarPage(session = session_top, inputId = "navbar", selected = "Deconvolution")
      })

      return(rv_input)
    }
  )
}
