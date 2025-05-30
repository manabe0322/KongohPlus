#' anameth_ui
#'
#' @description The function to create the ui module for the analysis method
anameth_ui <- function(id){
  ns <- NS(id)

  tabPanel("7. Set the analysis method",
           titlePanel("7. Set the analysis method"),
           br(),
           actionButton(ns("act_anameth_new"), "New"),
           actionButton(ns("act_anameth_edit"), "Edit"),
           actionButton(ns("act_anameth_copy"), "Copy"),
           actionButton(ns("act_anameth_delete"), "Delete"),
           br(),
           br(),
           dataTableOutput(ns("display_info_amameth"))
  )
}

#' create_dt_info_anameth
#'
#' @description The function to create the data.table for information on the analysis methods
#' @param path_pack The package path
create_dt_info_anameth <- function(path_pack){
  fn_anameth <- list.files(paste0(path_pack, "/extdata/anameth"))
  n_anameth <- length(fn_anameth)

  if(n_anameth != 0){
    info_kit <- info_stutter <- info_param <- info_dropin <- rep("", n_anameth)
    for(i in 1:n_anameth){
      dt_anameth_1 <- fread(paste0(path_pack, "/extdata/anameth/", fn_anameth[i]))
      info_kit[i] <- dt_anameth_1$Name[dt_anameth_1$Category == "Kit"]
      info_stutter[i] <- dt_anameth_1$Name[dt_anameth_1$Category == "Stutter data"]
      info_param[i] <- dt_anameth_1$Name[dt_anameth_1$Category == "Parameter data"]
      info_dropin[i] <- dt_anameth_1$Name[dt_anameth_1$Category == "Drop-in data"]
    }

    dt_info_anameth <- data.table("Analysis method" = gsub(".{4}$", "", fn_anameth),
                                  "Kit" = info_kit,
                                  "Stutter data" = info_stutter,
                                  "Parameter data" = info_param,
                                  "Drop-in data" = info_dropin)
    return(dt_info_anameth)
  }else{
    return(NULL)
  }
}

#' create_default_dt_anameth
#'
#' @description The function to create the default data.table for a new analysis method
#' @param path_pack The package path of KongohPlus
#' @param kit_name A kit name
#' @param stutter_name The name of the stutter data
#' @param param_name The name of the parameter data
#' @param dropin_name The name of the drop-in data
create_default_dt_anameth <- function(path_pack, kit_name, stutter_name, param_name, dropin_name){
  dt_kit <- fread(paste0(path_pack, "/extdata/kit/", kit_name, ".csv"))
  markers_auto <- unique(dt_kit[["Marker"]][dt_kit[["Chromosome"]] == "Autosomal"])
  n_mk <- length(markers_auto)

  load(paste0(path_pack, "/extdata/param/", param_name, ".RData"))
  load(paste0(path_pack, "/extdata/dropin/", dropin_name, ".RData"))

  default_dt_anameth <- data.table(Category = c("Kit", "Stutter data", "Parameter data", "Drop-in data",
                                                rep("AT", n_mk),
                                                rep("Genotype filter", 6),
                                                rep("LSAE", n_mk),
                                                rep("Peak height", 2),
                                                rep("Drop-in", 3)),
                                   Subcategory = c(rep("", 4),
                                                   markers_auto,
                                                   "Maximum BSR", "Maximum FSR", "Maximum DSR", "Maximum M2SR", "Maximum number of drop-in", "Maximum height of drop-in",
                                                   names(expect_ae),
                                                   "CV for allele", "CV for stutter",
                                                   "Frequency", "Shape parameter", "Rate parameter"),
                                   Name = c(kit_name, stutter_name, param_name, dropin_name, rep("", 2 * n_mk + 11)),
                                   Value = c(rep(NA, 4),
                                             rep(100, n_mk),
                                             0.3, 0.2, 0.05, 0.1, 2, 150,
                                             expect_ae,
                                             cv_allele, cv_stutter,
                                             freq_dropin, shape_di, rate_di
                                             ))
  return(default_dt_anameth)
}

#' anameth_server
#'
#' @description The function to create the server module for the analysis method
anameth_server <- function(id, path_pack){
  moduleServer(
    id,
    function(input, output, session){
      # Define the initial reactive values
      rv_anameth <- reactiveValues()
      rv_anameth$dt_info_anameth <- create_dt_info_anameth(path_pack)

      # Define a new analysis method
      observeEvent(input$act_anameth_new, {
        kit_all <- list.files(paste0(path_pack, "/extdata/kit"))
        kit_all <- gsub(".{4}$", "", kit_all)

        stutter_all <- list.files(paste0(path_pack, "/extdata/stutter"))
        stutter_all <- gsub(".{6}$", "", stutter_all)

        param_all <- list.files(paste0(path_pack, "/extdata/param"))
        param_all <- gsub(".{6}$", "", param_all)

        dropin_all <- list.files(paste0(path_pack, "/extdata/dropin"))
        dropin_all <- gsub(".{6}$", "", dropin_all)

        showModal(modalDialog(
          title = "Define a new analysis method",
          textInput(session$ns("anameth_name"), label = "Enter a name of the analysis method"),
          br(),
          selectInput(session$ns("kit_name"), label = "Select kit", choices = kit_all, selected = kit_all[1], selectize = FALSE),
          br(),
          selectInput(session$ns("stutter_name"), label = "Select stutter data", choices = stutter_all, selected = stutter_all[1], selectize = FALSE),
          br(),
          selectInput(session$ns("param_name"), label = "Select parameter data", choices = param_all, selected = param_all[1], selectize = FALSE),
          br(),
          selectInput(session$ns("dropin_name"), label = "Select drop-in data", choices = dropin_all, selected = dropin_all[1], selectize = FALSE),
          footer = tagList(
            actionButton(session$ns("act_add_anameth_new"), "Add the analysis method"),
            modalButton("Close")
          )
        ))
      })

      observeEvent(input$anameth_name, {
        anameth_name <- input$anameth_name
        anameth_all <- list.files(paste0(path_pack, "/extdata/anameth"))
        anameth_all <- gsub(".{4}$", "", anameth_all)

        if(is.element(anameth_name, anameth_all)){
          showFeedbackDanger(inputId = "anameth_name", text = "A file with the same name already exists.")
          disable("act_add_anameth_new")
        }else if(anameth_name == ""){
          hideFeedback("anameth_name")
          disable("act_add_anameth_new")
        }else{
          hideFeedback("anameth_name")
          enable("act_add_anameth_new")
        }
      })

      observeEvent(input$act_add_anameth_new, {
        anameth_name <- input$anameth_name
        kit_name <- input$kit_name
        stutter_name <- input$stutter_name
        param_name <- input$param_name
        dropin_name <- input$dropin_name

        if(all(!is.null(kit_name), !is.null(stutter_name), !is.null(param_name), !is.null(dropin_name))){
          load(paste0(path_pack, "/extdata/stutter/", stutter_name, ".RData"))
          load(paste0(path_pack, "/extdata/param/", param_name, ".RData"))
          load(paste0(path_pack, "/extdata/dropin/", dropin_name, ".RData"))

          if(kit_name == name_kit_stutter && kit_name == name_kit_param && kit_name == name_kit_dropin){
            dt_tmp <- data.table("Analysis method" = anameth_name,
                                 "Kit" = kit_name,
                                 "Stutter data" = stutter_name,
                                 "Parameter data" = param_name,
                                 "Drop-in data" = dropin_name)
            rv_anameth$dt_info_anameth <- rbind(rv_anameth$dt_info_anameth, dt_tmp)

            default_dt_anameth <- create_default_dt_anameth(path_pack, kit_name, stutter_name, param_name, dropin_name)
            write.csv(default_dt_anameth, paste0(path_pack, "/extdata/anameth/", anameth_name, ".csv"), row.names = FALSE)

          }
        }

        removeModal()
      })

      # Copy an analysis method
      observeEvent(input$act_anameth_copy, {
        index_select <- input$display_info_amameth_rows_selected
        if(isTruthy(index_select)){
          # Update information on the analysis method
          dt_info_anameth <- rv_anameth$dt_info_anameth
          dt_info_anameth_copy <- copy(dt_info_anameth)
          dt_info_anameth_copy <- dt_info_anameth_copy[index_select, ]
          anameth_select <- dt_info_anameth_copy[["Analysis method"]]
          dt_info_anameth_copy[["Analysis method"]] <- paste0(anameth_select, "_copy")
          dt_info_anameth <- rbind(dt_info_anameth, dt_info_anameth_copy)
          rv_anameth$dt_info_anameth <- dt_info_anameth

          # Make a copy of the selected analysis method
          dt_anameth_copy <- fread(paste0(path_pack, "/extdata/anameth/", anameth_select, ".csv"))
          write.csv(dt_anameth_copy, paste0(path_pack, "/extdata/anameth/", dt_info_anameth_copy[["Analysis method"]], ".csv"), row.names = FALSE)
        }else{
          showModal(modalDialog(title = "Error", "Select an analysis method!", easyClose = TRUE, footer = NULL))
        }
      })

      # Delete an analysis method
      observeEvent(input$act_anameth_delete, {
        index_select <- input$display_info_amameth_rows_selected
        if(isTruthy(index_select)){
          dt_info_anameth <- rv_anameth$dt_info_anameth
          file.remove(paste0(path_pack, "/extdata/anameth/", dt_info_anameth[["Analysis method"]][index_select], ".csv"))
          rv_anameth$dt_info_anameth <- dt_info_anameth[-index_select, ]
        }else{
          showModal(modalDialog(title = "Error", "Select an analysis method!", easyClose = TRUE, footer = NULL))
        }
      })

      # Display information on the analysis method
      output$display_info_amameth <- renderDataTable({
        dt_info_anameth <- rv_anameth$dt_info_anameth
        datatable(
          dt_info_anameth,
          selection = list(mode = "single", target = "row"),
          options = list(dom = "t", iDisplayLength = nrow(dt_info_anameth), ordering = FALSE, autoWidth = TRUE, scrollY = "650px", fixedHeader = TRUE),
          rownames = FALSE
        )
      })
    }
  )
}
