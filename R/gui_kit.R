#' kit_ui
#'
#' @description The function to create the ui module for the kit
kit_ui <- function(id){
  ns <- NS(id)

  tabPanel("Kit",
           titlePanel("Kit"),
           br(),
           sidebarLayout(
             sidebarPanel(width = 3,
                          actionButton(ns("act_kit_new"), label = "New"),
                          br(),
                          br(),
                          actionButton(ns("act_kit_edit"), label = "Edit"),
                          br(),
                          br(),
                          actionButton(ns("act_kit_delete"), label = "Delete")
             ),
             mainPanel(width = 9,
                       dataTableOutput(ns("dt_info_kit"))
             )
           )
  )
}

#' create_dt_info_kit
#'
#' @description The function to create the data.table for information on the kit
#' @param path_pack The package path
create_dt_info_kit <- function(path_pack){
  fn_kit <- list.files(paste0(path_pack, "/extdata/kit"))
  n_kit <- length(fn_kit)

  if(n_kit != 0){
    info_auto_mk <- info_sex_mk <- rep("", n_kit)
    for(i in 1:n_kit){
      dt_kit_1 <- fread(paste0(path_pack, "/extdata/kit/", fn_kit[i]))
      dt_kit_1 <- change_data_type_kit(dt_kit_1)
      info_auto_mk[i] <- length(unique(dt_kit_1$Marker[dt_kit_1$Chr == "Auto"]))
      info_sex_mk[i] <- length(unique(dt_kit_1$Marker[dt_kit_1$Chr == "Sex"]))
    }

    dt_info_kit <- data.table(Kit = gsub(".{4}$", "", fn_kit),
                              N_Auto_Mk = info_auto_mk,
                              N_Sex_Mk = info_sex_mk)
    return(dt_info_kit)
  }else{
    return(NULL)
  }
}

#' kit_server
#'
#' @description The function to create the server module for the kit
kit_server <- function(id, path_pack){
  moduleServer(
    id,
    function(input, output, session){

      ######################################
      # Define the initial reactive values #
      ######################################

      rv_kit <- reactiveValues()
      rv_kit$dt_info_kit <- create_dt_info_kit(path_pack)



      # New
      ## side-panel
        # Edit marker -> set marker name, and Chr (Auto or sex)
        # Add marker  -> set marker name, and Chr (Auto or sex)
        # Delete marker
      ## main-panel
        # data.table
          # Marker
          # Chr (Auto or sex)


      ##############################################
      # Display information on the analysis method #
      ##############################################

      output$dt_info_kit <- renderDataTable({
        dt_info_kit <- rv_kit$dt_info_kit
        datatable(
          dt_info_kit,
          colnames = c("Kit", "Number of autosomal markers", "Number of sex-chromosomal markers"),
          selection = list(mode = "single", target = "row"),
          options = list(dom = "t", iDisplayLength = nrow(dt_info_kit), ordering = FALSE, autoWidth = TRUE),
          rownames = FALSE
        )
      })
    }
  )
}
