#' sr_ui
#'
#' @description The function to create the ui module for the stutter ratio
sr_ui <- function(id){
  ns <- NS(id)

  # Analysis method毎に以下の機能をつける
  # ローカス毎にデータテーブル
  # 値を修正するモーダル
  # 値を追加するモーダル
  # 値を消す機能
  # alcorにあるデータにおける各SRのmodel名は、マニュアルで追加・修正したものについては"Manual"等に変更する
  # 今後: ファイルを取り込み、一括で上書きする機能をつける
  # 今後: 内挿・外挿により推定する機能をつける

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
    marker <- dt_alcor_combined$Marker[1]
    dt_info_sr <- create_dt_info_sr(dt_alcor_combined, marker)
  }else{
    dt_info_sr <- NULL
  }
  return(dt_info_sr)
}

#' create_dt_info_sr
#'
#' @description The function to create the data.table for information on the stutter ratio
create_dt_info_sr <- function(dt_alcor_combined, marker){
  index_mk <- which(dt_alcor_combined$Marker == marker)
  dt_info_sr <- data.table(Allele = dt_alcor_combined$Allele[index_mk],
                           Average_BSR = dt_alcor_combined$Average_BSR[index_mk],
                           Average_FSR = dt_alcor_combined$Average_FSR[index_mk],
                           Average_DSR = dt_alcor_combined$Average_DSR[index_mk])
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



      output$select_param_sr <- renderUI({
        param_all <- list.files(paste0(path_pack, "/extdata/param"))
        param_all <- gsub(".{6}$", "", param_all)
        selectInput(session$ns("select_param_sr"), label = h4("Parameter data"), choices = param_all, selected = param_all[1], selectize = FALSE)
      })

      observeEvent(input$select_param_sr, {
        load(paste0(path_pack, "/extdata/param/", input$select_param_sr, ".Rdata"))
        markers_sr <- unique(dt_alcor_combined$Marker)

        output$select_marker_sr <- renderUI({
          selectInput(session$ns("select_marker_sr"), label = h4("Marker"), choices = markers_sr, selected = markers_sr[1], selectize = FALSE)
        })
      })

      observeEvent(input$select_marker_sr, {
        load(paste0(path_pack, "/extdata/param/", input$select_param_sr, ".Rdata"))
        rv_sr$dt_info_sr <- create_dt_info_sr(dt_alcor_combined, input$select_marker_sr)
      })

      output$dt_info_sr <- renderDataTable({
        dt_info_sr <- rv_sr$dt_info_sr
        if(!is.null(dt_info_sr)){
          datatable(
            dt_info_sr,
            colnames = c("Allele", "Average back stutter ratio", "Average forward stutter ratio", "Average double-back stutter ratio"),
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
