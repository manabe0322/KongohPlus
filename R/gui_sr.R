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
                          uiOutput(ns("anameth_sr")),
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
                       dataTableOutput(ns("dt_sr"))
             )
           )
  )
}

#' sr_server
#'
#' @description The function to create the server module for the stutter ratio
sr_server <- function(id){
  moduleServer(
    id,
    function(input, output, session){
      output$anameth_sr <- renderUI({
        anameth_all <- list.files(paste0(path_pack, "/extdata/anameth"))
        anameth_all <- gsub(".{4}$", "", anameth_all)
        selectInput(session$ns("anameth_sr"), label = h4("Analysis method"), choices = anameth_all, selected = anameth_all[1])
      })
    }
  )
}
