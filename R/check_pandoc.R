#' check_pandoc
#'
#' @description The function check pandoc environment
check_pandoc <- function(){
  message_pandoc <- ""

  p <- find_pandoc()

  if(is.null(p$dir)){
    message_pandoc <- "Pandoc not found.<br>Please run this application from RStudio or install Pandoc manually.<br><br>How to install Pandoc (for Windows):<br>1. Visit https://pandoc.org/installing.html<br>2. Click 'Download the latest installer'<br>3. Download the Windows installer (*.msi)<br>4. Double-click the downloaded file and follow the installation wizard.<br><br>"
  }

  return(message_pandoc)
}
