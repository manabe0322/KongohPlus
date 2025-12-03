#' check_pdf_env
#'
#' @description The function check PDF environment
check_pdf_env <- function(){
  message_pdf_env <- ""
  message_no_latex <- "LaTeX environment not found.<br>To enable PDF report generation, please install TinyTeX or another LaTeX distribution.<br><br>How to install TinyTeX:<br>&nbsp; &nbsp; install.packages('tinytex')<br>&nbsp; &nbsp; tinytex::install_tinytex()<br><br>"
  message_no_pandoc <- "Pandoc not found.<br>Please run this application from RStudio or install Pandoc manually.<br><br>How to install Pandoc (for Windows):<br>1. Visit https://pandoc.org/installing.html<br>2. Click 'Download the latest installer'<br>3. Download the Windows installer (*.msi)<br>4. Double-click the downloaded file and follow the installation wizard.<br><br>"

  has_tinytex_pkg <- requireNamespace("tinytex", quietly = TRUE)
  has_xelatex <- nzchar(Sys.which("xelatex"))

  if(has_tinytex_pkg){

    if(!tinytex::is_tinytex() && !has_xelatex){
      message_pdf_env <- paste0(message_pdf_env, message_no_latex)
    }
  }

  else{

    if(!has_xelatex){
      message_pdf_env <- paste0(message_pdf_env, message_no_latex)
    }
  }

  p <- find_pandoc()
  if(is.null(p$dir)){
    message_pdf_env <- paste0(message_pdf_env, message_no_pandoc)
  }

  return(message_pdf_env)
}
