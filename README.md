# KongohPlus

## Initial setup

1. Ensure that R (>= 4.5.2) is installed. It is available from the [R Development Core Team website](http://www.R-project.org).

2. Begin an R session.

3. Execute the following command in R to install required packages.

``` r
install.packages('https://github.com/manabe0322/KongohPlus/releases/download/v0.24.0/KongohPlus_v0.24.0.zip', repos = NULL, type = 'win.binary')
install.packages(c("shiny", "cmaes", "combinat", "data.table", "dplyr", "DT", "ggplot2", "gridExtra", "gtools", "htmltools", "kableExtra", "knitr", "magrittr", "parallel", "Rcpp", "rhandsontable", "rmarkdown", "shinyFeedback", "shinyFiles", "shinyjs", "shinythemes", "stringr", "waiter"))
```

4. Execute the following commands in R to launch the software interface.

``` r
library(KongohPlus)
KongohPlus()
```
## Settings Required to Enable Report Output

Install TinyTeX or another LaTeX distribution to enable PDF report generation. Execute the following command in R to install  TinyTeX.

``` r
install.packages('tinytex')
tinytex::install_tinytex()
```

In addition, run this application from RStudio or install Pandoc manually. 

How to install Pandoc (for Windows):

1. Visit https://pandoc.org/installing.html

2. Click 'Download the latest installer'

3. Download the Windows installer (*.msi)

4. Double-click the downloaded file and follow the installation wizard.
