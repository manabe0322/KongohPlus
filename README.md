# KongohPlus

## Initial setup

1. Ensure that R (>= 4.4.3) is installed. It is available from the [R Development Core Team website](http://www.R-project.org).

2. Begin an R session.

3. Execute the following command in R to install required packages.

``` r
install.packages('https://github.com/manabe0322/KongohPlus/releases/download/v0.21.4/KongohPlus_0.21.7.zip', repos = NULL, type = 'win.binary')
install.packages(c("shiny", "combinat", "data.table", "dplyr", "DT", "GenSA", "ggplot2", "gridExtra", "gtools", "htmltools", "kableExtra", "knitr", "magrittr", "parallel", "Rcpp", "rhandsontable", "rmarkdown", "shinyFeedback", "shinyFiles", "shinyjs", "shinythemes", "stringr", "waiter"))
```

4. Execute the following commands in R to launch the software interface.

``` r
library(KongohPlus)
KongohPlus()
```
