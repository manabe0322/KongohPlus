# KongohPlus

## Intended Use

KongohPlus is primarily intended for research, development, and evaluation in forensic DNA mixture interpretation.

## Use in Forensic Casework

Before using KongohPlus in forensic casework, users are encouraged to evaluate and, where appropriate, validate its performance for their intended application and to ensure that analysts have received adequate training.

Users are responsible for determining whether KongohPlus is suitable for their intended application.

## Training

Training programs for users in Japan are planned to support the appropriate use of KongohPlus, including its underlying principles, validation, and interpretation of results.

Further details will be announced when available.

## Disclaimer

KongohPlus is provided without warranty, as described in the GNU General Public License v3.0. The software and its outputs should not be regarded as a substitute for professional judgment. Users are responsible for appropriately evaluating, implementing, and interpreting KongohPlus for their intended applications.

## Initial setup

1. Ensure that R (>= 4.5.2) is installed. It is available from the [R Development Core Team website](http://www.R-project.org).

2. Begin an R session.

3. Execute the following command in R to install required packages.

``` r
install.packages('https://github.com/manabe0322/KongohPlus/releases/download/v1.0.0/KongohPlus_v1.0.0/.zip', repos = NULL, type = 'win.binary')
install.packages(c("shiny", "cmaes", "combinat", "data.table", "dplyr", "DT", "ggplot2", "gridExtra", "gtools", "htmltools", "knitr", "magrittr", "parallel", "Rcpp", "rhandsontable", "rmarkdown", "shinyFeedback", "shinyFiles", "shinyjs", "shinythemes", "stringr", "waiter"))
```

4. Execute the following commands in R to launch the software interface.

``` r
library(KongohPlus)
KongohPlus()
```
## Settings Required to Enable Report Output

Run this application from RStudio or install Pandoc manually. 

How to install Pandoc (for Windows):

1. Visit https://pandoc.org/installing.html

2. Click 'Download the latest installer'

3. Download the Windows installer (*.msi)

4. Double-click the downloaded file and follow the installation wizard.
