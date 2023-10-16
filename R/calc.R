#' change_data_type_csp
#'
#' @description The function to change data type of dt_csp_init
#' @param dt_csp_init The initial data.table of a crime stain profile
change_data_type_csp <- function(dt_csp_init){
  names(dt_csp_init) <- gsub(" ", "", names(dt_csp_init))

  name_dt_csp <- names(dt_csp_init)
  col_char <- c("SampleName", "Marker", name_dt_csp[grep("Allele", name_dt_csp)])
  dt_csp_init[, (col_char) := lapply(.SD, as.character), .SDcols = col_char]
  col_numeric <- c(name_dt_csp[grep("Size", name_dt_csp)], name_dt_csp[grep("Height", name_dt_csp)])
  dt_csp_init[, (col_numeric) := lapply(.SD, as.numeric), .SDcols = col_numeric]

  return(dt_csp_init)
}

#' make_gt_comb
#'
#' @description The function to make genotype combinations
#' @param hnc The hypothesized number of contributors
#' @param peaks
#' @param heights
make_gt_comb <- function(hnc, peaks, heights){

  n_peak <- length(peaks)

  # Genotypes of one contributor
  gts <- combinations(n = n_peak, r = 2, v = peaks, repeats.allowed = TRUE)

  n_gt <- nrow(gts)

  # Indices of genotypes of each contributor
  id_gt_comb <- permutations(n = n_gt, r = hnc, repeats.allowed = TRUE)

  n_gt_comb <- nrow(id_gt_comb)

  gt_comb <- matrix(0, n_gt_comb, 2 * hnc)

  for(i in 1:hnc){
    gt_comb[, c(2 * (i - 1) + 1, 2 * i)] <- gts[id_gt_comb[, i], ]
  }

  return(gt_comb)
}

#' make_eph_base
#'
#' @description The function to make the basis of the expected peak heights
#'
make_eph_base <- function(){

}

#' calc_tbar
#'
#' @description The function to calculate average peak heights per locus
calc_tbar <- function(){

}

#' calc_mmw
#'
#' @description The function to calculate the meam molecular weight
calc_mmw <- function(){

}

# Load csp
#dt_csp_init <- fread("C:/Users/manab/Documents/RStudio/KongohPlus/inst/extdata/examples/mix7_CSP.csv")
dt_csp_init <- fread("D:/RStudio_GitHub/KongohPlus/inst/extdata/examples/mix7_CSP.csv")

# Change data type of dt_csp_init
dt_csp_init <- change_data_type_csp(dt_csp_init)

# Remove sex chromosomal markers
dt_csp <- dt_csp_init[which(!is.element(dt_csp_init[, Marker], c("Yindel", "AMEL", "DYS391"))), , ]

