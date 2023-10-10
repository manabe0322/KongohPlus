#' change_data_type_exp
#'
#' @description The function to change data type of dt_exp_init
#' @param dt_exp_init The initial data.table of experimental data
change_data_type_exp <- function(dt_exp_init){

  # Remove spaces of the names
  names(dt_exp_init) <- gsub(" ", "", names(dt_exp_init))

  name_dt_exp <- names(dt_exp_init)
  col_char <- c("SampleName", "Marker", name_dt_exp[grep("Allele", name_dt_exp)])
  dt_exp_init[, (col_char) := lapply(.SD, as.character), .SDcols = col_char]
  col_numeric <- c(name_dt_exp[grep("Size", name_dt_exp)], name_dt_exp[grep("Height", name_dt_exp)])
  dt_exp_init[, (col_numeric) := lapply(.SD, as.numeric), .SDcols = col_numeric]

  # Replace "" to NA
  dt_exp_init[which(dt_exp_init == "", arr.ind = TRUE)] <- NA

  # Return
  return(dt_exp_init)
}

#' change_data_type_gt
#'
#' @description The function to change data type of dt_gt_init
#' @param dt_gt_init The initial data.table of genotype data
change_data_type_gt <- function(dt_gt_init){

  # Remove spaces of the names
  names(dt_gt_init) <- gsub(" ", "", names(dt_gt_init))

  name_dt_gt <- names(dt_gt_init)
  col_char <- c("SampleName", "Marker", name_dt_gt[grep("Allele", name_dt_gt)])
  dt_gt_init[, (col_char) := lapply(.SD, as.character), .SDcols = col_char]

  # Replace "" to NA
  dt_gt_init[which(dt_gt_init == "", arr.ind = TRUE)] <- NA

  # Return
  return(dt_gt_init)
}

#' change_data_type_seq
#'
#' @description The function to change data type of dt_seq_init
#' @param dt_seq_init The initial data.table of sequence data
change_data_type_seq <- function(dt_seq_init){

  # Remove spaces of the names
  names(dt_seq_init) <- gsub(" ", "", names(dt_seq_init))

  name_dt_seq <- names(dt_seq_init)
  col_char <- c("Marker", "Allele", "RepeatRegion", "FlankingRegion")
  dt_seq_init[, (col_char) := lapply(.SD, as.character), .SDcols = col_char]
  col_numeric <- c("Count")
  dt_seq_init[, (col_numeric) := lapply(.SD, as.numeric), .SDcols = col_numeric]

  # Replace "" to NA
  dt_seq_init[which(dt_seq_init == "", arr.ind = TRUE)] <- NA

  # Return
  return(dt_seq_init)
}

#' pick_al_size
#'
#' @description The function to pick up allele sizes
#' @param dt_exp_one The data.table of the experimental data of one sample
#' @param dt_gt_one The data.table of the genotype data of one sample
#' @param at_all_mk Analytical thresholds of all markers
#' @param dye_all_mk Dyes of all markers
pick_al_size <- function(dt_exp_one, dt_gt_one, at_all_mk, dye_all_mk){

  # Define objects
  al_mk <- al_peak <- al_size <- al_dye <- numeric(0)

  # Names of dt_exp_one
  name_dt_exp <- names(dt_exp_one)

  # Extract data of one sample
  markers <- dt_exp_one[, Marker]
  peak_all_mk <- dt_exp_one[, grep("Allele", name_dt_exp), with = FALSE]
  size_all_mk <- dt_exp_one[, grep("Size", name_dt_exp), with = FALSE]
  height_all_mk <- dt_exp_one[, grep("Height", name_dt_exp), with = FALSE]

  # Number of markers
  n_mk <- length(markers)

  # Column indices of dt_gt_one
  index_gt_allele <- grep("Allele", names(dt_gt_one))

  for(i in 1:n_mk){

    # Extract data of one marker
    mk <- markers[i]
    peaks <- as.character(peak_all_mk[i, ])
    sizes <- as.numeric(size_all_mk[i, ])
    heights <- as.numeric(height_all_mk[i, ])

    # Remove NA
    index_not_na <- which(!is.na(peaks))
    peaks <- peaks[index_not_na]
    sizes <- sizes[index_not_na]
    heights <- heights[index_not_na]

    # Remove peaks below AT
    index_above_at <- which(heights >= at_all_mk[names(at_all_mk) == mk])
    peaks <- peaks[index_above_at]
    sizes <- sizes[index_above_at]
    heights <- heights[index_above_at]

    # Extract the genotype of one marker
    gt <- as.character(dt_gt_one[dt_gt_one[, Marker] == mk, index_gt_allele, with = FALSE])

    # Extract information on the allele
    index_al <- which(is.element(peaks, gt))
    al_mk <- c(al_mk, rep(mk, length(index_al)))
    al_peak <- c(al_peak, peaks[index_al])
    al_size <- c(al_size, sizes[index_al])
    al_dye <- c(al_dye, rep(dye_all_mk[names(dye_all_mk) == mk], length(index_al)))
  }

  # Make the output object
  al_size_info <- rbind(al_peak, al_size, al_dye)
  colnames(al_size_info) <- al_mk

  # Return
  return(al_size_info)
}

#' arrange_exp_data
#'
#' @description The function to arrange experimental data
arrange_exp_data <- function(dt_exp_init, dt_gt_init, par_mk, sex_mk, at_all_mk, dye_all_mk, sr_consider){

  # Names of dt_exp_init
  name_dt_exp <- names(dt_exp_init)

  # Column indices of dt_exp_init
  index_exp_peak <- grep("Allele", name_dt_exp)
  index_exp_size <- grep("Size", name_dt_exp)
  index_exp_height <- grep("Height", name_dt_exp)

  # Extract sample names
  sample_names <- unique(dt_exp_init[, SampleName])

  # Names of dt_gt_init
  name_dt_gt <- names(dt_gt_init)

  # Column indices of dt_gt_init
  index_gt_allele <- grep("Allele", name_dt_gt)

  # Define output objects
  exp_list <- gt_list <- list()

  for(i in 1:length(sample_names)){

    # Extract data of one sample
    dt_exp_one <- dt_exp_init[dt_exp_init[, SampleName] == sample_names[i], ]
    markers <- dt_exp_one[, Marker]
    peak_all_mk <- dt_exp_one[, index_exp_peak, with = FALSE]
    size_all_mk <- dt_exp_one[, index_exp_size, with = FALSE]
    height_all_mk <- dt_exp_one[, index_exp_height, with = FALSE]
    dt_gt_one <- dt_gt_init[dt_gt_init[, SampleName] == sample_names[i], ]

    # Pick up allele sizes
    al_size_info <- pick_al_size(dt_exp_one, dt_gt_one, at_all_mk, dye_all_mk)
    al_size <- as.numeric(al_size_info[2, ])
    al_dye <- al_size_info[3, ]

  }

}


##############
# Load files #
##############

# Load experimental data
#dt_exp_init <- fread("C:/Users/manab/Documents/RStudio/KongohPlus/inst/extdata/examples/expData_kakeiken.csv")
dt_exp_init <- fread("D:/RStudio_GitHub/KongohPlus/inst/extdata/examples/expData_kakeiken.csv")

# Change data type of dt_exp_init
dt_exp_init <- change_data_type_exp(dt_exp_init)

# Load genotype data
#dt_gt_init <- fread("C:/Users/manab/Documents/RStudio/KongohPlus/inst/extdata/examples/gtData_kakeiken.csv")
dt_gt_init <- fread("D:/RStudio_GitHub/KongohPlus/inst/extdata/examples/gtData_kakeiken.csv")

# Change data type of dt_gt_init
dt_gt_init <- change_data_type_gt(dt_gt_init)

# Load sequence data
#dt_seq_init <- fread("C:/Users/manab/Documents/RStudio/KongohPlus/inst/extdata/examples/Sequence-data_GF_Example.csv")
dt_seq_init <- fread("D:/RStudio_GitHub/KongohPlus/inst/extdata/examples/Sequence-data_GF_Example.csv")

# Change data type of dt_seq_init
dt_seq_init <- change_data_type_seq(dt_seq_init)

# Check files
  # sample names
  # markers


############
# Settings #
############

# Markers
markers <- unique(dt_exp_init[, Marker])

# Sex chromosomal markers
sex_markers <- c("Yindel", "AMEL", "DYS391")

# Markers for parameters
par_markers <- setdiff(markers, sex_markers)

# Analytical thresholds
at_all_mk <- c(rep(12, 5), rep(20, 6), rep(13, 4), rep(16, 5), rep(19, 4))
names(at_all_mk) <- markers

# Dyes
dye_all_mk <- c(rep("B", 5), rep("G", 6), rep("Y", 4), rep("R", 5), rep("P", 4))
names(dye_all_mk) <- markers

# Number of markers
#n_mk <- length(par_mk)

# Consideration of stutters
sr_consider <- matrix(TRUE, nl, 4)
sr_consider[, 4] <- FALSE
sr_consider[17, 4] <- TRUE
sr_consider[19, 4] <- TRUE


