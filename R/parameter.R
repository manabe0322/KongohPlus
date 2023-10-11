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

#' make_cand_peak
#'
#' @description The function to make detectable peaks in the crime stain profile under a genotype combination
#' @param gt_comb A genotype combination
#' @param sr_consider Information on consideration of back, forward, double-back, and minus 2-nt stutters
make_cand_peak <- function(gt_comb, sr_consider){

  # Define the output object
  cand_peak <- gt_comb

  # Whether back stutters are considered or not
  if(sr_consider[1]){
    cand_peak <- c(cand_peak, round(gt_comb - 1, 1))
  }

  # Whether forward stutters are considered or not
  if(sr_consider[2]){
    cand_peak <- c(cand_peak, round(gt_comb + 1, 1))
  }

  # whether double-back stutters are considered or not
  if(sr_consider[3]){
    cand_peak <- c(cand_peak, round(gt_comb - 2, 1))
  }

  # Whether minus 2-nt stutters are considered or not
  if(sr_consider[4]){
    cand_peak <- c(cand_peak, round(gt_comb[gt_comb %% 1 < 0.15] - 0.8, 1), round(gt_comb[gt_comb %% 1 >= 0.15] - 0.2, 1))
  }

  # Return
  return(sort(unique(cand_peak)))
}

#' remove_same_bin
#'
#' @description The function to remove two or more peaks which are located in the same bin
#' @param peaks Vector of peaks
#' @param sizes Vector of sizes
#' @param heights Vector of heights
remove_same_bin <- function(peaks, sizes, heights){

  # Count the number of peaks in each
  count_peak <- table(peaks)

  # Indices of two or more peaks which are located in the same bin
  index_multi <- which(count_peak != 1)

  # Remove two or more peaks which are located in the same bin
  if(length(index_multi) > 0){
    index_cut <- which(is.element(peaks, names(count_peak)[index_multi]))
    peaks <- peaks[-index_cut]
    sizes <- sizes[-index_cut]
    heights <- heights[-index_cut]
  }

  # Return
  return(list(peaks, sizes, heights))
}

#' remove_pu
#'
#' @description The function to remove pull-up peaks
#' @param peaks Vector of peaks
#' @param sizes Vector of sizes
#' @param heights Vector of heights
#' @param gt Genotype of one locus
#' @param dye Dye of one locus
#' @param al_size Allele sizes generated by the function "pick_al_size"
#' @param al_dye Allele dyes generated by the function "pick_al_size"
#' @param max_dist_pu The maximum distance between a pull-up peak and the parent allele
remove_pu <- function(peaks, sizes, heights, gt, dye, al_size, al_dye, max_dist_pu = 1){

  # Indices of alleles
  index_al <- which(is.element(peaks, gt))

  # Extract the sizes of the parent alleles
  index_parent <- which(al_dye != dye)
  sizes_pu <- al_size[index_parent]

  # Indices of pull-up peaks
  index_pu <- numeric(0)

  for(i in 1:length(peaks)){

    # Peaks other than alleles
    if(!is.element(i, index_al)){

      # Extract a size
      sz <- sizes[i]

      # Update indices of pull-up peaks
      if(length(which(abs(sz - sizes_pu) <= max_dist_pu)) != 0){
        index_pu <- c(index_pu, i)
      }
    }
  }

  # Remove pull-up peaks
  if(length(index_pu) > 0){
    peaks <- peaks[-index_pu]
    sizes <- sizes[-index_pu]
    heights <- heights[-index_pu]
  }

  # Return
  return(list(peaks, sizes, heights))
}

#' arrange_exp_data
#'
#' @description The function to arrange experimental data
#' @param dt_exp_init The initial data.table of experimental data
#' @param dt_gt_init The initial data.table of genotype data
#' @param par_markers Markers for parameters
#' @param at_all_mk Analytical thresholds of all markers
#' @param dye_all_mk Dyes of all markers
#' @param sr_consider_all_mk Information on consideration of back, forward, double-back, and minus 2-nt stutters in all markers
#' @param max_dist_pu The maximum distance between a pull-up peak and the parent allele
arrange_exp_data <- function(dt_exp_init, dt_gt_init, par_markers, at_all_mk, dye_all_mk, sr_consider_all_mk, max_dist_pu = 1){

  # Define output objects
  exp_list <- gt_list <- list()

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

  # Number of markers for parameters
  n_par_mk

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

    # Define objects for one sample
    exp_list_one <- list()
    gt_one <- matrix(NA, nrow = n_par_mk, ncol = 2)
    rownames(gt_one) <- par_markers
    colnames(gt_one) <- c("Allele1", "Allele2")

    for(j in 1:n_par_mk){

      # Extract data of one marker
      mk <- par_markers[j]
      index_mk <- which(markers == mk)
      peaks <- as.character(peak_all_mk[index_mk, ])
      sizes <- as.numeric(size_all_mk[index_mk, ])
      heights <- as.numeric(height_all_mk[index_mk, ])

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

      # Change data type of peaks (from character to numeric)
      peaks <- as.numeric(peaks)

      # Extract the genotype of one marker
      gt <- as.character(dt_gt_one[dt_gt_one[, Marker] == mk, index_gt_allele, with = FALSE])

      # Change data type of the genotype (from character to numeric)
      gt <- as.numeric(gt)

      # Assign the genotype
      gt_one[j, ] <- gt

      # Make detectable peaks in the profile under the genotype
      cand_peak <- make_cand_peak(gt, c(sr_consider_all_mk[rownames(sr_consider_all_mk) == mk, ]))

      # Indices of used peaks
      index_use <- which(is.element(peaks, cand_peak))

      if(length(index_use) > 0){

        # Extract used data
        peaks_use <- peaks[index_use]
        sizes_use <- sizes[index_use]
        heights_use <- heights[index_use]

        # Remove two or more peaks which are located in the same bin
        tmp <- remove_same_bin(peaks_use, sizes_use, heights_use)
        peaks_use <- tmp[[1]]
        sizes_use <- tmp[[2]]
        heights_use <- tmp[[3]]

        # Remove pull-up peaks
        tmp <- remove_pu(peaks, sizes, heights, gt, dye_all_mk[names(dye_all_mk) == mk], al_size, al_dye, max_dist_pu)
        peaks_use <- tmp[[1]]
        sizes_use <- tmp[[2]]
        heights_use <- tmp[[3]]

        # Assign
        info_peak <- rbind(peaks_use, sizes_use, heights_use)
        rownames(info_peak) <- c("Peak", "Size", "Height")
        exp_list_one[[j]] <- info_peak
      }else{
        exp_list_one[[j]] <- NULL
      }
    }

    # Assign
    names(exp_list_one) <- mk
    exp_list[[i]] <- exp_list_one
    gt_list[[i]] <- gt_one
  }

  # Return
  return(list(exp_list, gt_list))
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

# Consideration of stutters
sr_consider_all_mk <- matrix(TRUE, length(par_markers), 4)
rownames(sr_consider_all_mk) <- par_markers
colnames(sr_consider_all_mk) <- c("BS", "FS", "DS", "M2S")
sr_consider_all_mk[, 4] <- FALSE
sr_consider_all_mk[17, 4] <- TRUE
sr_consider_all_mk[19, 4] <- TRUE

############
# Analysis #
############

tmp <- arrange_exp_data(dt_exp_init, dt_gt_init, par_markers, at_all_mk, dye_all_mk, sr_consider_all_mk, max_dist_pu = 1)
exp_list <- tmp[[1]]
gt_list <- tmp[[2]]

