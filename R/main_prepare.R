#' extract_auto_mk
#'
#' @description The function to extract autosomal markers from the data.table of the information on the kit
#' @param dt_kit The data.table of the information on the kit
extract_auto_mk <- function(dt_kit){
  setkey(dt_kit, Sex_Chr)
  dt_kit_auto <- dt_kit[""]
  return(dt_kit_auto[, Marker])
}

#' extract_peak_data
#'
#' @description The function to extract peak data
#' @param dt_csp The data.table of the crime stain profiles
#' @param markers_auto A vector of the autosomal markers
extract_peak_data <- function(dt_csp, markers_auto){
  peaks_list <- list()
  sizes_list <- list()
  heights_list <- list()

  n_mk <- length(markers_auto)

  names_csp <- names(dt_csp)
  index_allele <- grep("Allele", names_csp)
  index_size <- grep("Size", names_csp)
  index_height <- grep("Height", names_csp)

  sn_all_csp <- unique(dt_csp[, SampleName])
  setkey(dt_csp, SampleName)
  for(i in 1:length(sn_all_csp)){
    peaks_list_s1 <- list()
    sizes_list_s1 <- list()
    heights_list_s1 <- list()

    dt_csp_s1 <- dt_csp[sn_all_csp[i]]
    setkey(dt_csp_s1, Marker)
    for(j in 1:n_mk){
      dt_csp_s1_mk1 <- dt_csp_s1[markers_auto[j]]
      peaks_list_s1[[j]] <- as.numeric(dt_csp_s1_mk1[, index_allele, with = FALSE])
      sizes_list_s1[[j]] <- as.numeric(dt_csp_s1_mk1[, index_size, with = FALSE])
      heights_list_s1[[j]] <- as.numeric(dt_csp_s1_mk1[, index_height, with = FALSE])
    }

    names(peaks_list_s1) <- markers_auto
    names(sizes_list_s1) <- markers_auto
    names(heights_list_s1) <- markers_auto

    peaks_list[[i]] <- peaks_list_s1
    sizes_list[[i]] <- sizes_list_s1
    heights_list[[i]] <- heights_list_s1
  }

  names(peaks_list) <- sn_all_csp
  names(sizes_list) <- sn_all_csp
  names(heights_list) <- sn_all_csp
  return("peaks_list" = peaks_list,
         "sizes_list" = sizes_list,
         "heights_list" = heights_list)
}

#' make_gt_comb
#'
#' @description The function to make genotype combinations
#' @param hnc The hypothesized number of contributors
#' @param peaks A vector of peaks
make_gt_comb <- function(hnc, peaks){
  # Add drop-out
  peaks <- c(peaks, 99)

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
