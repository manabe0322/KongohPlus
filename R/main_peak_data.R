#' update_heights_list_mk
#'
#' @description The function to update information on heights of a marker combining information on the peaks in each replicate
#' @param peaks_all_rep The vector of peaks of a marker which are observed across all replicates
#' @param peaks_list_mk The list of peaks per replicate in a marker
#' @param heights_list_mk The list of heights per replicate in a marker
update_heights_list_mk <- function(peaks_all_rep, peaks_list_mk, heights_list_mk){
  n_peaks <- length(peaks_all_rep)
  n_rep <- length(peaks_list_mk)
  height_list_mk_new <- vector(mode = "list", length = n_rep)
  for(i in 1:n_rep){
    heights_new <- rep(0, n_peaks)
    index_heights <- match(peaks_list_mk[[i]], peaks_all_rep)
    heights_new[index_heights] <- heights_list_mk[[i]]
    height_list_mk_new[[i]] <- heights_new
  }
  return(height_list_mk_new)
}

#' simple_linear_regression
#'
#' @description The function to perform simple linear regression
#' @param x The vector of data of x-axis
#' @param y The vector of data of y-axis
simple_linear_regression <- function(x, y){
  slr <- lm(y~x)
  param <- coef(slr)
  param <- param[c(2, 1)]
  names(param) <- c("slope", "intercept")
  return(param)
}

#' make_peak_data
#'
#' @description The function to make peak data
#' @param dt_csp The data.table of the crime stain profiles
#' @param dt_kit_auto The data.table of the kit information for autosomal STRs
#' @param markers_use A vector of the autosomal markers used
#' @param ats_use A vector of the analytical thresholds for each marker
make_peak_data <- function(dt_csp, dt_kit_auto, markers_use, ats_use){
  # Column indices for alleles, sizes, and heights
  names_csp <- names(dt_csp)
  index_allele <- grep("Allele", names_csp)
  index_height <- grep("Height", names_csp)

  name_csp_all <- unique(dt_csp[, SampleName])
  n_rep <- length(name_csp_all)
  n_mk <- length(markers_use)

  setindex(dt_kit_auto, Marker)
  setindex(dt_csp, Marker)

  peaks_list <- sizes_list <- heights_list <- vector(mode = "list", length = n_mk)

  for(i in 1:n_mk){
    mk <- markers_use[i]
    at <- ats_use[mk]
    dt_kit_auto_mk <- dt_kit_auto[Marker == mk]
    dt_csp_mk <- dt_csp[Marker == mk]
    setindex(dt_csp_mk, SampleName)
    peaks_list_mk <- heights_list_mk <- vector(mode = "list", length = n_rep)
    for(j in 1:n_rep){
      # Extract information on peaks of a marker in one replicate
      rep1 <- name_csp_all[j]
      dt_csp_mk_rep1 <- dt_csp_mk[SampleName == rep1]
      peaks <- as.numeric(dt_csp_mk_rep1[, index_allele, with = FALSE])
      peaks <- peaks[!is.na(peaks)]
      heights <- as.numeric(dt_csp_mk_rep1[, index_height, with = FALSE])
      heights <- heights[!is.na(heights)]

      # Keep peaks above AT
      if(length(peaks) != 0){
        index_above_at <- which(heights >= at)
        peaks <- peaks[index_above_at]
        heights <- heights[index_above_at]
      }

      peaks_list_mk[[j]] <- peaks
      heights_list_mk[[j]] <- heights
    }

    # Update information on heights of a marker combining information on the peaks in each replicate
    peaks_all_rep <- sort(unique(unlist(peaks_list_mk)))
    height_list_mk_new <- update_heights_list_mk(peaks_all_rep, peaks_list_mk, heights_list_mk)

    # Extract information on sizes
    alleles_kit <- as.numeric(dt_kit_auto_mk[["Allele"]])
    sizes_kit <- dt_kit_auto_mk[["Size"]]
    sizes_all_rep <- rep(0, length(peaks_all_rep))
    if(length(peaks_all_rep) != 0){
      index_sizes <- match(peaks_all_rep, alleles_kit)
      index_sizes_na <- which(is.na(index_sizes))
      index_sizes_exist <- which(!is.na(index_sizes))
      sizes_all_rep[index_sizes_exist] <- sizes_kit[index_sizes[index_sizes_exist]]

      # Interpolation or extrapolation to determine the sizes without information
      if(length(index_sizes_na) != 0){
        rep_uni_len <- dt_kit_auto_mk[["Repeat Unit Length"]][1]
        alleles_scale_change <- change_var_al_plot(alleles_kit, rep_uni_len)
        param_tmp <- simple_linear_regression(alleles_scale_change, sizes_kit)
        slope <- param_tmp["slope"]
        intercept <- param_tmp["intercept"]
        peaks_scale_change <- change_var_al_plot(peaks_all_rep[index_sizes_na], rep_uni_len)
        sizes_all_rep[index_sizes_na] <- slope * peaks_scale_change + intercept
      }
    }

    # Determine the size for the drop-out allele
    size_do <- median(sizes_kit)

    # Add information on the drop-out peak
    peaks_list[[i]] <- c(peaks_all_rep, 99)
    sizes_list[[i]] <- c(sizes_all_rep, size_do)
    heights_all_rep_list <- vector(mode = "list", length = n_rep)
    for(j in 1:n_rep){
      heights_all_rep_list[[j]] <- c(height_list_mk_new[[j]], 0)
    }
    heights_list[[i]] <- heights_all_rep_list
  }

  names(peaks_list) <- markers_use
  names(sizes_list) <- markers_use
  names(heights_list) <- markers_use
  return(list("peaks_list" = peaks_list,
              "sizes_list" = sizes_list,
              "heights_list" = heights_list))
}
