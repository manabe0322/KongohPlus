#' arrange_exp_data
#'
#' @description The function to arrange the experimental data
#' @param dt_exp The data.table of the experimental data
#' @param dt_gt The data.table of the genotype data
#' @param markers_auto Autosomal markers
#' @param sr_consider_auto Information on consideration of back, forward, double-back, and minus 2-nt stutters in autosomal markers
#' @param at_auto Analytical thresholds for autosomal markers
arrange_exp_data <- function(dt_exp, dt_gt, markers_auto, sr_consider_auto, at_auto){
  peak_list <- height_list <- gt_list <- list()

  name_dt_exp <- names(dt_exp)
  index_exp_peak <- grep("Allele", name_dt_exp)
  index_exp_height <- grep("Height", name_dt_exp)
  sample_names <- unique(dt_exp[, SampleName])

  name_dt_gt <- names(dt_gt)
  index_gt_allele <- grep("Allele", name_dt_gt)

  n_mk_auto <- length(markers_auto)

  setindex(dt_exp, SampleName)
  setindex(dt_gt, SampleName)

  for(i in 1:length(sample_names)){
    dt_exp_one <- dt_exp[SampleName == sample_names[i]]
    dt_gt_one <- dt_gt[SampleName == sample_names[i]]
    setindex(dt_exp_one, Marker)
    setindex(dt_gt_one, Marker)

    # Define objects for one sample
    peak_list_one <- height_list_one <- list()
    gt_one <- matrix(NA, nrow = n_mk_auto, ncol = 2)
    rownames(gt_one) <- markers_auto
    colnames(gt_one) <- c("Allele1", "Allele2")

    for(j in 1:n_mk_auto){
      # Extract data of one marker
      mk <- markers_auto[j]
      dt_exp_one_mk <- dt_exp_one[Marker == mk]
      peaks <- as.character(dt_exp_one_mk[, index_exp_peak, with = FALSE])
      heights <- as.numeric(dt_exp_one_mk[, index_exp_height, with = FALSE])

      # Remove NA
      index_not_na <- which(!is.na(peaks))
      if(length(index_not_na) == 0){
        peak_list_one[[j]] <- character(0)
        height_list_one[[j]] <- numeric(0)
      }else{
        peaks <- peaks[index_not_na]
        heights <- heights[index_not_na]

        # Remove peaks below AT
        index_above_at <- which(heights >= at_auto[mk])
        if(length(index_above_at) == 0){
          peak_list_one[[j]] <- character(0)
          height_list_one[[j]] <- numeric(0)
        }else{
          peaks <- peaks[index_above_at]
          heights <- heights[index_above_at]

          # Extract the genotype of one marker
          dt_gt_one_mk <- dt_gt_one[Marker == mk]
          gt <- as.character(dt_gt_one_mk[, index_gt_allele, with = FALSE])
          gt_one[j, ] <- gt

          # Make detectable peaks in the profile under the genotype
          cand_peak <- make_cand_peak(as.numeric(gt), c(sr_consider_auto[mk, ]))

          index_use <- which(is.element(peaks, as.character(cand_peak)))
          if(length(index_use) > 0){
            peak_list_one[[j]] <- peaks_adopt <- peaks[index_use]
            height_list_one[[j]] <- heights[index_use]
          }else{
            peak_list_one[[j]] <- character(0)
            height_list_one[[j]] <- numeric(0)
          }
        }
      }
    }

    names(peak_list_one) <- names(height_list_one) <- markers_auto
    peak_list[[i]] <- peak_list_one
    height_list[[i]] <- height_list_one
    gt_list[[i]] <- gt_one
  }

  names(peak_list) <- names(height_list) <- names(gt_list) <- sample_names
  return(list("peak_list" = peak_list,
              "height_list" = height_list,
              "gt_list" = gt_list))
}
