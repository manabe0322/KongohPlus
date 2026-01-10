#' change_var_al_plot
#'
#' @description The function to change the format of the variant alleles for plot
#' @param pk A peak number
#' @param rep_unit_len The length of the repeat unit
change_var_al_plot <- function(pk, rep_unit_len){
  return(round(pk %/% 1 + 10 * pk %% 1 / rep_unit_len, 3))
}

#' create_peaks_plot
#'
#' @description The function to create the format of peaks for plot
#' @param peaks A numeric vector of peaks
#' @param rep_unit_len The length of the repeat unit
#' @param q_distance The distance of the position Q from the maximum peak number
create_peaks_plot <- function(peaks, rep_unit_len, q_distance = 3){
  peaks_plot <- change_var_al_plot(peaks, rep_unit_len)
  pk_q <- max(setdiff(peaks_plot, c(99, 991, 992, 993, 994))) + q_distance

  count_q <- 0
  if(any(peaks_plot == 991)){
    peaks_plot[peaks_plot == 991] <- pk_q + count_q
    count_q <- count_q + 1
  }
  if(any(peaks_plot == 992)){
    peaks_plot[peaks_plot == 992] <- pk_q + count_q
    count_q <- count_q + 1
  }
  if(any(peaks_plot == 993)){
    peaks_plot[peaks_plot == 993] <- pk_q + count_q
    count_q <- count_q + 1
  }
  if(any(peaks_plot == 994)){
    peaks_plot[peaks_plot == 994] <- pk_q + count_q
    count_q <- count_q + 1
  }

  return(peaks_plot)
}

#' calc_eph
#'
#' @description The function to calculate the expected peak heights under the selected genotype combination
#' @param gt_select_list A list of genotypes in each contributor
#' @param temps The template amount of each contributor
#' @param deg The degradation parameters of each contributor
#' @param rep_ae_one The replication amplification efficiency in one replicate
#' @param ae The locus-specific amplification efficiencies in each marker
#' @param dt_average_sr The data.table for the average stutter ratios
#' @param dt_size The data.table for sizes
#' @param base_size The base size of the profile
calc_eph <- function(gt_select_list, temps, deg, rep_ae_one, ae, dt_average_sr, dt_size, base_size){
  # Extract SR
  al_contrib <- unlist(gt_select_list)
  index_sr <- match(al_contrib, dt_average_sr$Allele)
  bsrs <- dt_average_sr[["Average BSR"]][index_sr]
  bsrs[is.na(bsrs)] <- 0
  fsrs <- dt_average_sr[["Average FSR"]][index_sr]
  fsrs[is.na(fsrs)] <- 0
  dsrs <- dt_average_sr[["Average DSR"]][index_sr]
  dsrs[is.na(dsrs)] <- 0
  m2srs <- dt_average_sr[["Average M2SR"]][index_sr]
  m2srs[is.na(m2srs)] <- 0

  # Extract sizes
  index_size <- match(al_contrib, as.numeric(dt_size$Allele))
  sizes <- dt_size$Size[index_size]
  sizes[is.na(sizes)] <- median(dt_size$Size)

  eph_tmp <- ae * as.vector(sapply(temps * rep_ae_one, rep, 2)) * as.vector(sapply(deg, rep, 2))^((sizes - base_size) / 100)
  eph_al <- eph_tmp / (1 + bsrs + fsrs + dsrs + m2srs)
  eph_bs <- eph_tmp * bsrs / (1 + bsrs + fsrs + dsrs + m2srs)
  eph_fs <- eph_tmp * fsrs / (1 + bsrs + fsrs + dsrs + m2srs)
  eph_ds <- eph_tmp * dsrs / (1 + bsrs + fsrs + dsrs + m2srs)
  eph_m2s <- eph_tmp * m2srs / (1 + bsrs + fsrs + dsrs + m2srs)

  bs_contrib <- calc_bs(al_contrib)
  fs_contrib <- calc_fs(al_contrib)
  ds_contrib <- calc_ds(al_contrib)
  m2s_contrib <- calc_m2s(al_contrib)
  peaks_eph <- sort(unique(c(al_contrib, bs_contrib, fs_contrib, ds_contrib, m2s_contrib, 991, 992, 993, 994))) # 991, 992, 993, 994: temporary names of allele Q for each contributor
  ephs <- rep(0, length(peaks_eph))
  names(ephs) <- peaks_eph

  hnc <- length(al_contrib) / 2
  for(i in 1:hnc){
    for(j in 1:2){
      index <- 2 * (i - 1) + j

      if(al_contrib[index] == 99){
        if(ephs[peaks_eph == 990 + i] == 0){
          ephs[peaks_eph == 990 + i] <- eph_al[index]
        }
      }else{
        ephs[peaks_eph == al_contrib[index]] <- ephs[peaks_eph == al_contrib[index]] + eph_al[index]
        ephs[peaks_eph == bs_contrib[index]] <- ephs[peaks_eph == bs_contrib[index]] + eph_bs[index]
        ephs[peaks_eph == fs_contrib[index]] <- ephs[peaks_eph == fs_contrib[index]] + eph_fs[index]
        ephs[peaks_eph == ds_contrib[index]] <- ephs[peaks_eph == ds_contrib[index]] + eph_ds[index]
        ephs[peaks_eph == m2s_contrib[index]] <- ephs[peaks_eph == m2s_contrib[index]] + eph_m2s[index]
      }
    }
  }

  ephs <- ephs[ephs != 0]

  return(ephs)
}
