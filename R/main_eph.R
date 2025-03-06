

calc_eph <- function(gt_select_list, tph, mp, deg, ae, dt_average_sr, dt_size, mean_size){
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

  eph_tmp <- ae * tph * as.vector(sapply(mp, rep, 2)) * as.vector(sapply(deg, rep, 2))^((sizes - mean_size) / mean_size)

  bs_contrib <- calc_bs(al_contrib)
  fs_contrib <- calc_fs(al_contrib)
  ds_contrib <- calc_ds(al_contrib)
  m2s_contrib <- calc_m2s(al_contrib)
  peaks_eph <- sort(unique(c(al_contrib, bs_contrib, fs_contrib, ds_contrib, m2s_contrib)))
  ephs <- rep(0, length(peaks_eph))
  names(ephs) <- peaks_eph

  eph_al <- eph_tmp / (1 + bsrs + fsrs + dsrs + m2srs)
  eph_bs <- eph_tmp * bsrs / (1 + bsrs + fsrs + dsrs + m2srs)
  eph_fs <- eph_tmp * fsrs / (1 + bsrs + fsrs + dsrs + m2srs)
  eph_ds <- eph_tmp * dsrs / (1 + bsrs + fsrs + dsrs + m2srs)
  eph_m2s <- eph_tmp * m2srs / (1 + bsrs + fsrs + dsrs + m2srs)

  for(i in 1:length(eph_al)){
    ephs[peaks_eph == al_contrib[i]] <- ephs[peaks_eph == al_contrib[i]] + eph_al[i]
  }

  for(i in 1:length(eph_bs)){
    ephs[peaks_eph == bs_contrib[i]] <- ephs[peaks_eph == bs_contrib[i]] + eph_bs[i]
  }

  for(i in 1:length(eph_fs)){
    ephs[peaks_eph == fs_contrib[i]] <- ephs[peaks_eph == fs_contrib[i]] + eph_fs[i]
  }

  for(i in 1:length(eph_ds)){
    ephs[peaks_eph == ds_contrib[i]] <- ephs[peaks_eph == ds_contrib[i]] + eph_ds[i]
  }

  for(i in 1:length(eph_m2s)){
    ephs[peaks_eph == m2s_contrib[i]] <- ephs[peaks_eph == m2s_contrib[i]] + eph_m2s[i]
  }

  ephs <- ephs[ephs != 0]

  return(ephs)
}
