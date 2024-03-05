simulate_oph_eph <- function(pae_actual, tph, gt_s1, markers_auto, sr_consider_auto,
                             ae_alpha, phr_myu, phr_sigma, bsr_sigma, fsr_sigma, dsr_sigma, m2sr_myu, m2sr_sigma,
                             dt_alcor_bsr, dt_alcor_fsr, dt_alcor_dsr,
                             n_simu = 2){
  n_mk <- length(markers_auto)

  # Initial randomly generated OPH
  rand_oph_start <- 0.5 * tph * n_mk * rdirichlet(n = n_simu, alpha = tph * ae_alpha)
  rand_oph_auto <- vector(mode = "list", length = n_mk)
  names(rand_oph_auto) <- markers_auto

  # Initial EPH
  eph_start <- 0.5 * tph * n_mk * pae_actual
  eph_auto <- vector(mode = "list", length = n_mk)
  names(eph_auto) <- markers_auto
  for(i in 1:n_mk){
    # Determine peaks which are expected to be observed
    gt <- as.numeric(gt_s1[i, ])
    sr_consider <- sr_consider_auto[i, ]
    cand_peak <- make_cand_peak(gt, sr_consider)
    n_peak <- length(cand_peak)

    rand_oph_al_1 <- matrix(0, n_simu, n_peak)
    rand_oph_al_2 <- rand_oph_bs_1 <- rand_oph_bs_2 <- rand_oph_fs_1 <- rand_oph_fs_2 <- rand_oph_ds_1 <- rand_oph_ds_2 <- rand_oph_m2s_1 <- rand_oph_m2s_2 <- rand_oph_al_1

    # Randomly generate peak height ratios
    rand_oph_start_mk <- rand_oph_start[, i]
    rand_phr <- rlnorm(n = n_simu, meanlog = phr_myu[i], sdlog = sqrt(phr_sigma[i] / rand_oph_start_mk))
    rand_tap_al_1 <- rand_oph_start_mk * 2 / (1 + rand_phr)
    rand_tap_al_2 <- rand_oph_start_mk * 2 * rand_phr / (1 + rand_phr)

    # Extract average back stutter ratios
    mk <- markers_auto[i]
    dt_alcor_bsr_mk <- dt_alcor_bsr[mk]
    parent_bsr <- dt_alcor_bsr_mk[, Allele]
    average_bsr <- dt_alcor_bsr_mk[, Average_SR]
    index_bsr_1 <- match(gt[1], parent_bsr)
    if(!is.na(index_bsr_1)){
      average_bsr_al_1 <- average_bsr[index_bsr_1]
    }else{
      average_bsr_al_1 <- 0
    }
    index_bsr_2 <- match(gt[2], parent_bsr)
    if(!is.na(index_bsr_2)){
      average_bsr_al_2 <- average_bsr[index_bsr_2]
    }else{
      average_bsr_al_2 <- 0
    }

    # Extract average forward stutter ratios
    dt_alcor_fsr_mk <- dt_alcor_fsr[mk]
    parent_fsr <- dt_alcor_fsr_mk[, Allele]
    average_fsr <- dt_alcor_fsr_mk[, Average_SR]
    index_fsr_1 <- match(gt[1], parent_fsr)
    if(!is.na(index_fsr_1)){
      average_fsr_al_1 <- average_fsr[index_fsr_1]
    }else{
      average_fsr_al_1 <- 0
    }
    index_fsr_2 <- match(gt[2], parent_fsr)
    if(!is.na(index_fsr_2)){
      average_fsr_al_2 <- average_fsr[index_fsr_2]
    }else{
      average_fsr_al_2 <- 0
    }

    # Extract average double-back stutter ratios
    dt_alcor_dsr_mk <- dt_alcor_dsr[mk]
    parent_dsr <- dt_alcor_dsr_mk[, Allele]
    average_dsr <- dt_alcor_dsr_mk[, Average_SR]
    index_dsr_1 <- match(gt[1], parent_dsr)
    if(!is.na(index_dsr_1)){
      average_dsr_al_1 <- average_dsr[index_dsr_1]
    }else{
      average_dsr_al_1 <- 0
    }
    index_dsr_2 <- match(gt[2], parent_dsr)
    if(!is.na(index_dsr_2)){
      average_dsr_al_2 <- average_dsr[index_dsr_2]
    }else{
      average_dsr_al_2 <- 0
    }

    # Randomly generate stutter ratios
    rand_bsr_al_1 <- rlnorm(n = n_simu, meanlog = log(average_bsr_al_1), sdlog = sqrt(bsr_sigma[i] / rand_tap_al_1))
    rand_bsr_al_2 <- rlnorm(n = n_simu, meanlog = log(average_bsr_al_2), sdlog = sqrt(bsr_sigma[i] / rand_tap_al_2))
    rand_fsr_al_1 <- rlnorm(n = n_simu, meanlog = log(average_fsr_al_1), sdlog = sqrt(fsr_sigma[i] / rand_tap_al_1))
    rand_fsr_al_2 <- rlnorm(n = n_simu, meanlog = log(average_fsr_al_2), sdlog = sqrt(fsr_sigma[i] / rand_tap_al_2))
    rand_dsr_al_1 <- rlnorm(n = n_simu, meanlog = log(average_dsr_al_1), sdlog = sqrt(dsr_sigma[i] / rand_tap_al_1))
    rand_dsr_al_2 <- rlnorm(n = n_simu, meanlog = log(average_dsr_al_2), sdlog = sqrt(dsr_sigma[i] / rand_tap_al_2))
    if(sr_consider["M2S"]){
      rand_m2sr_al_1 <- rlnorm(n = n_simu, meanlog = log(m2sr_myu[i]), sdlog = sqrt(m2sr_sigma[i] / rand_tap_al_1))
      rand_m2sr_al_2 <- rlnorm(n = n_simu, meanlog = log(m2sr_myu[i]), sdlog = sqrt(m2sr_sigma[i] / rand_tap_al_2))
    }else{
      rand_m2sr_al_1 <- rand_m2sr_al_2 <- rep(0, n_simu)
    }

    # Calculate observed peak heights
    rand_tmp_al_1 <- 1 + rand_bsr_al_1 + rand_fsr_al_1 + rand_dsr_al_1 + rand_m2sr_al_1
    rand_tmp_al_2 <- 1 + rand_bsr_al_2 + rand_fsr_al_2 + rand_dsr_al_2 + rand_m2sr_al_2
    rand_oph_al_1[, match(gt[1], cand_peak)] <- rand_tap_al_1 / rand_tmp_al_1
    rand_oph_al_2[, match(gt[2], cand_peak)] <- rand_tap_al_2 / rand_tmp_al_2
    bs <- calc_bs(gt)
    rand_oph_bs_1[, match(bs[1], cand_peak)] <- rand_tap_al_1 * rand_bsr_al_1 / rand_tmp_al_1
    rand_oph_bs_2[, match(bs[2], cand_peak)] <- rand_tap_al_2 * rand_bsr_al_2 / rand_tmp_al_2
    fs <- calc_fs(gt)
    rand_oph_fs_1[, match(fs[1], cand_peak)] <- rand_tap_al_1 * rand_fsr_al_1 / rand_tmp_al_1
    rand_oph_fs_2[, match(fs[2], cand_peak)] <- rand_tap_al_2 * rand_fsr_al_2 / rand_tmp_al_2
    ds <- calc_ds(gt)
    rand_oph_ds_1[, match(ds[1], cand_peak)] <- rand_tap_al_1 * rand_dsr_al_1 / rand_tmp_al_1
    rand_oph_ds_2[, match(ds[2], cand_peak)] <- rand_tap_al_2 * rand_dsr_al_2 / rand_tmp_al_2
    if(sr_consider["M2S"]){
      m2s <- calc_m2s(gt)
      rand_oph_m2s_1[, match(m2s[1], cand_peak)] <- rand_tap_al_1 * rand_m2sr_al_1 / rand_tmp_al_1
      rand_oph_m2s_2[, match(m2s[2], cand_peak)] <- rand_tap_al_2 * rand_m2sr_al_2 / rand_tmp_al_2
    }
    rand_oph_auto_mk <- rand_oph_al_1 + rand_oph_al_2 + rand_oph_bs_1 + rand_oph_bs_2 + rand_oph_fs_1 + rand_oph_fs_2 + rand_oph_ds_1 + rand_oph_ds_2 + rand_oph_m2s_1 + rand_oph_m2s_2
    rand_oph_auto[[i]] <- as.vector(rand_oph_auto_mk)

    # Calculate expected peak heights
    eph_start_mk <- eph_start[i]
    eph_al_1 <- rep(0, n_peak)
    eph_al_2 <- eph_bs_1 <- eph_bs_2 <- eph_fs_1 <- eph_fs_2 <- eph_ds_1 <- eph_ds_2 <- eph_m2s_1 <- eph_m2s_2 <- eph_al_1
    tmp_al_1 <- 1 + average_bsr_al_1 + average_fsr_al_1 + average_dsr_al_1
    tmp_al_2 <- 1 + average_bsr_al_2 + average_fsr_al_2 + average_dsr_al_2
    if(sr_consider["M2S"]){
      tmp_al_1 <- tmp_al_1 + m2sr_myu[i]
      tmp_al_2 <- tmp_al_2 + m2sr_myu[i]
    }
    eph_al_1[match(gt[1], cand_peak)] <- eph_start_mk / tmp_al_1
    eph_al_2[match(gt[2], cand_peak)] <- eph_start_mk / tmp_al_2
    eph_bs_1[match(bs[1], cand_peak)] <- eph_start_mk * average_bsr_al_1 / tmp_al_1
    eph_bs_2[match(bs[2], cand_peak)] <- eph_start_mk * average_bsr_al_2 / tmp_al_2
    eph_fs_1[match(fs[1], cand_peak)] <- eph_start_mk * average_fsr_al_1 / tmp_al_1
    eph_fs_2[match(fs[2], cand_peak)] <- eph_start_mk * average_fsr_al_2 / tmp_al_2
    eph_ds_1[match(ds[1], cand_peak)] <- eph_start_mk * average_dsr_al_1 / tmp_al_1
    eph_ds_2[match(ds[2], cand_peak)] <- eph_start_mk * average_dsr_al_2 / tmp_al_2
    if(sr_consider["M2S"]){
      eph_m2s_1[match(m2s[1], cand_peak)] <- eph_start_mk * m2sr_myu[i] / tmp_al_1
      eph_m2s_2[match(m2s[2], cand_peak)] <- eph_start_mk * m2sr_myu[i] / tmp_al_2
    }
    eph_auto_mk <- eph_al_1 + eph_al_2 + eph_bs_1 + eph_bs_2 + eph_fs_1 + eph_fs_2 + eph_ds_1 + eph_ds_2 + eph_m2s_1 + eph_m2s_2
    eph_auto[[i]] <- as.vector(sapply(eph_auto_mk, rep, n_simu))
  }
  return(list("rand_oph_auto" = rand_oph_auto,
              "eph_auto" = eph_auto))
}

simulate_oph_eph_all_tph <- function(pae_data_use, tph_data_use, gt_use, markers_auto, sr_consider_auto,
                                     dt_result_ae, dt_result_phr, dt_result_bsr, dt_result_fsr, dt_result_dsr, dt_result_m2sr,
                                     dt_alcor_bsr, dt_alcor_fsr, dt_alcor_dsr,
                                     n_simu = 2){
  # Extract parameters
  ae_alpha <- dt_result_ae[, Alpha]
  phr_myu <- dt_result_phr[, Myu]
  phr_sigma <- dt_result_phr[, Sigma]
  bsr_sigma <- dt_result_bsr[, Sigma]
  fsr_sigma <- dt_result_fsr[, Sigma]
  dsr_sigma <- dt_result_dsr[, Sigma]
  m2sr_myu <- dt_result_m2sr[, Myu]
  m2sr_sigma <- dt_result_m2sr[, Sigma]
  n_mk <- length(markers_auto)

  # Prepare data.table of the allele repeat correction to extract data of each marker
  setkey(dt_alcor_bsr, Marker)
  setkey(dt_alcor_fsr, Marker)
  setkey(dt_alcor_dsr, Marker)

  # Define objects
  n_use <- length(tph_data_use)
  rand_oph_auto_all_tph <- vector(mode = "list", length = n_use)
  eph_auto_all_tph <- vector(mode = "list", length = n_use)

  for(i in 1:n_use){
    tmp <- simulate_oph_eph(pae_data_use[, i], tph_use[i], gt_use[[i]], markers_auto, sr_consider_auto,
                            ae_alpha, phr_myu, phr_sigma, bsr_sigma, fsr_sigma, dsr_sigma, m2sr_myu, m2sr_sigma,
                            dt_alcor_bsr, dt_alcor_fsr, dt_alcor_dsr,
                            n_simu)
    rand_oph_auto_all_tph[[i]] <- tmp$rand_oph_auto
    eph_auto_all_tph[[i]] <- tmp$eph_auto
  }
  return(list("rand_oph_auto_all_tph" = rand_oph_auto_all_tph,
              "eph_auto_all_tph" = eph_auto_all_tph))
}
