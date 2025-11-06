#' create_dt_info_anameth
#'
#' @description The function to create the data.table for information on the analysis methods
#' @param path_pack The package path
create_dt_info_anameth <- function(path_pack){
  fn_anameth <- list.files(paste0(path_pack, "/extdata/anameth"))
  n_anameth <- length(fn_anameth)

  if(n_anameth != 0){
    info_kit <- info_stutter <- info_param <- info_dropin <- rep("", n_anameth)
    for(i in 1:n_anameth){
      dt_anameth_1 <- fread(paste0(path_pack, "/extdata/anameth/", fn_anameth[i]))
      info_kit[i] <- dt_anameth_1$Name[dt_anameth_1$Category == "Kit"]
      info_stutter[i] <- dt_anameth_1$Name[dt_anameth_1$Category == "Stutter data"]
      info_param[i] <- dt_anameth_1$Name[dt_anameth_1$Category == "Parameter data"]
      info_dropin[i] <- dt_anameth_1$Name[dt_anameth_1$Category == "Drop-in data"]
    }

    dt_info_anameth <- data.table("Analysis method" = gsub(".{4}$", "", fn_anameth),
                                  "Kit" = info_kit,
                                  "Stutter data" = info_stutter,
                                  "Parameter data" = info_param,
                                  "Drop-in data" = info_dropin)
    return(dt_info_anameth)
  }else{
    return(NULL)
  }
}

#' create_dt_anameth_sr
#'
#' @description The function to create the data.table for stutter ratios as a part of the analysis method
#' @param list_dt_average_sr A list of the data.table for the average stutter ratios
#' @param markers_auto A vector of the autosomal markers
create_dt_anameth_sr <- function(list_dt_average_sr, markers_auto){
  vec_nrow_dt_average_sr <- sapply(list_dt_average_sr, nrow)
  nrow_sr_all <- sum(vec_nrow_dt_average_sr)
  sr_subcategory <- rep("", nrow_sr_all)
  sr_allele <- rep("", nrow_sr_all)
  sr_alcor_bsr <- rep("", nrow_sr_all)
  sr_bsr <- rep("", nrow_sr_all)
  sr_alcor_fsr <- rep("", nrow_sr_all)
  sr_fsr <- rep("", nrow_sr_all)
  sr_alcor_dsr <- rep("", nrow_sr_all)
  sr_dsr <- rep("", nrow_sr_all)
  sr_alcor_m2sr <- rep("", nrow_sr_all)
  sr_m2sr <- rep("", nrow_sr_all)

  index_sr_start <- 1
  for(i in 1:length(markers_auto)){
    mk <- markers_auto[i]
    dt_average_sr <- list_dt_average_sr[[mk]]

    index_assign <- index_sr_start:(index_sr_start + vec_nrow_dt_average_sr[i] - 1)
    sr_subcategory[index_assign] <- mk
    sr_allele[index_assign] <- dt_average_sr[["Allele"]]
    sr_alcor_bsr[index_assign] <- dt_average_sr[["Corrected allele (BSR)"]]
    sr_bsr[index_assign] <- dt_average_sr[["Average BSR"]]
    sr_alcor_fsr[index_assign] <- dt_average_sr[["Corrected allele (FSR)"]]
    sr_fsr[index_assign] <- dt_average_sr[["Average FSR"]]
    sr_alcor_dsr[index_assign] <- dt_average_sr[["Corrected allele (DSR)"]]
    sr_dsr[index_assign] <- dt_average_sr[["Average DSR"]]
    sr_alcor_m2sr[index_assign] <- dt_average_sr[["Corrected allele (M2SR)"]]
    sr_m2sr[index_assign] <- dt_average_sr[["Average M2SR"]]

    index_sr_start <- index_sr_start + vec_nrow_dt_average_sr[i]
  }

  dt_anameth_sr <- data.table(Category = rep("SR", nrow_sr_all),
                              Subcategory = sr_subcategory,
                              Name = rep("", nrow_sr_all),
                              Value = rep(NA, nrow_sr_all),
                              Allele = sr_allele,
                              "Corrected allele (BSR)" = sr_alcor_bsr,
                              "Average BSR" = sr_bsr,
                              "Corrected allele (FSR)" = sr_alcor_fsr,
                              "Average FSR" = sr_fsr,
                              "Corrected allele (DSR)" = sr_alcor_dsr,
                              "Average DSR" = sr_dsr,
                              "Corrected allele (M2SR)" = sr_alcor_m2sr,
                              "Average M2SR" = sr_m2sr)

  return(dt_anameth_sr)
}

#' create_default_dt_anameth
#'
#' @description The function to create the default data.table for a new analysis method
#' @param path_pack The package path of KongohPlus
#' @param kit_name A kit name
#' @param stutter_name The name of the stutter data
#' @param param_name The name of the parameter data
#' @param dropin_name The name of the drop-in data
create_default_dt_anameth <- function(path_pack, kit_name, stutter_name, param_name, dropin_name){
  dt_kit <- fread(paste0(path_pack, "/extdata/kit/", kit_name, ".csv"))
  markers_auto <- unique(dt_kit[["Marker"]][dt_kit[["Chromosome"]] == "Autosomal"])
  n_mk <- length(markers_auto)

  load(paste0(path_pack, "/extdata/stutter/", stutter_name, ".RData"))
  load(paste0(path_pack, "/extdata/param/", param_name, ".RData"))
  load(paste0(path_pack, "/extdata/dropin/", dropin_name, ".RData"))

  dt_expect_ae <- data_est_par$dt_expect_ae
  dt_est_gamma <- data_est_par$dt_est_gamma
  dt_est_dropin <- data_est_dropin$dt_est_dropin

  default_dt_anameth <- data.table(Category = c("Kit", "Stutter data", "Parameter data", "Drop-in data",
                                                rep("AT", n_mk),
                                                rep("Genotype filter", 6),
                                                rep("LSAE", n_mk),
                                                rep("Peak height", 2),
                                                rep("Drop-in", 3)),
                                   Subcategory = c(rep("", 4),
                                                   markers_auto,
                                                   "Maximum BSR", "Maximum FSR", "Maximum DSR", "Maximum M2SR", "Maximum number of drop-in", "Maximum height of drop-in",
                                                   dt_expect_ae$Marker,
                                                   "Scale parameter (allele)", "Scale parameter (stutter)",
                                                   "Frequency", "Shape parameter", "Scale parameter"),
                                   Name = c(kit_name, stutter_name, param_name, dropin_name, rep("", 2 * n_mk + 11)),
                                   Value = c(rep(NA, 4),
                                             rep(100, n_mk),
                                             0.3, 0.2, 0.05, 0.1, 2, 150,
                                             dt_expect_ae$LSAE,
                                             dt_est_gamma$Value,
                                             dt_est_dropin$Value[dt_est_dropin$Parameter == "Frequency"],
                                             dt_est_dropin$Value[dt_est_dropin$Parameter == "Shape parameter"],
                                             dt_est_dropin$Value[dt_est_dropin$Parameter == "Scale parameter"]
                                   ),
                                   Allele = rep(NA, 2 * n_mk + 15),
                                   "Corrected allele (BSR)" = rep(NA, 2 * n_mk + 15),
                                   "Average BSR" = rep(NA, 2 * n_mk + 15),
                                   "Corrected allele (FSR)" = rep(NA, 2 * n_mk + 15),
                                   "Average FSR" = rep(NA, 2 * n_mk + 15),
                                   "Corrected allele (DSR)" = rep(NA, 2 * n_mk + 15),
                                   "Average DSR" = rep(NA, 2 * n_mk + 15),
                                   "Corrected allele (M2SR)" = rep(NA, 2 * n_mk + 15),
                                   "Average M2SR" = rep(NA, 2 * n_mk + 15))

  default_dt_anameth_sr <- create_dt_anameth_sr(data_est_sr$list_dt_average_sr, markers_auto)

  default_dt_anameth <- rbind(default_dt_anameth, default_dt_anameth_sr)

  return(default_dt_anameth)
}

#' compare_list_dt_average_sr
#'
#' @description The function to compare 'list_dt_average_sr_origin' (before editing) with 'list_dt_average_sr' (after editing)
#' @param list_dt_average_sr_origin A list of the data.table for average stutter ratios in each marker (before editing)
#' @param list_dt_average_sr A list of the data.table for average stutter ratios in each marker (after editing)
compare_list_dt_average_sr <- function(list_dt_average_sr_origin, list_dt_average_sr){
  markers <- names(list_dt_average_sr_origin)
  bool_edited <- FALSE

  for(i in 1:length(markers)){
    mk <- markers[i]
    dt_average_sr_origin <- list_dt_average_sr_origin[[mk]]
    dt_average_sr <- list_dt_average_sr[[mk]]

    if(!setequal(dt_average_sr_origin$Allele, dt_average_sr$Allele)){
      bool_edited <- TRUE
      break
    }else if(any(dt_average_sr_origin[["Average BSR"]] != dt_average_sr[["Average BSR"]])){
      bool_edited <- TRUE
      break
    }else if(any(dt_average_sr_origin[["Average FSR"]] != dt_average_sr[["Average FSR"]])){
      bool_edited <- TRUE
      break
    }else if(any(dt_average_sr_origin[["Average DSR"]] != dt_average_sr[["Average DSR"]])){
      bool_edited <- TRUE
      break
    }else if(any(dt_average_sr_origin[["Average M2SR"]] != dt_average_sr[["Average M2SR"]])){
      bool_edited <- TRUE
      break
    }
  }

  return(bool_edited)
}

#' create_list_dt_average_sr
#'
#' @description The function to create the object 'list_dt_average_sr'
#' @param dt_anameth A data.table for the analysis method
#' @param markers_auto A vector of the autosomal markers
create_list_dt_average_sr <- function(dt_anameth, markers_auto){
  n_mk <- length(markers_auto)
  list_dt_average_sr <- vector(mode = "list", length = n_mk)
  names(list_dt_average_sr) <- markers_auto
  for(i in 1:n_mk){
    mk <- markers_auto[i]
    index_target <- sort(intersect(which(dt_anameth$Category == "SR"), which(dt_anameth$Subcategory == mk)))

    dt_average_sr <- data.table("Allele" = dt_anameth[["Allele"]][index_target],
                                "Corrected allele (BSR)" = dt_anameth[["Corrected allele (BSR)"]][index_target],
                                "Average BSR" = dt_anameth[["Average BSR"]][index_target],
                                "Corrected allele (FSR)" = dt_anameth[["Corrected allele (FSR)"]][index_target],
                                "Average FSR" = dt_anameth[["Average FSR"]][index_target],
                                "Corrected allele (DSR)" = dt_anameth[["Corrected allele (DSR)"]][index_target],
                                "Average DSR" = dt_anameth[["Average DSR"]][index_target],
                                "Corrected allele (M2SR)" = dt_anameth[["Corrected allele (M2SR)"]][index_target],
                                "Average M2SR" = dt_anameth[["Average M2SR"]][index_target])
    list_dt_average_sr[[i]] <- dt_average_sr
  }

  return(list_dt_average_sr)
}
