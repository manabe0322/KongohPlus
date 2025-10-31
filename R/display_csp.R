#' create_dt_display_csp
#'
#' @description The function to create the data.table for the displayed crime stain profile
#' @param dt_csp The data.table of the crime stain profile
#' @export
create_dt_display_csp <- function(dt_csp){
  sn_csp_all <- unique(dt_csp[, SampleName])
  n_rep <- length(sn_csp_all)
  markers <- unique(dt_csp[, Marker])
  n_mk <- length(markers)
  index_alleles <- grep("Allele", names(dt_csp))
  index_heights <- grep("Height", names(dt_csp))

  display_mat <- NULL

  for(i in 1:n_mk){
    mk <- markers[i]
    dt_csp_mk <- dt_csp[dt_csp[, Marker] == mk]
    dt_alleles <- dt_csp_mk[, index_alleles, with = FALSE]
    alleles_all <- unique(unlist(dt_alleles))
    alleles_all <- alleles_all[!is.na(alleles_all)]
    n_al_all <- length(alleles_all)

    options(warn = -1)
    alleles_all <- sort(alleles_all)
    alleles_all <- alleles_all[order(as.numeric(alleles_all))]
    options(warn = 0)

    if(n_al_all == 0){
      display_mk <- matrix(NA, 1, 2 + n_rep)
      display_mk[1, 1] <- mk
    }else{
      display_mk <- matrix(NA, n_al_all, 2)
      display_mk[1, 1] <- mk
      display_mk[, 2] <- alleles_all

      for(j in 1:n_rep){
        dt_csp_mk_s1 <- dt_csp_mk[dt_csp_mk[, SampleName] == sn_csp_all[j]]
        dt_alleles_s1 <- dt_csp_mk_s1[, index_alleles, with = FALSE]
        dt_heights_s1 <- dt_csp_mk_s1[, index_heights, with = FALSE]

        alleles_s1 <- unlist(dt_alleles_s1)
        heights_s1 <- unlist(dt_heights_s1)
        index_tmp_s1_1 <- which(!is.na(alleles_s1))
        alleles_s1 <- alleles_s1[index_tmp_s1_1]
        heights_s1 <- heights_s1[index_tmp_s1_1]

        index_tmp_s1_2 <- match(alleles_s1, alleles_all)
        display_heights_s1 <- rep(NA, n_al_all)
        display_heights_s1[index_tmp_s1_2] <- heights_s1

        display_mk <- cbind(display_mk, display_heights_s1)
      }
    }

    display_mat <- rbind(display_mat, display_mk)
  }

  # --- Create the data.table --- #

  dt_display_csp <- as.data.table(display_mat)
  names_display_csp <- c("Marker", "Allele", rep(NA, n_rep))
  names_display_csp[3:(n_rep + 2)] <- paste0("Height", 1:n_rep)
  names(dt_display_csp) <- names_display_csp

  # --- Change data type for each column --- #

  col_char <- c("Marker", "Allele")
  dt_display_csp[, (col_char) := lapply(.SD, as.character), .SDcols = col_char]
  col_numeric <- names_display_csp[grep("Height", names_display_csp)]
  dt_display_csp[, (col_numeric) := lapply(.SD, as.numeric), .SDcols = col_numeric]

  return(dt_display_csp)
}

#' define_rep_groups
#'
#' @description The function to define groups for replicated profiles
#' @param group_names A character vector for group names
define_rep_groups <- function(group_names){
  groups <- list()

  for(i in 1:length(group_names)){
    groups[[i]] <- paste0("Height", i)
  }

  names(groups) <- group_names

  return(groups)
}

#' create_sketch_csp
#'
#' @description The function to create sketch used for displaying csp
#' @param dt_display_csp The object created by the function 'create_dt_display_csp'
#' @param groups The object created by the function 'define_rep_groups'
create_sketch_csp <- function(dt_display_csp, groups){

  # --- Check if all group columns exist in dt_display_csp --- #

  grp_cols <- unlist(groups, use.names = FALSE)
  if(!all(grp_cols %in% names(dt_display_csp))) {
    missing <- setdiff(grp_cols, names(dt_display_csp))
    stop("The following columns are not found in the data: ", paste(missing, collapse = ", "))
  }

  top_id_th <- list()
  top_id_th[[1]] <- tags$th("Marker", rowspan = 2, style = "background:#f3f3f3; text-align:center; font-weight:bold;")
  top_id_th[[2]] <- tags$th("Allele", rowspan = 2, style = "background:#e0e0e0; text-align:center; font-weight:bold;")

  top_group_th <- list()
  for(i in 1:length(groups)){
    if(i %% 2 == 0){
      top_group_th[[i]] <- tags$th(names(groups)[i], colspan = 1, style = "background:#e0e0e0; text-align:center; font-weight:bold;")
    }else{
      top_group_th[[i]] <- tags$th(names(groups)[i], colspan = 1, style = "background:#f3f3f3; text-align:center; font-weight:bold;")
    }
  }

  second_th <- lapply(grp_cols, function(nm) tags$th(nm))

  sketch <- withTags(
    tags$table(
      class = "display",
      tags$thead(
        do.call(tags$tr, unname(c(top_id_th, top_group_th))),
        do.call(tags$tr, unname(second_th))
      )
    )
  )

  return(sketch)
}
