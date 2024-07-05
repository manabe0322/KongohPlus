#' create_dt_info_anameth
#'
#' @description The function to create the data.table for information on the analysis methods
#' @param path_pack The package path
create_dt_info_anameth <- function(path_pack){
  fn_anameth <- list.files(paste0(path_pack, "/extdata/anameth"))
  n_anameth <- length(fn_anameth)

  info_kit <- info_param <- info_at <- info_allele_prob <- rep("", n_anameth)
  for(i in 1:n_anameth){
    dt_anameth_1 <- fread(paste0(path_pack, "/extdata/anameth/", fn_anameth[i]))

    # Set information on the kit
    info_kit[i] <- dt_anameth_1$Value[dt_anameth_1$Category == "Kit"]

    # Set information on the parameters
    info_param[i] <- dt_anameth_1$Value[dt_anameth_1$Category == "Parameter"]

    # Set information on the AT
    ats <- as.numeric(dt_anameth_1$Value[dt_anameth_1$Category == "AT"])
    at_unique <- unique(ats)
    if(length(at_unique) == 1){
      info_at[i] <- as.character(at_unique)
    }else{
      info_at[i] <- paste0(min(at_unique), "-", max(at_unique))
    }

    # Set information on the allele probability
    info_allele_prob[i] <- dt_anameth_1$Value[dt_anameth_1$Subcategory == "Method"]
  }

  dt_info_anameth <- data.table(Name = gsub(".{4}$", "", fn_anameth),
                                Kit = info_kit,
                                Parameter = info_param,
                                AT = info_at,
                                Allele_Probability = info_allele_prob)
  return(dt_info_anameth)
}
