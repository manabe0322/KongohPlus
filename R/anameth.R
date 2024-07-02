#' create_dt_info_anameth
#'
#' @description The function to create the data.table for information on the analysis methods
create_dt_info_anameth <- function(path_pack){
  fn_anameth <- list.files(paste0(path_pack, "/extdata/anameth"))
  n_anameth <- length(fn_anameth)

  info_kit <- info_param <- info_at <- info_allele_prob <- info_dropin <- rep("", n_anameth)
  for(i in 1:n_anameth){
    dt_anameth_1 <- fread(paste0(path_pack, "/extdata/anameth/", fn_anameth[i]))

    # Set information on the kit
    info_kit[i] <- dt_anameth_1$Name[dt_anameth_1$Folder == "kit"]

    # Set information on the parameters
    info_param[i] <- dt_anameth_1$Name[dt_anameth_1$Folder == "param"]

    # Set information on the AT
    fn_at_1 <- paste0(dt_anameth_1$Name[dt_anameth_1$Folder == "at"], ".csv")
    dt_at_1 <- fread(paste0(path_pack, "/extdata/at/", fn_at_1))
    at_unique <- unique(dt_at_1[, AT])
    if(length(at_unique) == 1){
      info_at[i] <- as.character(at_unique)
    }else{
      info_at[i] <- paste0(min(at_unique), "-", max(at_unique))
    }

    # Set information on the allele probability
    fn_gtc_prob_1 <- paste0(dt_anameth_1$Name[dt_anameth_1$Folder == "gtc_prob"], ".csv")
    dt_gtc_prob_1 <- fread(paste0(path_pack, "/extdata/gtc_prob/", fn_gtc_prob_1))
    info_allele_prob[i] <- dt_gtc_prob_1$Value[dt_gtc_prob_1$Parameter == "Method"]

    # Set information on the drop-in

  }

  dt_info_anameth <- data.table(Name = gsub(".csv", "", fn_anameth),
                                Kit = info_kit,
                                Parameter = info_param,
                                AT = info_at,
                                Allele_Probability = info_allele_prob,
                                Dropin = info_dropin)
  return(dt_info_anameth)
}
