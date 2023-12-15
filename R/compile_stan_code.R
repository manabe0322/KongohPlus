#' compile_stan_code
#'
#' @description The function to compile stan codes
#' @param path_pack The package path
compile_stan_code <- function(path_pack){
  stan_models <- list.files(paste0(path_pack, "/data"))
  if(!is.element("model_ae.rds", stan_models)){
    model_ae <- cmdstan_model(paste0(path_pack, "/stan/model_ae.stan"))
    saveRDS(model_ae, paste0(path_pack, "/data/model_ae.rds"))
  }
  if(!is.element("model_sr_gen.rds", stan_models)){
    model_sr_gen <- cmdstan_model(paste0(path_pack, "/stan/model_sr_gen.stan"))
    saveRDS(model_sr_gen, paste0(path_pack, "/data/model_sr_gen.rds"))
  }
  if(!is.element("model_sr_together.rds", stan_models)){
    model_sr_together <- cmdstan_model(paste0(path_pack, "/stan/model_sr_together.stan"))
    saveRDS(model_sr_together, paste0(path_pack, "/data/model_sr_together.rds"))
  }
  if(!is.element("model_sr_uni.rds", stan_models)){
    model_sr_uni <- cmdstan_model(paste0(path_pack, "/stan/model_sr_uni.stan"))
    saveRDS(model_sr_uni, paste0(path_pack, "/data/model_sr_uni.rds"))
  }
}
