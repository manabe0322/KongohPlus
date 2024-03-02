#' compile_stan_code
#'
#' @description The function to compile stan codes
#' @param path_pack The package path
compile_stan_code <- function(path_pack){
  stan_models <- list.files(paste0(path_pack, "/model"))
  if(!is.element("model_ae.rds", stan_models)){
    model_ae <- cmdstan_model(paste0(path_pack, "/stan/model_ae.stan"))
    saveRDS(model_ae, paste0(path_pack, "/model/model_ae.rds"))
  }
  if(!is.element("model_sr_gen.rds", stan_models)){
    model_sr_gen <- cmdstan_model(paste0(path_pack, "/stan/model_sr_gen.stan"))
    saveRDS(model_sr_gen, paste0(path_pack, "/model/model_sr_gen.rds"))
  }
  if(!is.element("model_sr_together.rds", stan_models)){
    model_sr_together <- cmdstan_model(paste0(path_pack, "/stan/model_sr_together.stan"))
    saveRDS(model_sr_together, paste0(path_pack, "/model/model_sr_together.rds"))
  }
  if(!is.element("model_sr_uni.rds", stan_models)){
    model_sr_uni <- cmdstan_model(paste0(path_pack, "/stan/model_sr_uni.stan"))
    saveRDS(model_sr_uni, paste0(path_pack, "/model/model_sr_uni.rds"))
  }
  if(!is.element("model_phr.rds", stan_models)){
    model_phr <- cmdstan_model(paste0(path_pack, "/stan/model_phr.stan"))
    saveRDS(model_phr, paste0(path_pack, "/model/model_phr.rds"))
  }
  if(!is.element("model_gamma.rds", stan_models)){
    model_gamma <- cmdstan_model(paste0(path_pack, "/stan/model_gamma.stan"))
    saveRDS(model_gamma, paste0(path_pack, "/model/model_gamma.rds"))
  }
  if(!is.element("model_main.rds", stan_models)){
    model_main <- cmdstan_model(paste0(path_pack, "/stan/model_main.stan"))
    saveRDS(model_main, paste0(path_pack, "/model/model_main.rds"))
  }
}
