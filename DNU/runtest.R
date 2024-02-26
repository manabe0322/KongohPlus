path_pack <- "D:/RStudio_GitHub/KongohPlus"

dt_csp <- fread(paste0(path_pack, "/inst/extdata/examples/csp_example.csv"))
dt_csp <- change_data_type_csp(dt_csp)

dt_ref <- fread(paste0(path_pack, "/inst/extdata/examples/ref_example.csv"))
dt_ref <- change_data_type_ref(dt_ref)

dt_af <- fread(paste0(path_pack, "/inst/extdata/examples/Allele-count_GF_Japanese.csv"))
dt_af <- change_data_type_af(dt_af)

dt_kit <- fread(paste0(path_pack, "/inst/extdata/kit/GlobalFiler.csv"))

dt_at <- fread(paste0(path_pack, "/inst/extdata/at/at_GF_default.csv"))

dt_param <- fread(paste0(path_pack, "/inst/extdata/param/param_GF_default.csv"))

dt_alcor <- fread(paste0(path_pack, "/inst/extdata/alcor/alcor_GF_default.csv"))

dt_gtc_fltr <- fread(paste0(path_pack, "/inst/extdata/gtc_fltr/gtc_fltr_GF_default.csv"))

dt_gtc_prob <- fread(paste0(path_pack, "/inst/extdata/gtc_prob/gtc_prob_GF_default.csv"))

dt_anameth <- fread(paste0(path_pack, "/inst/extdata/anameth/GF_default.csv"))

hnc <- 2
warmup <- 1200
sampling <- 300
chains <- 4

stan_models <- list.files(paste0(path_pack, "/model"))
if(!is.element("model_main.rds", stan_models)){
  model_main <- cmdstan_model(paste0(path_pack, "/inst/stan/model_main.stan"))
  saveRDS(model_main, paste0(path_pack, "/inst/model/model_main.rds"))
}
