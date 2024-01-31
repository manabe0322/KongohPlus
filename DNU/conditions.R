path_pack <- "D:/RStudio_GitHub/KongohPlus"

dt_csp <- fread(paste0(path_pack, "/inst/extdata/examples/csp_example.csv"))
dt_csp <- change_data_type_csp(dt_csp)

dt_ref <- fread(paste0(path_pack, "/inst/extdata/examples/ref_example.csv"))
dt_ref <- change_data_type_ref(dt_ref)

dt_af <- fread(paste0(path_pack, "/inst/extdata/examples/Allele-count_GF_Japanese.csv"))
dt_af <- change_data_type_af(dt_af)

dt_parameter <- fread(paste0(path_pack, "/inst/extdata/parameter/param_GF_default.csv"))

dt_alcor <- fread(paste0(path_pack, "/inst/extdata/corrected_allele/alcor_GF_default.csv"))

dt_kit <- fread(paste0(path_pack, "/inst/extdata/kit/GlobalFiler.csv"))

dt_analysis_method <- fread(paste0(path_pack, "/inst/extdata/analysis_method/GF_default.csv"))
names <- dt_analysis_method[, Name]
values <- dt_analysis_method[, Value]
index_at <- grep("AT_", names)
dt_at <- data.table(Marker = gsub("AT_", "", names[index_at]),
                    AT = as.numeric(values[index_at]))
bsr_fltr <- as.numeric(values[names == "Filter of back stutter ratio"])
fsr_fltr <- as.numeric(values[names == "Filter of forward stutter ratio"])
dsr_fltr <- as.numeric(values[names == "Filter of double-back stutter ratio"])
m2sr_fltr <- as.numeric(values[names == "Filter of minus 2-nt stutter ratio"])
phr_fltr <- as.numeric(values[names == "Filter of peak height ratio"])
st <- as.numeric(values[names == "Stochastic threshold"])

hnc <- 2

stan_models <- list.files(paste0(path_pack, "/model"))
if(!is.element("model_main.rds", stan_models)){
  model_main <- cmdstan_model(paste0(path_pack, "/inst/stan/model_main.stan"))
  saveRDS(model_main, paste0(path_pack, "/inst/model/model_main.rds"))
}

fit_main$summary()
fit_main_summary <- fit_main$summary()
fit_main$draws() %>% mcmc_trace() -> traceplot

draws <- fit_main$draws()
sampling_data <- as_draws_matrix(draws)
#colnames(sampling_data)
median_mp <- matrix(0, hnc, 4)
for(i in 1:4){
  sampling_data_chain1 <- sampling_data[(1000 * (i - 1) + 1):(1000 * i), ]
  for(j in 1:hnc){
    median_mp[j, i] <- median(sampling_data_chain1[, 2 + j])
  }
}

sampling_data_mp <- matrix(0, 4000, 2)
sampling_data_mp[1:1000, 1] <- sampling_data[1:1000, 4]
sampling_data_mp[1:1000, 2] <- sampling_data[1:1000, 3]
sampling_data_mp[1001:2000, 1] <- sampling_data[1001:2000, 4]
sampling_data_mp[1001:2000, 2] <- sampling_data[1001:2000, 3]
sampling_data_mp[2001:3000, 1] <- sampling_data[2001:3000, 3]
sampling_data_mp[2001:3000, 2] <- sampling_data[2001:3000, 4]
sampling_data_mp[3001:4000, 1] <- sampling_data[3001:4000, 3]
sampling_data_mp[3001:4000, 2] <- sampling_data[3001:4000, 4]


sampling_data_deg <- matrix(0, 4000, 2)
sampling_data_deg[1:1000, 1] <- sampling_data[1:1000, 6]
sampling_data_deg[1:1000, 2] <- sampling_data[1:1000, 5]
sampling_data_deg[1001:2000, 1] <- sampling_data[1001:2000, 6]
sampling_data_deg[1001:2000, 2] <- sampling_data[1001:2000, 5]
sampling_data_deg[2001:3000, 1] <- sampling_data[2001:3000, 5]
sampling_data_deg[2001:3000, 2] <- sampling_data[2001:3000, 6]
sampling_data_deg[3001:4000, 1] <- sampling_data[3001:4000, 5]
sampling_data_deg[3001:4000, 2] <- sampling_data[3001:4000, 6]
